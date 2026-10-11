#!/usr/bin/env python3
"""Step 08_3: post-fix an EMBLmyGFF3 flat file for what EMBLmyGFF3 cannot write.

Usage: 08_3_fix_embl.py in.embl annotation.gff3 out.embl

(exon and, for the mitogenome, mRNA features are already left out by
EMBLmyGFF3's translation_gff_feature_to_embl_feature.json, see 08_3_ena_embl.sh)

  1. /anticodon: EMBLmyGFF3 quotes it and wraps it mid-word ("seq:cg\ng");
     INSDC wants it unquoted. Rewritten unquoted, wrapped after a comma.
  2. Identical features (same key, location and locus_tag) are kept once: two
     TransDecoder ORFs on one transcript give two identical mRNA features. The
     mRNA kept is marked partial at an end if any of the ORFs is. The UTR
     features of such transcripts are dropped, since each ORF's UTRs cover the
     other ORF's CDS once the two mRNAs are one.
  3. Loci whose CDS have different eggNOG /gene names ("Features sharing
     locus_tag are associated with gene qualifiers with different values") lose
     /gene on their CDS: the isoforms are often different genes (e.g. AIFM3 and
     PPM1D), so one name for all would mislabel some. /product stays.
  4. Partial ends from the GFF (start_range / end_range on gene, mRNA, CDS and
     tRNA, set by 08_1 from codons, TransDecoder ORF type and tRNAscan-SE
     truncation) become < / > in the
     location. EMBLmyGFF3 only marks CDS without start/stop codon, so 5'-partial
     ORFs that happen to begin with ATG looked complete.
  5. Circular entries: a gene across the origin (mitochondrial nad2) comes out
     as gene 16676..17710 and CDS join(1..81,16676..>17710); both become
     join(16676..17710,1..81).
Stops if a GFF attribute has a word longer than an EMBL line (see read_partials).
"""
import re
import sys
from collections import Counter, defaultdict
from urllib.parse import unquote

FT = "FT                   "  # qualifier indent (21 chars)


def read_features(lines):
    """Split the FT lines of one entry into features (lists of lines)."""
    feats = []
    for line in lines:
        if line[5] != " ":  # 'FT   key ...' starts a new feature
            feats.append([line])
        else:
            feats[-1].append(line)
    return feats


def parse(feat):
    """Return key, location string, and qualifiers as [name, value, lines]."""
    key = feat[0][5:21].strip()
    loc = feat[0][21:].rstrip("\n")
    quals = []
    i = 1
    while i < len(feat) and not feat[i][21:].startswith("/"):
        loc += feat[i][21:].rstrip("\n")  # location wrapped over several lines
        i += 1
    open_quote = False
    for line in feat[i:]:
        text = line[21:].rstrip("\n")
        if not open_quote and text.startswith("/"):
            name, _, value = text[1:].partition("=")
            quals.append([name, value, [line]])
        else:
            quals[-1][1] += (" " if open_quote else "") + text
            quals[-1][2].append(line)
        open_quote = (quals[-1][1].count('"') % 2) == 1
    return key, loc, quals


def gff_id(quals):
    """The GFF ID that EMBLmyGFF3 keeps as /note="ID:..."."""
    for name, value, _ in quals:
        if name == "note" and value.startswith('"ID:'):
            return value[4:-1]
    return None


def wrap(first, text):
    """Lines of a location or qualifier, wrapped after commas at 80 columns."""
    out = []
    prefix = first
    while len(text) > 59:
        cut = text.rfind(",", 0, 59) + 1
        out.append(prefix + text[:cut] + "\n")
        text = text[cut:]
        prefix = FT
    out.append(prefix + text + "\n")
    return out


def mark_partial(loc, low, high):
    """Add < to the lowest and > to the highest coordinate of a location."""
    nums = list(re.finditer(r"[<>]?\d+", loc))
    first, last = nums[0], nums[-1]
    if high and not last.group().startswith(">"):
        loc = loc[:last.start()] + ">" + loc[last.start():]
    if low and not first.group().startswith("<"):
        loc = loc[:first.start()] + "<" + loc[first.start():]
    return loc


def read_partials(gff):
    """GFF ID -> (low end partial, high end partial) from start_range/end_range.

    Also stops if a product/gene/note has a word longer than an EMBL line: it
    would be split mid-word and read back with a space in it."""
    partial = {}
    with open(gff) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            for text in re.findall(r"(?:^|;)(?:product|gene|Note)=([^;]*)", line.rstrip("\n").split("\t")[-1]):
                for word in unquote(text).split():
                    if len(word) <= 55:
                        continue
                    sys.exit("ERROR: word longer than an EMBL line in GFF, rename it in 08_1: " + word)
            col = line.rstrip("\n").split("\t")
            if len(col) < 9 or col[2] not in ("gene", "mRNA", "CDS", "tRNA"):
                continue
            attr = dict(a.split("=", 1) for a in col[8].split(";") if "=" in a)
            low, high = "start_range" in attr, "end_range" in attr
            if low or high:
                old = partial.get(attr["ID"], (False, False))
                partial[attr["ID"]] = (old[0] or low, old[1] or high)
    return partial


def across_origin(loc, length):
    """join(1..x,...,y..length) of a + strand feature -> join(y..length,...,1..x)."""
    parts = re.findall(r"<?(\d+)\.\.>?(\d+)", loc)
    if loc.startswith("complement") or len(parts) < 2:
        return None
    if parts[0][0] != "1" or parts[-1][1] != str(length):
        return None
    return "join(%s)" % ",".join("%s..%s" % p for p in parts[1:] + parts[:1])


def entries(path):
    """Yield (header lines, FT lines, trailer lines) for each entry."""
    head, ft, tail = [], [], []
    with open(path) as fh:
        for line in fh:
            if line.startswith("FT"):
                ft.append(line)
            elif ft:
                tail.append(line)
            else:
                head.append(line)
            if line.startswith("//"):
                yield head, ft, tail
                head, ft, tail = [], [], []


def main(inp, gff, out):
    partial = read_partials(gff)

    # pass 1: CDS /gene names per locus_tag, and transcripts carrying 2+ ORFs
    names = defaultdict(set)
    multi_orf = set()
    for _, ft, _ in entries(inp):
        mrnas = defaultdict(list)
        for feat in read_features(ft):
            key, loc, quals = parse(feat)
            q = {n: v for n, v, _ in quals}
            if key == "CDS" and "gene" in q:
                names[q["locus_tag"]].add(q["gene"])
            if key == "mRNA":
                mrnas[(loc, q.get("locus_tag"))].append(gff_id(quals))
        for ids in mrnas.values():
            if len(ids) > 1:
                multi_orf.update(ids)
                # the one mRNA kept is open at an end if any of its ORFs is
                ends = [partial.get(i, (False, False)) for i in ids]
                for i in ids:
                    partial[i] = (any(e[0] for e in ends), any(e[1] for e in ends))
    conflict = {tag for tag, n in names.items() if len(n) > 1}

    # pass 2: rewrite
    stats = Counter()
    with open(out, "w") as oh:
        for head, ft, tail in entries(inp):
            oh.writelines(head)
            seen = set()
            circular = "; circular;" in head[0]
            length = int(re.search(r"(\d+) BP\.", head[0]).group(1))
            wrapped = {}     # locus_tag -> joined CDS location, for the gene
            if circular:
                for feat in read_features(ft):
                    key, loc, quals = parse(feat)
                    new = across_origin(loc, length) if key == "CDS" else None
                    if new:
                        wrapped[{n: v for n, v, _ in quals}["locus_tag"]] = new
            for feat in read_features(ft):
                key, loc, quals = parse(feat)
                q = {n: v for n, v, _ in quals}
                fid = gff_id(quals)
                if key in ("5'UTR", "3'UTR") and fid and fid.rsplit(".utr", 1)[0] in multi_orf:
                    stats["UTR of 2-ORF transcript dropped"] += 1
                    continue
                if key != "source":
                    ident = (key, loc, q.get("locus_tag"))
                    if ident in seen:
                        stats["duplicate %s dropped" % key] += 1
                        continue
                    seen.add(ident)
                if key in ("gene", "mRNA", "CDS", "tRNA") and q.get("locus_tag") in wrapped:
                    cds_loc = wrapped[q["locus_tag"]]
                    if key == "gene":
                        first = re.match(r"(\d+)\.\.", loc).group(1)
                        last = re.findall(r"(\d+)\)$", cds_loc)[0]
                        loc = "join(%s..%d,1..%s)" % (first, length, last)
                    else:
                        loc = cds_loc
                    stats["%s joined across the origin" % key] += 1
                    lines = wrap("FT   %-16s" % key, loc)
                elif key in ("gene", "mRNA", "CDS", "tRNA") and fid in partial:
                    new = mark_partial(loc, *partial[fid])
                    if new != loc:
                        stats["%s partial end added" % key] += 1
                        loc = new
                    lines = wrap("FT   %-16s" % key, loc)
                else:
                    lines = [feat[0]] + feat[1:len(feat) - sum(len(l) for _, _, l in quals)]
                for name, value, qlines in quals:
                    if name == "anticodon":
                        lines += wrap(FT, "/anticodon=" + value.replace('"', "").replace(" ", ""))
                        stats["anticodon unquoted"] += 1
                    elif key == "CDS" and name == "gene" and q["locus_tag"] in conflict:
                        stats["CDS /gene removed (conflicting locus)"] += 1
                    else:
                        lines += qlines
                oh.writelines(lines)
            oh.writelines(tail)
    print("loci with conflicting CDS /gene: %d" % len(conflict))
    for k, v in sorted(stats.items()):
        print("%-40s %d" % (k, v))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])
