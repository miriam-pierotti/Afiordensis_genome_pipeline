#!/usr/bin/env python3
"""Step 08_1 (protein-coding genes): add eggNOG and InterProScan results to the
genome-projected TransDecoder GFF and make the attributes ENA-friendly.

Usage: 08_1_protein_coding.py genome.gff3 orf_type.gff3 eggnog.annotations
       interpro.tsv go-basic.obo genome.fasta out.gff3 report.tsv add_gene(true/false)

- genome.gff3: TransDecoder ORFs projected to genome coordinates (08_1 step 2).
  Its own UTRs are kept as they are.
- orf_type.gff3: transcript-level TransDecoder GFF, source of the ORF type
  (complete/5prime_partial/...); the projected GFF mislabels some ORFs with the
  type of another ORF of the same gene.
- eggNOG and InterProScan were run on the proteins, so the functional
  annotation (product, gene, EC_number, Dbxref) goes on the CDS only. gene
  and mRNA keep just ID/Parent; Ontology_term, em_* and TransDecoder 'Name'
  strings are dropped everywhere.
- Product name, in order of preference:
    1. eggNOG Description (first sentence), unless it is a GO term name (current
       GO release; exact synonyms only with GO-style wording, as many are protein
       names), GO-style wording for old terms, a sentence, uninformative or names
       an organism. Domain/family names become '<domain>-containing protein' /
       '<family> protein'; locus names in RENAME are replaced
    2. best (lowest e-value) informative Pfam hit from InterProScan, written
       the same way
    3. 'X family protein' from an eggNOG 'Belongs to the X family' sentence
    4. 'hypothetical protein'
  Text is tidied first: COG/Pfam/TIGRFAM prefixes, bracketed organism names,
  quotes and underscores removed, '+' signs lost by eggNOG restored.
  Every decision is logged in the report.
- Dbxref on the CDS: the InterPro and Pfam accessions of the protein in
  interpro.tsv (the INSDC db_xref names InterPro and PFAM).
- Marks partial CDS with partial=true and start_range/end_range (NCBI GFF3
  convention) on the CDS, mRNA and, when it shares the open end, the gene.
  A CDS is 5'/3' partial when TransDecoder says so (ORF type) or when it does
  not begin with ATG / end with a stop codon in the genome.
"""
import re
import sys
from collections import defaultdict
from urllib.parse import quote, unquote

(gff_in, orf_gff, ann_tsv, ipr_tsv, go_obo, genome_fa,
 gff_out, report_out, add_gene) = sys.argv[1:10]
add_gene = add_gene.lower() == "true"

SENTENCE = re.compile(
    r"^(belongs to|catalyzes|catalyses|involved in|plays|binds|required|may |"
    r"probably|possibly|component of|functions|acts |has |is |participates|"
    r"mediates|regulates|part of|this |in the |forms |transports|converts|"
    r"an? (?!disintegrin)|the |stores |promotes |assists |responsible |control of |major role |"
    r"core component |catalytic subunit of |electron carrier (?!protein$)|specifically |"
    r"essential |important |together with )|"
    r"\b(is|are|which|that|catalyzes|catalyses|accelerates?|prepares)\b",
    re.IGNORECASE)
UNINFORMATIVE = re.compile(
    r"(domain of unknown function|family of unknown function|\bDUF\d*|"
    r"protein of unknown function|uncharacteri[sz]ed|hypothetical|unknown function|"
    r"orthologous group|^(function|activity|protein|domain|family|receptor|complement)$)",
    re.IGNORECASE)
# GO-style names that are not in the current GO release (old/renamed terms)
GO_STYLE = re.compile(
    r"(\b(activity|binding|process)$)|^((positive|negative) )?regulation of|"
    r"^(cellular )?response to|^establishment of|^(activation|inhibition|modulation) (of|by)|"
    r"^it is involved in|\bactivity\b(?!-)|signal+ing pathway|\binvolved in\b|\bvia\b|"
    r"\bin response to\b|\bprocess\b|^transcription\b(?! factor)|templated transcription$|"
    r"\bunwinding$|\bprocessing$|\bsignal+ing by\b|\bDNA incision\b|cleavage, endonucleolytic|^recruitment of", re.IGNORECASE)
# a GO exact synonym only counts as a GO term when it has GO-style wording:
# many synonyms are protein names (cytochrome P450, thioredoxin, globin, ephrin)
GO_WORD = re.compile(
    r"\b(activity|binding|process|pathway|signal+ing|transport|transposition|"
    r"regulation|response|mediated|metabolism|metabolic|biosynthe\w*|catabol\w*|"
    r"assembly|secretion|homeostasis|checkpoint|transcription|import|export|"
    r"localization|organi[sz]ation|scrambling|complex|cleavage|constituent|"
    r"involved|development|differentiation)\b", re.IGNORECASE)
# organism names have no place in a product name; when in brackets they are
# removed, anywhere else the description is rejected
ORGANISM = (r"Drosophila|yeast|mouse|human|Arabidopsis|zebrafish|Xenopus|chicken|"
            r"S\. ?cerevisiae|S\. ?pombe|C\. ?elegans|E\. ?coli")
ORG_BRACKET = re.compile(r"\s*\((%s)\b[^)]*\)|,\s*(%s)\)" % (ORGANISM, ORGANISM))
ORG_ANY = re.compile(r"\b(%s)\b" % ORGANISM)
# eggNOG names that are clone/locus names, replaced with a descriptive name
RENAME = {"k02a2.6-like": "retrotransposon-derived protein",
          # 66-character word: longer than an EMBL line (59), so it gets split
          # mid-word and read back with a space; UniProt name of CDKAL1 instead
          "trna (n(6)-l-threonylcarbamoyladenosine(37)-c(2))-methylthiotransferase":
              "Threonylcarbamoyladenosine tRNA methylthiotransferase"}
GO_STYLE_SOFT = re.compile(
    r"\b(signaling|signalling|pathway|transport|development|differentiation|"
    r"morphogenesis|proliferation|migration|secretion|homeostasis|organization|"
    r"localization|biosynthesis|metabolism|catabolism|import|export|adhesion|"
    r"apoptosis|cascade)$", re.IGNORECASE)
PROTEIN_WORD = re.compile(
    r"\b(protein|factor|ligand|regulators?|suppressors?|receptor|enzyme|subunit)\b",
    re.IGNORECASE)
EC_OK = re.compile(r"^\d+\.(\d+|-)\.(\d+|-)\.(\d+|-|n\d+)$")
# attributes kept on each feature type (functional ones are added to CDS below)
KEEP = {"ID", "Parent"}
STOPS = {"TAA", "TAG", "TGA"}
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def enc(value):
    # GFF3 reserved characters (; = & , %) get percent-encoded; so does '+',
    # which some parsers read as a space
    return quote(value, safe=" :/()[]'-._|")


def parse_attrs(field):
    attrs = []
    for item in field.strip().strip(";").split(";"):
        if item:
            key, _, value = item.partition("=")
            attrs.append([key.strip(), value])
    return attrs


def get(attrs, key):
    for k, v in attrs:
        if k == key:
            return v
    return None


def read_go_names(path):
    names, synonyms = set(), set()
    in_term = False
    with open(path) as fh:
        for line in fh:
            if line.startswith("["):
                in_term = line.startswith("[Term]")
            elif in_term and line.startswith("name: "):
                names.add(line[6:].strip().lower())
            elif in_term and line.startswith("synonym: ") and " EXACT " in line:
                synonyms.add(line.split('"')[1].lower())
    return names, synonyms


def is_go(d):
    d = d.lower()
    return d in GO_NAMES or (d in GO_SYNONYMS and bool(GO_WORD.search(d)))


def tidy(d):
    """Clean up eggNOG/Pfam text: database prefixes, organism tags, quotes,
    '+' signs that eggNOG lost, underscores, repeated words."""
    d = re.sub(r"\s*Add an annotation\s*", " ", d)
    d = re.sub(r"^(COGs?\s+)?(K\s+)?(COG|KOG)\d+\s+|^Pfam:|^TIGRFAM\s+\S+\s+", "", d)
    d = ORG_BRACKET.sub(lambda m: ")" if m.group(0).startswith(",") else "", d)
    d = re.sub(r"\(\s?\)", "(+)", d)
    d = re.sub(r"\((\d) \)", r"(\1+)", d)
    d = re.sub(r"\s+-like\b", "-like", d)
    d = re.sub(r"(?<!\w)'([A-Za-z][^']*)'", r"\1", d)
    d = d.replace("_", "-")
    d = re.sub(r"\b([A-Za-z]{4,}) \1\b",
               lambda m: m.group(0) if m.group(1).lower() == "kinase" else m.group(1),
               d, flags=re.IGNORECASE)
    return re.sub(r"\s{2,}", " ", d).strip().rstrip(",").strip()


def as_protein(d, pfam_names=(), from_pfam=False):
    """Domain/family names -> '<domain>-containing protein' / '<family> protein'."""
    d = re.sub(r"\s*\(\w+ copies\)", "", d).strip().rstrip(".")
    # trailing abbreviation after a domain name: 'X domain (ABC)' -> 'X domain'
    d = re.sub(r"\b(domain|repeat|motif)(s?)\s+\([^()]*\)$", r"\1\2", d, flags=re.IGNORECASE)
    d = re.sub(r"\b(repeat|domain)s$", r"\1", d)
    if from_pfam:
        d = re.sub(r"\s+protein$", "", d)
    if re.search(r"zinc finger$", d, re.IGNORECASE):
        return d + " protein"
    if re.search(r"(domain|repeat|motif|finger|box|region)$", d, re.IGNORECASE):
        return d + "-containing protein"
    if re.search(r"\b[NC][- ]termin(al|us)$", d):
        return d + " domain-containing protein"
    if re.search(r"family$", d, re.IGNORECASE):
        return d + " protein"
    # a bare Pfam short name (e.g. CAP_GLY, TLD) listed for this protein by eggNOG
    if from_pfam or d.lower() in pfam_names:
        return d + " domain-containing protein"
    return d


def eggnog_product(desc, pfam_names):
    if not desc or desc == "-":
        return None, "no eggNOG description"
    d = unquote(desc).strip()
    # first sentence, not splitting at abbreviations such as 'C. elegans'
    d = re.split(r"(?<![\s(,][A-Z])(?<!^[A-Z])\. ", d)[0].rstrip(".").strip()
    if d.lower() in RENAME:
        return RENAME[d.lower()], "eggNOG Description (renamed)"
    if d.startswith("Pfam:"):
        pfam_names = set(pfam_names) | {d[5:].lower().replace("_", "-")}
    d = tidy(d)
    if is_go(d):
        return None, "eggNOG description is a GO term"
    if GO_STYLE.search(d) or (GO_STYLE_SOFT.search(d) and not PROTEIN_WORD.search(d)):
        return None, "eggNOG description is a GO-style term"
    if SENTENCE.search(d):
        return None, "eggNOG description is a sentence"
    if UNINFORMATIVE.search(d):
        return None, "eggNOG description uninformative"
    if ORG_ANY.search(d):
        return None, "eggNOG description names an organism"
    p = as_protein(d, pfam_names)
    if len(p) > 100:
        return None, "eggNOG description too long"
    return p, "eggNOG Description"


def pfam_product(desc):
    """Pfam family/domain description -> protein name."""
    d = tidy(re.sub(r"\s*\(\d+ copies\)", "", desc).strip().rstrip("."))
    if UNINFORMATIVE.search(d) or d.lower() in GO_NAMES:
        return None
    p = as_protein(d, from_pfam=True)
    return p if len(p) <= 100 else None


def family_product(desc):
    """'Belongs to the X family' -> 'X family protein' (used when nothing else fits)."""
    m = re.match(r"belongs to the (.+? (super)?family)\b", unquote(desc), re.IGNORECASE)
    if not m:
        return None
    p = tidy(m.group(1)) + " protein"
    return p if len(p) <= 100 and not UNINFORMATIVE.search(p) else None


def read_fasta(path):
    seqs, name, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name:
        seqs[name] = "".join(buf)
    return seqs


GO_NAMES, GO_SYNONYMS = read_go_names(go_obo)

# eggNOG annotations table
ann = {}
header = None
with open(ann_tsv) as fh:
    for line in fh:
        if line.startswith("##"):
            continue
        if line.startswith("#query"):
            header = line.lstrip("#").rstrip("\n").split("\t")
            continue
        if header is None or not line.strip():
            continue
        row = dict(zip(header, line.rstrip("\n").split("\t")))
        ann[row["query"]] = row

# Pfam hits from InterProScan: protein -> [(evalue, accession, description)],
# and the InterPro/Pfam accessions as db_xrefs: protein -> {"InterPro": [...], "PFAM": [...]}
pfam = defaultdict(list)
xrefs = defaultdict(lambda: {"InterPro": [], "PFAM": []})
with open(ipr_tsv) as fh:
    for line in fh:
        c = line.rstrip("\n").split("\t")
        if len(c) > 8 and c[3] == "Pfam":
            ev = float(c[8]) if c[8] not in ("-", "") else 1.0
            pfam[c[0]].append((ev, c[4], c[5]))
            for db, acc in (("InterPro", c[11] if len(c) > 11 else "-"), ("PFAM", c[4])):
                if acc not in ("-", "") and acc not in xrefs[c[0]][db]:
                    xrefs[c[0]][db].append(acc)

with open(gff_in) as fh:
    rows = [l.rstrip("\n").split("\t") for l in fh]
feats = [r for r in rows if len(r) == 9 and not r[0].startswith("#")]


# The TransDecoder projection writes features in a different order on each run
# (Perl hash order). Put them in a fixed order so the outputs are the same every
# run: genes by position, then each gene's mRNAs by ID, each followed by its parts.
def natural(s):
    return [int(t) if t.isdigit() else t for t in re.split(r"(\d+)", s)]


TYPE_ORDER = {"exon": 0, "CDS": 1, "five_prime_UTR": 2, "three_prime_UTR": 3}
genes, mrnas, parts = [], defaultdict(list), defaultdict(list)
for r in feats:
    attrs = parse_attrs(r[8])
    if r[2] == "gene":
        genes.append(r)
    elif r[2] == "mRNA":
        mrnas[get(attrs, "Parent")].append(r)
    else:
        parts[get(attrs, "Parent").split(",")[0]].append(r)
genes.sort(key=lambda r: (natural(r[0]), int(r[3]), int(r[4]), natural(get(parse_attrs(r[8]), "ID"))))
ordered = []
for g in genes:
    ordered.append(g)
    for m in sorted(mrnas.pop(get(parse_attrs(g[8]), "ID")), key=lambda r: natural(get(parse_attrs(r[8]), "ID"))):
        ordered.append(m)
        ordered += sorted(parts.pop(get(parse_attrs(m[8]), "ID")), key=lambda r: (int(r[3]), TYPE_ORDER[r[2]]))
if mrnas or parts or len(ordered) != len(feats):
    sys.exit("ERROR: features without a gene/mRNA parent in %s" % gff_in)
feats = ordered
rows = [["##gff-version 3"]] + ordered

# TransDecoder ORF type (complete/5prime_partial/...), from the transcript-level
# TransDecoder GFF (the genome projection mislabels some ORFs with the Name of
# another ORF of the same gene)
orf_type = {}
with open(orf_gff) as fh:
    for line in fh:
        c = line.rstrip("\n").split("\t")
        if len(c) == 9 and c[2] == "mRNA":
            m = re.search(r"ORF type:(\w+)", unquote(c[8]))
            orf_type[get(parse_attrs(c[8]), "ID")] = m.group(1) if m else None

# first pass: CDS structure and mRNA/gene coordinates
cds = defaultdict(list)
mrna_info, gene_span = {}, {}
for r in feats:
    attrs = parse_attrs(r[8])
    if r[2] == "CDS":
        cds[get(attrs, "Parent").split(",")[0]].append((int(r[3]), int(r[4]), r[7]))
    elif r[2] == "mRNA":
        mid = get(attrs, "ID")
        mrna_info[mid] = (r[0], int(r[3]), int(r[4]), r[6], get(attrs, "Parent"))
    elif r[2] == "gene":
        gene_span[get(attrs, "ID")] = (int(r[3]), int(r[4]))

# partial CDS
genome = read_fasta(genome_fa)
partial = {}          # mRNA ID -> (low end open, high end open) in genome coordinates
n5 = n3 = n_td_only = 0
for mid, parts in cds.items():
    seqid, _, _, strand, _ = mrna_info[mid]
    parts.sort()
    s = "".join(genome[seqid][a - 1:b] for a, b, _ in parts).upper()
    if strand == "-":
        s = s.translate(COMP)[::-1]
    first_phase = (parts[-1] if strand == "-" else parts[0])[2]
    t = orf_type.get(mid) or ""
    no_start = first_phase not in ("0", ".") or not s.startswith("ATG")
    no_stop = s[-3:] not in STOPS
    p5 = no_start or t in ("5prime_partial", "internal")
    p3 = no_stop or t in ("3prime_partial", "internal")
    if (p5 and not no_start) or (p3 and not no_stop):
        n_td_only += 1
    if p5 or p3:
        n5 += p5
        n3 += p3
        low, high = (p5, p3) if strand == "+" else (p3, p5)
        partial[mid] = (low, high, parts[0][0], parts[-1][1])

gene_partial = {}
for mid, (low, high, _, _) in partial.items():
    _, ms, me, _, gid = mrna_info[mid]
    gs, ge = gene_span[gid]
    gl, gh = gene_partial.get(gid, (False, False))
    gene_partial[gid] = (gl or (low and ms == gs), gh or (high and me == ge))


def partial_attrs(low, high, start, end):
    out = []
    if low or high:
        out.append(["partial", "true"])
    if low:
        out.append(["start_range", ".,%d" % start])
    if high:
        out.append(["end_range", "%d,." % end])
    return out


# product decisions per mRNA (= protein ID)
final = {}
report = []
dbx = {mid: ",".join("%s:%s" % (db, acc) for db in ("InterPro", "PFAM") for acc in xrefs[mid][db]) or None
       for mid in mrna_info}
for mid in mrna_info:
    row = ann.get(mid)
    desc = row.get("Description", "-") if row else "-"
    pref = row.get("Preferred_name", "-") if row else "-"
    pfam_names = {x.lower().replace("_", "-") for x in (row.get("PFAMs", "-") if row else "-").split(",")}
    product, reason = eggnog_product(desc, pfam_names) if row else (None, "no eggNOG hit")
    used_pfam = "-"
    if product is None:
        for ev, acc, pdesc in sorted(pfam.get(mid, [])):
            product = pfam_product(pdesc)
            if product:
                used_pfam = "%s %s" % (acc, pdesc)
                reason += "; Pfam " + acc
                break
    if product is None and row:
        product = family_product(desc)
        if product:
            reason += "; eggNOG 'Belongs to' family"
    ecs = [e for e in row.get("EC", "-").split(",") if EC_OK.match(e)] if row else []
    final[mid] = {
        "product": product or "hypothetical protein",
        "gene": pref if add_gene and pref and pref != "-" else None,
        "ec": ecs,
        "dbxref": dbx.get(mid),
    }
    report.append([mid, pref, unquote(desc), used_pfam, final[mid]["product"], reason,
                   orf_type.get(mid) or "-",
                   "5'" * (mid in partial and (partial[mid][0] if mrna_info[mid][3] == "+" else partial[mid][1]))
                   + "3'" * (mid in partial and (partial[mid][1] if mrna_info[mid][3] == "+" else partial[mid][0]))
                   or "-"])

# second pass: write the GFF
with open(gff_out, "w") as out:
    for cols in rows:
        if len(cols) != 9 or cols[0].startswith("#"):
            out.write("\t".join(cols) + "\n")
            continue
        ftype = cols[2]
        old = parse_attrs(cols[8])
        attrs = [kv for kv in old if kv[0] in KEEP]
        if ftype == "gene":
            gid = get(old, "ID")
            if gid in gene_partial:
                attrs += partial_attrs(*gene_partial[gid], int(cols[3]), int(cols[4]))
        elif ftype == "mRNA":
            mid = get(old, "ID")
            if mid in partial:
                attrs += partial_attrs(*partial[mid][:2], int(cols[3]), int(cols[4]))
        elif ftype == "CDS":
            mid = get(old, "Parent").split(",")[0]
            info = final[mid]
            attrs.append(["product", enc(info["product"])])
            if info["gene"]:
                attrs.append(["gene", enc(info["gene"])])
            if info["ec"]:
                attrs.append(["EC_number", ",".join(info["ec"])])
            if info["dbxref"]:
                attrs.append(["Dbxref", info["dbxref"]])
            if mid in partial:
                attrs += partial_attrs(*partial[mid])
        cols[8] = ";".join(f"{k}={v}" if v != "" else k for k, v in attrs)
        out.write("\t".join(cols) + "\n")

with open(report_out, "w") as fh:
    fh.write("mRNA_ID\tPreferred_name\teggNOG_Description\tPfam_used\tfinal_product\t"
             "reason\tTransDecoder_ORF_type\tpartial_end\n")
    for r in report:
        fh.write("\t".join(r) + "\n")

n = len(final)
print(f"mRNAs processed:                      {n}")
print(f"product from eggNOG Description:      {sum(1 for r in report if r[5].startswith('eggNOG Description'))}")
print(f"product from Pfam:                    {sum(1 for r in report if r[3] != '-')}")
print(f"product from eggNOG 'Belongs to':     {sum(1 for r in report if 'Belongs to' in r[5])}")
print(f"  (eggNOG description was a GO term:  {sum(1 for r in report if r[5].startswith('eggNOG description is a GO'))} in total)")
print(f"set to 'hypothetical protein':        {sum(1 for v in final.values() if v['product'] == 'hypothetical protein')}")
print(f"with gene name:                       {sum(1 for v in final.values() if v['gene'])}")
print(f"with EC number:                       {sum(1 for v in final.values() if v['ec'])}")
print(f"partial CDS (5' / 3' / any):          {n5} / {n3} / {len(partial)}")
print(f"  partial only by TransDecoder type:  {n_td_only} (codons look complete)")
print(f"eggNOG entries not matched to mRNA:   {len(set(ann) - set(final))}")
