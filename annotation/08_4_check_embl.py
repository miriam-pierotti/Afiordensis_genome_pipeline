#!/usr/bin/env python3
"""Step 08_4: independent check of the ENA flat file against its inputs.

Usage: 08_4_check_embl.py flatfile.embl.gz nuclear.gff3 mito.gff3 genome.fasta mito.fasta
       transdecoder.pep eggnog.annotations interpro.tsv trna_hc.out

Parses the flat file with Biopython (not with the code that wrote it) and checks:
  sequences   every entry equals a FASTA sequence; gap features = N runs
  proteins    every CDS translates (table 1 nuclear, 4 mito) without internal
              stops; a CDS without start/stop codon is marked partial
  GFF         every gene/mRNA/CDS/UTR/rRNA/tRNA of the GFFs is in the EMBL with
              the same location, strand, partial ends, product, gene, EC_number,
              db_xref; GFF features missing from the EMBL are only the expected
              ones (identical isoform copies, the 2nd mRNA of 2-ORF transcripts
              and their UTRs)
  sources     nuclear proteins = TransDecoder proteins; /gene = eggNOG
              Preferred_name (or none); db_xref = InterProScan accessions;
              tRNAs = tRNAscan-SE high-confidence set
  format      unique locus_tag per gene, features inside their gene, no
              URL-encoded text, no double spaces, valid EC numbers and db_xrefs
Prints one line per check and exits with status 1 if any check fails.
"""
import gzip
import re
import sys
from collections import Counter, defaultdict
from urllib.parse import unquote

from Bio import SeqIO

(flatfile, nuclear_gff, mito_gff, genome_fa, mito_fa, pep_fa,
 eggnog_ann, interpro_tsv, trna_hc_out) = sys.argv[1:10]

fails = []


def check(name, problems, n):
    """Print one check line; problems is a list of examples."""
    status = "OK  " if not problems else "FAIL"
    print("%s %-55s %d checked, %d problems %s" % (status, name, n, len(problems), problems[:3] or ""))
    if problems:
        fails.append(name)


# ---------------- inputs
fasta = {}
for path in (genome_fa, mito_fa):
    for r in SeqIO.parse(path, "fasta"):
        fasta[r.id] = str(r.seq).upper()
mito_name = next(SeqIO.parse(mito_fa, "fasta")).id

gff, parent = {}, {}
for path in (nuclear_gff, mito_gff):
    for line in open(path):
        c = line.rstrip("\n").split("\t")
        if line.startswith("#") or len(c) != 9:
            continue
        a = {k: unquote(v) for k, v in (x.split("=", 1) for x in c[8].split(";") if "=" in x)}
        d = gff.setdefault(a["ID"], {"type": c[2], "seq": c[0], "strand": c[6], "parts": [], "attr": a})
        d["parts"].append((int(c[3]), int(c[4])))
        if "Parent" in a:
            parent[a["ID"]] = a["Parent"]
# mRNA/tRNA/rRNA locations are their exons; UTR pieces of one transcript are one feature
exons = defaultdict(list)
for i, g in gff.items():
    if g["type"] == "exon":
        exons[parent[i]] += g["parts"]
utr_group = defaultdict(list)
for i, g in gff.items():
    if g["type"] in ("mRNA", "tRNA", "rRNA") and exons[i]:
        g["parts"] = exons[i]
    if g["type"] in ("five_prime_UTR", "three_prime_UTR"):
        utr_group[(parent[i], g["type"])].append(i)
for ids in utr_group.values():
    allp = [p for i in ids for p in gff[i]["parts"]]
    for i in ids:
        gff[i]["parts"] = allp
        gff[i]["group"] = ids


def top(i):
    while i in parent:
        i = parent[i]
    return i


# mRNAs with the same location in one gene (2 ORFs on 1 transcript) become one EMBL
# mRNA, partial at an end if any of them is
same_mrna = defaultdict(list)
for i, g in gff.items():
    if g["type"] == "mRNA":
        same_mrna[(top(i), tuple(sorted(g["parts"])))].append(i)
two_orf = {i for v in same_mrna.values() if len(v) > 1 for i in v}
for ids in same_mrna.values():
    if len(ids) > 1:
        lo_any = any("start_range" in gff[i]["attr"] for i in ids)
        hi_any = any("end_range" in gff[i]["attr"] for i in ids)
        for i in ids:
            gff[i]["ends"] = (lo_any, hi_any)

pep = {r.id: str(r.seq).rstrip("*") for r in SeqIO.parse(pep_fa, "fasta")}
egg = {}
for line in open(eggnog_ann):
    if not line.startswith("#"):
        c = line.rstrip("\n").split("\t")
        egg[c[0]] = c
ipr = defaultdict(set)
for line in open(interpro_tsv):
    c = line.rstrip("\n").split("\t")
    ipr[c[0]].add("PFAM:" + c[4])
    if len(c) > 11 and c[11] not in ("-", ""):
        ipr[c[0]].add("InterPro:" + c[11])
hc = set()
for line in open(trna_hc_out).readlines()[3:]:
    c = [x.strip() for x in line.split("\t")]
    if "high confidence set" in c[-1]:
        hc.add("%s.trna%s" % (c[0], c[1]))

# ---------------- flat file
P = defaultdict(list)    # check name -> problems
N = Counter()            # check name -> items checked
KEY = {"five_prime_UTR": "5'UTR", "three_prime_UTR": "3'UTR"}
STARTS = {1: {"ATG"}, 4: {"ATG", "ATA", "ATT", "ATC", "TTG", "TTA", "CTG", "GTG"}}
seen_ids, tag_gene, gene_tag = set(), {}, {}
entries = 0
for rec in SeqIO.parse(gzip.open(flatfile, "rt"), "embl"):
    entries += 1
    s = str(rec.seq).upper()
    name = [k for k, v in fasta.items() if v == s]
    N["sequences: entry = FASTA sequence"] += 1
    if not name:
        P["sequences: entry = FASTA sequence"].append(len(s))
        continue
    name = name[0]
    is_mt = name == mito_name
    topology = rec.annotations.get("topology")
    N["sequences: topology (mito circular, rest linear)"] += 1
    if topology != ("circular" if is_mt else "linear"):
        P["sequences: topology (mito circular, rest linear)"].append((name, topology))
    gaps = [(int(f.location.start) + 1, int(f.location.end)) for f in rec.features if f.type == "gap"]
    N["sequences: gap features = N runs"] += 1
    if gaps != [(m.start() + 1, m.end()) for m in re.finditer("N+", s)]:
        P["sequences: gap features = N runs"].append(name)

    by_tag = defaultdict(list)
    for f in rec.features:
        q = f.qualifiers
        if f.type in ("source", "gap"):
            continue
        for k, vs in q.items():
            for v in vs:
                N["format: no URL-encoded text, no double spaces"] += 1
                if re.search(r"%[0-9A-Fa-f]{2}", v) or "  " in v:
                    P["format: no URL-encoded text, no double spaces"].append((k, v))
        for e in q.get("EC_number", []):
            N["format: EC numbers valid"] += 1
            if not re.fullmatch(r"\d+\.(\d+|-)\.(\d+|-)\.(n?\d+|-)", e):
                P["format: EC numbers valid"].append(e)
        for x in q.get("db_xref", []):
            N["format: db_xref valid"] += 1
            if not re.fullmatch(r"PFAM:PF\d{5}|InterPro:IPR\d{6}", x):
                P["format: db_xref valid"].append(x)
        tag = q.get("locus_tag", [None])[0]
        by_tag[tag].append(f)
        lo = type(f.location.start).__name__ == "BeforePosition"
        hi = type(f.location.end).__name__ == "AfterPosition"

        # proteins
        if f.type == "CDS":
            tt = int(q.get("transl_table", ["0"])[0])
            N["proteins: transl_table (1 nuclear, 4 mito)"] += 1
            if tt != (4 if is_mt else 1):
                P["proteins: transl_table (1 nuclear, 4 mito)"].append((name, tt))
                continue
            nt = f.location.extract(rec.seq)[int(q["codon_start"][0]) - 1:]
            nt = nt[:len(nt) // 3 * 3]
            prot = str(nt.translate(table=tt))
            p5, p3 = (lo, hi) if f.location.strand != -1 else (hi, lo)
            N["proteins: no internal stop"] += 1
            if "*" in prot.rstrip("*")[:-1] or "*" in prot[:-1]:
                P["proteins: no internal stop"].append(tag)
            N["proteins: no start/stop codon -> marked partial"] += 1
            if (not p5 and str(nt[:3]).upper() not in STARTS[tt]) or (not p3 and not prot.endswith("*")):
                P["proteins: no start/stop codon -> marked partial"].append(tag)

        # GFF comparison
        ids = [n[3:] for n in q.get("note", []) if n.startswith("ID:")]
        fid = ids[0] if ids else ("%s-%s" % ({"gene": "gene", "CDS": "cds", "tRNA": "trna", "rRNA": "rrna",
                                               "intron": "intron"}.get(f.type, f.type), q["gene"][0])
                                  if is_mt and "gene" in q else None)
        if is_mt and f.type == "intron":
            fid = [i for i, g in gff.items() if g["type"] == "intron" and g["seq"] == name
                   and g["parts"] == [(int(f.location.start) + 1, int(f.location.end))]]
            fid = fid[0] if fid else None
        g = gff.get(fid)
        N["GFF: EMBL feature found in GFF"] += 1
        if g is None:
            P["GFF: EMBL feature found in GFF"].append((name, f.type, fid))
            continue
        seen_ids.update(g.get("group", [fid]))
        eparts = [(int(p.start) + 1, int(p.end)) for p in f.location.parts]
        gparts = g["parts"]
        if is_mt and fid.startswith("gene-") and len(eparts) > 1:     # gene across the origin
            gparts = [p for i, x in gff.items() if x["type"] == "CDS" and top(i) == fid for p in x["parts"]]
        N["GFF: same location, strand, partial ends"] += 1
        same_loc = sorted(eparts) == sorted(gparts) and g["seq"] == name and \
            f.location.strand == (1 if g["strand"] == "+" else -1)
        same_partial = f.type not in ("gene", "mRNA", "CDS", "tRNA", "rRNA") or \
            (lo, hi) == g.get("ends", ("start_range" in g["attr"], "end_range" in g["attr"]))
        if not (same_loc and same_partial) or KEY.get(g["type"], g["type"]) != f.type:
            P["GFF: same location, strand, partial ends"].append((fid, eparts[:2], sorted(gparts)[:2]))
        for qual, attr in (("product", "product"), ("EC_number", "EC_number"),
                           ("db_xref", "Dbxref"), ("anticodon", "anticodon"), ("note", "note")):
            want = set(g["attr"][attr].split(",")) if attr in ("EC_number", "Dbxref") and attr in g["attr"] \
                else ({g["attr"][attr]} if attr in g["attr"] else set())
            have = set(q.get(qual, [])) if qual != "note" else {n for n in q.get("note", []) if not n.startswith(("ID:", "source:"))}
            if qual == "anticodon":
                have = {v.replace(" ", "") for v in have}
                want = {v.replace(" ", "") for v in want}
            N["GFF: same product, EC, db_xref, anticodon, note"] += 1
            if have != want:
                P["GFF: same product, EC, db_xref, anticodon, note"].append((fid, qual, have, want))
        if f.type == "CDS":
            N["GFF: same /gene (or dropped in a conflicting locus)"] += 1
            if "gene" in q and q["gene"][0] != g["attr"].get("gene"):
                P["GFF: same /gene (or dropped in a conflicting locus)"].append(fid)

        # sources
        if f.type == "CDS" and not is_mt:
            pid = fid[4:]
            N["sources: protein = TransDecoder protein"] += 1
            if prot.rstrip("*") != pep.get(pid):
                P["sources: protein = TransDecoder protein"].append(pid)
            e = egg.get(pid)
            N["sources: /gene = eggNOG Preferred_name"] += 1
            if "gene" in q and (e is None or q["gene"][0] != e[8]):
                P["sources: /gene = eggNOG Preferred_name"].append(pid)
            N["sources: db_xref = InterProScan accessions"] += 1
            if set(q.get("db_xref", [])) != ipr.get(pid, set()):
                P["sources: db_xref = InterProScan accessions"].append(pid)
        if f.type == "tRNA" and not is_mt:
            N["sources: tRNA in tRNAscan-SE high-confidence set"] += 1
            if fid not in hc:
                P["sources: tRNA in tRNAscan-SE high-confidence set"].append(fid)

    # format: one gene per locus_tag, features inside their gene
    for tag, fs in by_tag.items():
        genes = [f for f in fs if f.type == "gene"]
        N["format: one gene per locus_tag"] += 1
        if tag is None or len(genes) != 1 or tag in tag_gene:
            P["format: one gene per locus_tag"].append(tag)
            continue
        tag_gene[tag] = name
        gl = genes[0].location
        span = sorted((int(p.start), int(p.end)) for p in gl.parts)
        for f in fs:
            N["format: features inside their gene, same strand"] += 1
            inside = all(any(a <= int(p.start) and int(p.end) <= b for a, b in span) for p in f.location.parts)
            if not inside or f.location.strand != gl.strand:
                P["format: features inside their gene, same strand"].append((tag, f.type))

# GFF features not in the EMBL: only identical copies (same gene, type, location),
# the 2nd mRNA of a 2-ORF transcript and that transcript's UTRs
present = Counter((g["type"], top(i), tuple(sorted(g["parts"]))) for i, g in gff.items() if i in seen_ids)
missing = [i for i, g in gff.items() if i not in seen_ids and
           g["type"] in ("gene", "mRNA", "CDS", "tRNA", "rRNA", "five_prime_UTR", "three_prime_UTR", "intron")
           and not (gff[i]["seq"] == mito_name and g["type"] == "mRNA")]
explained = Counter()
unexplained = []
for i in missing:
    g = gff[i]
    if present[(g["type"], top(i), tuple(sorted(g["parts"])))]:
        explained["identical copy"] += 1
    elif g["type"].endswith("UTR") and parent[i] in two_orf:
        explained["UTR of 2-ORF transcript"] += 1
    else:
        unexplained.append(i)
check("GFF: every GFF feature in the EMBL (or expected not)", unexplained, len(gff))
print("     GFF features not in the EMBL, expected: %s" % dict(explained))

N["sequences: every FASTA sequence has an entry"] = len(fasta)
if entries != len(fasta):
    P["sequences: every FASTA sequence has an entry"].append("%d entries, %d sequences" % (entries, len(fasta)))
N["sources: every high-confidence tRNA in the EMBL"] = len(hc)
P["sources: every high-confidence tRNA in the EMBL"] = sorted(hc - seen_ids)
for name in sorted(N):
    check(name, P[name], N[name])

print("\n%s" % ("ALL CHECKS PASSED" if not fails else "FAILED: " + "; ".join(fails)))
sys.exit(1 if fails else 0)
