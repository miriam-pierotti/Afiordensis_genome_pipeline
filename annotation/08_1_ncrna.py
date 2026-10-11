#!/usr/bin/env python3
"""Step 08_1 (rRNA and tRNA genes): build ENA-ready rRNA and tRNA features from
the barrnap and tRNAscan-SE results and merge them with the protein-coding GFF.

Usage: 08_1_ncrna.py coding.gff3 barrnap.gff trna.out trna_hc.out trna.gff trna_hc.ss
       exclude_region exclude_ids out.gff3 overlaps.tsv

rRNA (barrnap):
- keeps complete hits only (drops every hit barrnap marks "(partial)")
- drops hits inside exclude_region (chr:start-end; the inverted rDNA units at the
  start of the AntFio1_12 array, a likely assembly fold-back artifact, see
  projects/black_coral/rdna_check/README.md)
- IDs <chromosome>.rrna<n>, numbered along each chromosome
- barrnap's 28S hit starts on top of the 5.8S: the 28S is trimmed to start
  after the 5.8S
tRNA (tRNAscan-SE):
- keeps the EukHighConfidenceFilter "high confidence set" (drops isotype
  mismatches, undetermined isotypes and unexpected anticodons)
- product=tRNA-Xxx (iMet -> tRNA-Met with an initiator note, SeC -> tRNA-Sec)
- a truncated end (tRNAscan-SE note trunc_start/trunc_end) is marked partial
  (partial=true, start_range/end_range), like partial CDS
- anticodon in the INSDC format (pos:<location>,aa:<amino acid>,seq:<codon>),
  from the genomic anticodon positions in the .ss file; commas are written as
  %2C so the GFF parser keeps the value in one piece
Both: gene -> rRNA/tRNA -> exon; loci whose IDs are in exclude_ids are removed;
ncRNAs overlapping a protein-coding CDS are reported; IDs must be unique; the
output is sorted by chromosome and locus start.
"""
import bisect
import re
import sys
from collections import defaultdict

(coding_gff, barrnap_gff, trna_out, trna_hc_out, trna_gff, trna_ss,
 exclude_region, exclude_ids, gff_out, overlap_out) = sys.argv[1:11]
exclude = {x for x in exclude_ids.split(",") if x}
ex_chr, ex_range = exclude_region.split(":")
ex_start, ex_end = map(int, ex_range.split("-"))

AA = {"iMet": "Met", "SeC": "Sec"}   # tRNAscan-SE isotype -> INSDC amino acid


def read_gff(path):
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) != 9:
                continue
            cols[3], cols[4] = int(cols[3]), int(cols[4])
            cols[8] = [list(kv.partition("=")[::2])
                       for kv in cols[8].strip().strip(";").split(";") if kv]
            rows.append(cols)
    return rows


def get(attrs, key):
    for k, v in attrs:
        if k == key:
            return v
    return None


def natural(s):
    return [int(t) if t.isdigit() else t for t in re.split(r"(\d+)", s)]


# ---- rRNA from barrnap: one locus = [gene, rRNA, exon]
rrna_loci = []
n_partial = n_region = 0
number = defaultdict(int)
for r in read_gff(barrnap_gff):
    if r[2] != "rRNA":
        continue
    if "(partial)" in get(r[8], "product"):    # e.g. product=5S ribosomal RNA (partial)
        n_partial += 1
        continue
    if r[0] == ex_chr and r[3] >= ex_start and r[4] <= ex_end:
        n_region += 1
        continue
    number[r[0]] += 1
    rid = "%s.rrna%d" % (r[0], number[r[0]])
    attrs = [["ID", rid], ["Parent", rid + ".gene"], ["product", get(r[8], "product")]]
    rrna_loci.append([
        [r[0], r[1], "gene", r[3], r[4], ".", r[6], ".", [["ID", rid + ".gene"]]],
        [r[0], r[1], "rRNA", r[3], r[4], r[5], r[6], ".", attrs],
        [r[0], r[1], "exon", r[3], r[4], ".", r[6], ".", [["ID", rid + ".exon1"], ["Parent", rid]]],
    ])

# trim 28S predictions that start on top of a 5.8S
s58 = [(l[1][0], l[1][6], l[1][3], l[1][4]) for l in rrna_loci if "5.8S" in get(l[1][8], "product")]
n_trim = 0
for locus in rrna_loci:
    r = locus[1]
    if "28S" not in get(r[8], "product"):
        continue
    for seq, strand, s, e in s58:
        if seq != r[0] or strand != r[6] or s > r[4] or e < r[3]:
            continue
        # only trim when the 5.8S sits at the 5' end of the 28S
        if strand == "+" and s - r[3] <= 50:
            new = (e + 1, r[4])
        elif strand == "-" and r[4] - e <= 50:
            new = (r[3], s - 1)
        else:
            print("WARNING: %s overlaps a 5.8S away from its 5' end, not trimmed" % get(r[8], "ID"))
            continue
        for f in locus:
            f[3], f[4] = new
        n_trim += 1

# ---- tRNA from tRNAscan-SE: high-confidence set only
def trna_notes(path):
    """tRNA ID -> Note column of a tRNAscan-SE .out file."""
    notes = {}
    with open(path) as fh:
        for line in fh.readlines()[3:]:
            c = [x.strip() for x in line.rstrip("\n").split("\t")]
            notes["%s.trna%s" % (c[0], c[1])] = c[-1]
    return notes


all_trna = trna_notes(trna_out)
hc = {i: n for i, n in trna_notes(trna_hc_out).items() if "high confidence set" in n}

anticodon = {}
current = None
with open(trna_ss) as fh:
    for line in fh:
        m = re.match(r"^(\S+) \(\d+-\d+\)", line)
        if m:
            current = m.group(1)
            continue
        m = re.search(r"Anticodon: (\w+) at \d+-\d+ \((\d+)-(\d+)\)", line)
        if m and current:
            anticodon[current] = (m.group(1), int(m.group(2)), int(m.group(3)))

trna_rows = defaultdict(list)
for r in read_gff(trna_gff):
    tid = get(r[8], "ID") if r[2] == "tRNA" else get(r[8], "Parent")
    if tid in hc:
        trna_rows[tid].append(r)
trna_loci = []
n_trunc = 0
for tid, rows in trna_rows.items():
    r = [x for x in rows if x[2] == "tRNA"][0]
    iso, codon = get(r[8], "isotype"), get(r[8], "anticodon")
    ss_codon, a, b = anticodon[tid]
    lo, hi = min(a, b), max(a, b)
    if ss_codon.upper() != codon.upper() or hi - lo != 2 or not r[3] <= lo <= hi <= r[4]:
        sys.exit("ERROR: anticodon of %s in the .ss file does not match the GFF" % tid)
    loc = "complement(%d..%d)" % (lo, hi) if r[6] == "-" else "%d..%d" % (lo, hi)
    attrs = [["ID", tid], ["Parent", tid + ".gene"],
             ["anticodon", "(pos:%s%%2Caa:%s%%2Cseq:%s)" % (loc, AA.get(iso, iso), codon.lower())],
             ["product", "tRNA-" + AA.get(iso, iso)]]
    if iso == "iMet":
        attrs.append(["note", "initiator methionine tRNA"])
    gene = [r[0], r[1], "gene", r[3], r[4], ".", r[6], ".", [["ID", tid + ".gene"]]]
    # tRNAscan-SE notes a truncated 5'/3' end as trunc_start:N / trunc_end:N
    p5, p3 = "trunc_start" in hc[tid], "trunc_end" in hc[tid]
    low, high = (p5, p3) if r[6] == "+" else (p3, p5)
    ends = [["partial", "true"]] if low or high else []
    ends += [["start_range", ".,%d" % r[3]]] if low else []
    ends += [["end_range", "%d,." % r[4]]] if high else []
    if ends:
        n_trunc += 1
    gene[8] += ends
    attrs += ends
    trna = r[:8] + [attrs]
    exons = [x[:8] + [[["ID", get(x[8], "ID")], ["Parent", tid]]] for x in rows if x[2] == "exon"]
    trna_loci.append([gene, trna] + exons)

# ---- protein-coding loci, grouped under their gene
coding = read_gff(coding_gff)
parent = {get(r[8], "ID"): (get(r[8], "Parent") or "").split(",")[0] or None for r in coding}


def root(i):
    while parent.get(i):
        i = parent[i]
    return i


coding_loci = defaultdict(list)
for r in coding:
    coding_loci[root(get(r[8], "ID") or (get(r[8], "Parent") or "").split(",")[0])].append(r)

# ---- IDs unique, excluded loci removed
all_ids = defaultdict(int)
for locus in list(coding_loci.values()) + rrna_loci + trna_loci:
    for r in locus:
        if r[2] != "CDS" and get(r[8], "ID"):    # CDS pieces share their ID
            all_ids[get(r[8], "ID")] += 1
clash = [i for i, n in all_ids.items() if n > 1]
if clash:
    sys.exit("ERROR: IDs used more than once: %s" % ", ".join(sorted(clash)[:10]))
loci, excluded = [], set()
for locus in list(coding_loci.values()) + rrna_loci + trna_loci:
    hit = {get(r[8], "ID") for r in locus} & exclude
    if hit:
        excluded |= hit
    else:
        loci.append(locus)
for x in sorted(exclude - excluded):
    sys.exit("ERROR: excluded ID %s not found" % x)

# ---- ncRNAs overlapping a CDS
products = {(get(r[8], "Parent") or ""): get(r[8], "product") for r in coding if r[2] == "CDS"}
cds = defaultdict(list)
for r in coding:
    if r[2] == "CDS":
        cds[r[0]].append((r[3], r[4], r[6], get(r[8], "Parent")))
max_len = defaultdict(int)
for seq in cds:
    cds[seq].sort()
    max_len[seq] = max(e - s + 1 for s, e, _, _ in cds[seq])
starts = {seq: [c[0] for c in v] for seq, v in cds.items()}
overlaps = {}
for locus in loci:
    for r in locus:
        if r[2] not in ("rRNA", "tRNA"):
            continue
        lst = cds.get(r[0], [])
        i = bisect.bisect_left(starts.get(r[0], []), r[3] - max_len[r[0]])
        while i < len(lst) and lst[i][0] <= r[4]:
            s, e, strand, mrna = lst[i]
            if e >= r[3]:
                key = (get(r[8], "ID"), mrna)
                bp = min(e, r[4]) - max(s, r[3]) + 1
                prev = overlaps.get(key)
                overlaps[key] = [get(r[8], "product"), r[0], r[3], r[4], r[6], strand,
                                 (prev[6] if prev else 0) + bp, products.get(mrna)]
            i += 1
with open(overlap_out, "w") as fh:
    fh.write("ncRNA_ID\tncRNA_product\tseqid\tstart\tend\tstrand\tmRNA_ID\t"
             "CDS_strand\tCDS_overlap_bp\tmRNA_product\n")
    for (nid, mrna), v in sorted(overlaps.items()):
        fh.write("\t".join(map(str, [nid] + v[:5] + [mrna] + v[5:])) + "\n")

# ---- write, sorted by chromosome and locus start
loci.sort(key=lambda l: (natural(l[0][0]), min(r[3] for r in l), max(r[4] for r in l),
                         natural(get(l[0][8], "ID"))))
with open(gff_out, "w") as out:
    out.write("##gff-version 3\n")
    for locus in loci:
        for r in locus:
            attrs = ";".join("%s=%s" % (k, v) if v != "" else k for k, v in r[8])
            out.write("\t".join(map(str, r[:8])) + "\t" + attrs + "\n")

print("rRNA: barrnap complete hits kept      %d (partial dropped %d, in %s dropped %d)"
      % (len(rrna_loci), n_partial, exclude_region, n_region))
print("rRNA: 28S trimmed after the 5.8S      %d" % n_trim)
print("tRNA: predicted / high confidence     %d / %d (truncated, marked partial: %d)"
      % (len(all_trna), len(trna_loci), n_trunc))
print("excluded loci (%s):  %d" % (exclude_ids, len(excluded)))
print("ncRNA-CDS overlaps (see report)       %d (same strand: %d)"
      % (len(overlaps), sum(1 for v in overlaps.values() if v[4] == v[5])))
