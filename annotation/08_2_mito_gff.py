#!/usr/bin/env python3
"""Step 08_2: make the mitogenome annotation (GeSeq 2.03) ENA-compatible.

Usage: 08_2_mito_gff.py geseq.gff3 mito.fasta out.gff3

Coordinates are taken unchanged from the GeSeq GFF; what changes:
- one gene -> mRNA -> exon/CDS model per protein-coding gene, so EMBLmyGFF3 joins
  the CDS pieces of the intron-containing genes (cox1, nad5) into one CDS;
  08_3_ena_embl.sh drops the mRNA and exon features from the EMBL, as usual for
  mitochondrial records. Introns are kept as intron features.
- CDS phase 0 (GeSeq writes phase 1 on every line, which would give
  /codon_start=2); every CDS is checked to translate (code 4) with a start
  codon, a terminal stop and no internal stop.
- gene -> tRNA and gene -> rRNA models.
- tRNA genes annotated twice (blatN and tRNAscan-SE at the same position) are
  kept once, from blatN.
- INSDC gene names and /product names (GENES below).
- nad2 crosses the origin (16676..17710 + 1..81): written as one gene/mRNA on
  the 3' part plus a CDS/exon piece at 1..81; 08_3_fix_embl.py turns it into
  join(16676..17710,1..81), which EMBLmyGFF3 cannot write.
"""
import re
import sys
from collections import defaultdict
from urllib.parse import unquote

geseq_gff, mito_fa, gff_out = sys.argv[1:4]

# GeSeq gene name -> (INSDC gene name, feature type, product)
GENES = {
    "rrn12": ("rrnS", "rRNA", "12S ribosomal RNA"),
    "rrn16": ("rrnL", "rRNA", "16S ribosomal RNA"),
    "trnM": ("trnM", "tRNA", "tRNA-Met"),
    "trnW": ("trnW", "tRNA", "tRNA-Trp"),
    "COX1": ("cox1", "CDS", "cytochrome c oxidase subunit I"),
    "COX2": ("cox2", "CDS", "cytochrome c oxidase subunit II"),
    "COX3": ("cox3", "CDS", "cytochrome c oxidase subunit III"),
    "CYTB": ("cob", "CDS", "cytochrome b"),
    "ATP6": ("atp6", "CDS", "ATP synthase F0 subunit 6"),
    "ATP8": ("atp8", "CDS", "ATP synthase F0 subunit 8"),
    "ND1": ("nad1", "CDS", "NADH dehydrogenase subunit 1"),
    "ND2": ("nad2", "CDS", "NADH dehydrogenase subunit 2"),
    "ND3": ("nad3", "CDS", "NADH dehydrogenase subunit 3"),
    "ND4": ("nad4", "CDS", "NADH dehydrogenase subunit 4"),
    "ND4L": ("nad4L", "CDS", "NADH dehydrogenase subunit 4L"),
    "ND5": ("nad5", "CDS", "NADH dehydrogenase subunit 5"),
    "ND6": ("nad6", "CDS", "NADH dehydrogenase subunit 6"),
}
# table 4 (mold/protozoan/coelenterate mitochondrial): TGA = Trp
STARTS4 = {"ATG", "ATA", "ATT", "ATC", "TTG", "TTA", "CTG", "GTG"}
STOPS4 = {"TAA", "TAG"}
COMP = str.maketrans("ACGT", "TGCA")

seqid, seq = None, []
for line in open(mito_fa):
    if line.startswith(">"):
        seqid = line[1:].split()[0]
    else:
        seq.append(line.strip().upper())
seq = "".join(seq)
length = len(seq)

# GeSeq features by gene name: {name: {type: [(start, end, strand, source)]}}
feats = defaultdict(lambda: defaultdict(list))
for line in open(geseq_gff):
    c = line.rstrip("\n").split("\t")
    if line.startswith("#") or len(c) != 9 or c[2] in ("region", "source"):
        continue
    if c[0] != seqid:
        sys.exit("ERROR: GFF seqid %s is not the FASTA sequence %s" % (c[0], seqid))
    gene = re.search(r"gene=([^;]+)", c[8]).group(1)
    feats[gene][c[2]].append((int(c[3]), int(c[4]), c[6], c[1]))

# drop tRNAscan-SE duplicates of blatN tRNAs (same position)
blatn = {(s, e, st) for g in feats.values() for s, e, st, src in g.get("tRNA", []) if src == "blatN"}
for name in list(feats):
    trnas = feats[name].get("tRNA", [])
    if trnas and all(src != "blatN" and (s, e, st) in blatn for s, e, st, src in trnas):
        print("dropped duplicate tRNA %s (same position as the blatN tRNA)" % name)
        del feats[name]
unknown = set(feats) - set(GENES)
if unknown:
    sys.exit("ERROR: GeSeq genes without an entry in GENES: %s" % ", ".join(sorted(unknown)))


def translate_ok(parts, strand):
    nt = "".join(seq[s - 1:e] for s, e in parts)
    if strand == "-":
        nt = nt.translate(COMP)[::-1]
    codons = [nt[i:i + 3] for i in range(0, len(nt) - 2, 3)]
    return (len(nt) % 3 == 0 and codons[0] in STARTS4 and codons[-1] in STOPS4
            and not any(c in STOPS4 for c in codons[:-1]))


rows = []
for name, f in feats.items():
    gname, ftype, product = GENES[name]
    gene_parts = sorted((s, e) for s, e, _, _ in f["gene"])
    strand = f["gene"][0][2]
    parts = sorted((s, e) for s, e, _, _ in f[ftype])
    # a feature across the origin: GeSeq gives the 3' part (ending at the sequence
    # end) and the 5' part (starting at 1) separately; order them 5' -> 3'
    wraps = len(parts) > 1 and parts[0][0] == 1 and parts[-1][1] == length
    if wraps:
        parts = parts[1:] + parts[:1]
        gene_span = parts[0]                      # 3' part; the rest is added in 08_3
    else:
        gene_span = (gene_parts[0][0], gene_parts[-1][1])
    gid = "gene-" + gname
    common = "gene=%s" % gname
    rows.append((gene_span[0], strand, [gene_span[0], gene_span[1], "gene", ".", "ID=%s;%s" % (gid, common)]))
    if ftype == "CDS":
        if not translate_ok(parts, strand):
            sys.exit("ERROR: %s CDS %s does not translate cleanly with code 4" % (gname, parts))
        rid = "rna-" + gname
        rows.append((gene_span[0], strand, [gene_span[0], gene_span[1], "mRNA", ".", "ID=%s;Parent=%s;%s" % (rid, gid, common)]))
        for n, (s, e) in enumerate(parts, 1):
            rows.append((gene_span[0], strand, [s, e, "exon", ".", "ID=exon-%s-%d;Parent=%s" % (gname, n, rid)]))
        for s, e in parts:
            rows.append((gene_span[0], strand, [s, e, "CDS", "0", "ID=cds-%s;Parent=%s;%s;product=%s" % (gname, rid, common, product)]))
    else:
        if len(parts) != 1:
            sys.exit("ERROR: %s %s has %d parts" % (gname, ftype, len(parts)))
        s, e = parts[0]
        rows.append((gene_span[0], strand, [s, e, ftype, ".", "ID=%s-%s;Parent=%s;%s;product=%s"
                                    % (ftype.lower(), gname, gid, common, product)]))
    for n, (s, e, _, _) in enumerate(sorted(f.get("intron", [])), 1):
        rows.append((gene_span[0], strand, [s, e, "intron", ".", "ID=intron-%s-%d;Parent=%s;%s" % (gname, n, gid, common)]))

with open(gff_out, "w") as out:
    out.write("##gff-version 3\n##sequence-region %s 1 %d\n" % (seqid, length))
    for _, strand, (s, e, ftype, phase, attrs) in sorted(rows, key=lambda r: r[0]):
        out.write("\t".join(map(str, [seqid, "GeSeq", ftype, s, e, ".", strand, phase, attrs])) + "\n")

count = defaultdict(int)
for _, _, r in rows:
    count[r[2]] += 1
print("mito features: " + ", ".join("%s %d" % kv for kv in sorted(count.items())))
