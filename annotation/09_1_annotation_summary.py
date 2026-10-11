#!/usr/bin/env python3
"""
Step 09_1: gene-level annotation summary of the annotation as submitted to ENA
(AntFio1_complete_ena.gff3, from 08_1), so reported numbers match the deposited
GFF3. Run after 08_1; then run 09_2_busco_genelevel.sh.

- Translates every mRNA's CDS from the genome and writes the longest protein
  per ENA gene (input for 09_2_busco_genelevel.sh).
- Links each ENA mRNA back to its TransDecoder protein (AGAT renamed the IDs)
  by identical CDS coordinates in the genome-projected TransDecoder GFF, and
  checks that the two protein sequences are identical.
- Counts eggNOG/InterProScan/KEGG/GO/Pfam per ENA gene with the same rules as
  the earlier 8_Functional_Annotation.R (now in _oldscripts), plus the
  rRNA/tRNA genes.

Writes rna/output/annotation_summary.csv and rna/output/final_stats/ena_genelevel/.
"""

import re
from collections import Counter, defaultdict
from pathlib import Path

BC = Path("/nfs/scratch/pierotmi/projects/black_coral")
FA_DIR = BC / "rna/output/functional_annotation"
ENA_GFF = FA_DIR / "AntFio1_complete_ena.gff3"
PROJ_GFF = FA_DIR / "AntFio1_transdecoder.genome.gff3"
PEP_IN = BC / "rna/output/clean_transdecoder/bambu.extended_annotations.fa.transdecoder_clean.pep"
EGGNOG_IN = BC / "rna/output/eggnog/eggnog.emapper.annotations"
INTERPRO_IN = BC / "rna/output/interpro/interpro.tsv"
GENOME_FA = BC / "ref_genome/4/output/final_assembly/AntFio1_genome.fasta"

OUT_DIR = BC / "rna/output/final_stats/ena_genelevel"
PEP_OUT = OUT_DIR / "ena_gene_longest_protein.pep"
SUMMARY_OUT = BC / "rna/output/annotation_summary.csv"
STATS_OUT = OUT_DIR / "ena_stats.txt"

CODON = {}
_b = "TCAG"
_aa = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
for i, a in enumerate(_b):
    for j, b in enumerate(_b):
        for k, c in enumerate(_b):
            CODON[a + b + c] = _aa[16 * i + 4 * j + k]
COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def iter_fasta(path):
    header, chunks = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header, chunks = line[1:], []
            else:
                chunks.append(line)
    if header is not None:
        yield header, "".join(chunks)


def attr(col9, key):
    m = re.search(rf"(?:^|;){key}=([^;]+)", col9)
    return m.group(1) if m else None


def read_gff(path):
    """mRNA -> gene, mRNA -> [(seqid, strand, start, end, phase)], ncRNA products."""
    mrna_gene, cds, ncrna = {}, defaultdict(list), Counter()
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            if f[2] == "mRNA":
                mrna_gene[attr(f[8], "ID")] = attr(f[8], "Parent")
            elif f[2] == "CDS":
                cds[attr(f[8], "Parent")].append(
                    (f[0], f[6], int(f[3]), int(f[4]), int(f[7]) if f[7] != "." else 0))
            elif f[2] in ("rRNA", "tRNA"):
                ncrna[(f[2], attr(f[8], "product"))] += 1
    return mrna_gene, cds, ncrna


def translate(parts, genome):
    seqid, strand = parts[0][0], parts[0][1]
    parts = sorted(parts, key=lambda p: p[2])
    nt = "".join(genome[seqid][s - 1:e] for _, _, s, e, _ in parts)
    if strand == "-":
        nt = nt.translate(COMP)[::-1]
        parts = parts[::-1]  # first CDS in transcription order carries the phase
    nt = nt[parts[0][4]:].upper()
    prot = "".join(CODON.get(nt[i:i + 3], "X") for i in range(0, len(nt) - 2, 3))
    return prot.rstrip("*")


def cds_key(parts):
    return tuple(sorted((s, st, a, b) for s, st, a, b, _ in parts))


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    ena_gene, ena_cds, ncrna = read_gff(ENA_GFF)
    _, proj_cds, _ = read_gff(PROJ_GFF)
    by_key = defaultdict(list)
    for pid, parts in proj_cds.items():
        by_key[cds_key(parts)].append(pid)

    # ENA mRNA -> TransDecoder protein ID(s) with the same CDS
    ena2td = {m: by_key.get(cds_key(p), []) for m, p in ena_cds.items()}
    unmatched = [m for m, t in ena2td.items() if not t]

    # Proteins from the genome, checked against TransDecoder's .pep
    genome = {h.split()[0]: s for h, s in iter_fasta(GENOME_FA)}
    td_pep = {h.split()[0]: s.rstrip("*") for h, s in iter_fasta(PEP_IN)}
    prot, mismatch = {}, []
    for m, parts in ena_cds.items():
        prot[m] = translate(parts, genome)
        for t in ena2td[m]:
            if td_pep[t] != prot[m]:
                mismatch.append(m)

    genes = sorted(set(ena_gene[m] for m in ena_cds))
    best = {}
    for m in ena_cds:
        g = ena_gene[m]
        if g not in best or len(prot[m]) > len(prot[best[g]]):
            best[g] = m
    with open(PEP_OUT, "w") as out:
        for g in genes:
            m = best[g]
            out.write(f">{g} mRNA={m} transdecoder={','.join(ena2td[m])} length={len(prot[m])}\n")
            for i in range(0, len(prot[m]), 60):
                out.write(prot[m][i:i + 60] + "\n")

    # Functional hits per TransDecoder protein (same rules as the old 8_Functional_Annotation.R)
    egg, og, kegg, ipr, pfam, go = set(), set(), set(), set(), set(), set()
    with open(EGGNOG_IN) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            egg.add(f[0])
            if f[4] not in ("", "-"):
                og.add(f[0])
            if f[11] not in ("", "-") or f[12] not in ("", "-"):
                kegg.add(f[0])
    with open(INTERPRO_IN) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            ipr.add(f[0])
            if f[3] == "Pfam":
                pfam.add(f[0])
            if len(f) > 13 and f[13] not in ("", "-"):
                go.add(f[0])

    gene_td = defaultdict(set)
    for m, t in ena2td.items():
        gene_td[ena_gene[m]].update(t)

    def n_genes(hit):
        return sum(1 for g in genes if gene_td[g] & hit)

    total = len(genes)
    n_ann = n_genes(kegg | go | pfam)
    rows = [
        ("Total", total),
        ("eggNOG/InterProScan hit", n_genes(egg | ipr)),
        ("eggNOG", n_genes(og)),
        ("KEGG", n_genes(kegg)),
        ("GO", n_genes(go)),
        ("Pfam", n_genes(pfam)),
        ("Annotated", n_ann),
        ("Unannotated", total - n_ann),
    ]
    with open(SUMMARY_OUT, "w") as out:
        out.write("Type,Number,Percentage\n")
        for name, n in rows:
            out.write(f"{name},{n},{round(100 * n / total, 2)}\n")

    stats = [
        f"ENA GFF3:                     {ENA_GFF}",
        f"Protein-coding genes:         {total}",
        f"mRNAs:                        {len(ena_cds)}",
        f"mRNAs without TransDecoder match: {len(unmatched)}",
        f"Proteins differing from .pep: {len(mismatch)}",
        f"Output fasta:                 {PEP_OUT}",
        "ncRNA genes:",
    ] + [f"  {t}\t{p}\t{n}" for (t, p), n in sorted(ncrna.items())] + [
        f"  rRNA total {sum(n for (t, _), n in ncrna.items() if t == 'rRNA')}, "
        f"tRNA total {sum(n for (t, _), n in ncrna.items() if t == 'tRNA')}",
    ]
    STATS_OUT.write_text("\n".join(stats) + "\n")
    print("\n".join(stats))
    print(SUMMARY_OUT.read_text())


if __name__ == "__main__":
    main()
