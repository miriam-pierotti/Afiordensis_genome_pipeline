#!/bin/bash
# Makes the mitogenome annotation (GeSeq, curated by Annabel) ENA-compatible:
# see 08_2_mito_gff.py for what changes. Fast, run with bash (no sbatch needed).

set -eo pipefail

ROOT=/nfs/scratch/$USER
PROJECT=black_coral
SCRIPTS=$ROOT/scripts/phd_blackcoral/rna_pipeline
assembly=$ROOT/projects/$PROJECT/ref_genome/4/output/final_assembly

mito_gff=$assembly/AntFio1_mito_geseq.gff3  # GeSeq 2.03 annotation (Annabel), seqid renamed to AntFio_mt
mito_fa=$assembly/AntFio1_mito.fasta
out_gff=$ROOT/projects/$PROJECT/rna/output/functional_annotation/AntFio1_mito_ena.gff3

for f in "$mito_gff" "$mito_fa"; do
	if [ ! -s "$f" ]; then
		echo "ERROR: missing or empty file: $f" >&2
		exit 1
	fi
done

python3 $SCRIPTS/08_2_mito_gff.py "$mito_gff" "$mito_fa" "$out_gff"
echo "Mito GFF for ENA: $out_gff"
