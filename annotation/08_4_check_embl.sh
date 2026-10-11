#!/bin/bash
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --partition=quicktest
#SBATCH --time=0-1:00:00
#SBATCH --job-name=check_embl
#SBATCH -o /nfs/scratch/pierotmi/stdout/check_embl.%j.out
#SBATCH -e /nfs/scratch/pierotmi/stdout/check_embl.%j.err
#SBATCH --mail-type=END,FAIL

# Checks the ENA flat file (08_3) against the GFFs and the tool outputs it was
# built from (see 08_4_check_embl.py). Run before webin-cli validation.
#
# Requires conda environment: emblmygff3_py310 (for Biopython)

set -eo pipefail

ROOT=/nfs/scratch/$USER
PROJECT=black_coral
SCRIPTS=$ROOT/scripts/phd_blackcoral/rna_pipeline
dir=$ROOT/projects/$PROJECT/rna/output
assembly=$ROOT/projects/$PROJECT/ref_genome/4/output/final_assembly

module purge
source activate emblmygff3_py310
python3 $SCRIPTS/08_4_check_embl.py \
	$ROOT/projects/$PROJECT/ena_submission/afiordensis_1.embl.gz \
	$dir/functional_annotation/AntFio1_complete_ena.gff3 \
	$dir/functional_annotation/AntFio1_mito_ena.gff3 \
	$assembly/AntFio1_genome.fasta \
	$assembly/AntFio1_mito.fasta \
	$dir/clean_transdecoder/bambu.extended_annotations.fa.transdecoder_clean.pep \
	$dir/eggnog/eggnog.emapper.annotations \
	$dir/interpro/interpro.tsv \
	$dir/trnascan-se/AntFio1_trna_hc.out
conda deactivate
