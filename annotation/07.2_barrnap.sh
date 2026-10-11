#!/bin/bash
#SBATCH --cpus-per-task=50
#SBATCH --mem=120G
#SBATCH --partition=bigmem
#SBATCH --time=0-15:00
#SBATCH --job-name=barrnap
#SBATCH -o /nfs/scratch/pierotmi/stdout/barrnap.%j.out
#SBATCH -e /nfs/scratch/pierotmi/stdout/barrnap.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=miriam.pierotti@vuw.ac.nz

# rRNA genes on the nuclear genome with barrnap (eukaryotic models).
# Output is read by 08_1_merge_annotations.sh.

set -eo pipefail

ROOT=/nfs/scratch/pierotmi
PROJECT=black_coral
VERSION=4
REF_GENOME=$ROOT/projects/$PROJECT/ref_genome/$VERSION/output/final_assembly
GENOME=$REF_GENOME/AntFio1_genome.fasta
OUTPUT_DIR=$ROOT/projects/$PROJECT/rna/output/barrnap

mkdir -p $OUTPUT_DIR
printf "Running barrnap on %s\n" "$GENOME"
source activate barrnap

barrnap --kingdom euk \
        --threads ${SLURM_CPUS_PER_TASK:-50} \
        --outseq $OUTPUT_DIR/AntFio1_rrna_seqs.fa \
        "$GENOME" > $OUTPUT_DIR/AntFio1_genome_rRNA.gff

conda deactivate

echo "rRNA predictions: $(grep -vc '^#' $OUTPUT_DIR/AntFio1_genome_rRNA.gff)"
