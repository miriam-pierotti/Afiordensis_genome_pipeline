#!/bin/bash
#SBATCH --cpus-per-task=50
#SBATCH --mem=120G
#SBATCH --partition=bigmem
#SBATCH --time=0-15:00
#SBATCH --job-name=trnascan-se
#SBATCH -o /nfs/scratch/pierotmi/stdout/trnascan-se.%j.out
#SBATCH -e /nfs/scratch/pierotmi/stdout/trnascan-se.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=miriam.pierotti@vuw.ac.nz

# tRNA genes on the nuclear genome with tRNAscan-SE 2.0 (eukaryotic models),
# then EukHighConfidenceFilter to keep high-confidence tRNA genes (removes
# tRNA-derived repeats). Output is read by 08_1_merge_annotations.sh.

set -eo pipefail

ROOT=/nfs/scratch/pierotmi
PROJECT=black_coral
VERSION=4
REF_GENOME=$ROOT/projects/$PROJECT/ref_genome/$VERSION/output/final_assembly
GENOME=$REF_GENOME/AntFio1_genome.fasta
OUTPUT_DIR=$ROOT/projects/$PROJECT/rna/output/trnascan-se

# Recreate the directory
rm -rf "$OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR"

source $(conda info --base)/etc/profile.d/conda.sh
conda activate trnascan-se

tRNAscan-SE -E -H --detail -Q \
            --thread ${SLURM_CPUS_PER_TASK:-50} \
            -o "$OUTPUT_DIR/AntFio1_trna.out" \
            -f "$OUTPUT_DIR/AntFio1_trna.ss" \
            -m "$OUTPUT_DIR/AntFio1_trna.stats" \
            -b "$OUTPUT_DIR/AntFio1_trna.bed" \
            -j "$OUTPUT_DIR/AntFio1_trna.gff" \
            -a "$OUTPUT_DIR/AntFio1_trna_seqs.fasta" \
            "$GENOME"

# Keep only high-confidence tRNA genes (removes tRNA-derived repeats)
EukHighConfidenceFilter -r \
    -i "$OUTPUT_DIR/AntFio1_trna.out" \
    -s "$OUTPUT_DIR/AntFio1_trna.ss" \
    -o "$OUTPUT_DIR" \
    -p AntFio1_trna_hc

conda deactivate

echo "tRNAs predicted: $(tail -n +4 $OUTPUT_DIR/AntFio1_trna.out | wc -l)"
echo "high confidence: $(tail -n +4 $OUTPUT_DIR/AntFio1_trna_hc.out | grep -c 'high confidence set')"
