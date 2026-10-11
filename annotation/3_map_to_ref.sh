#!/bin/bash
#SBATCH --cpus-per-task=128
#SBATCH --mem=256G
#SBATCH --partition=gpu
#SBATCH --time=0-12:00:00
#SBATCH --job-name=map_to_ref
#SBATCH -o /nfs/scratch/pierotmi/stdout/map_to_ref.%j.out
#SBATCH -e /nfs/scratch/pierotmi/stdout/map_to_ref.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=miriam.pierotti@vuw.ac.nz

# Load minimap2 conda environment
source activate minimap2

module load GCC/12.3.0
module load SAMtools/1.19.2

ROOT=/nfs/scratch/$USER
PROJECT=$1
TYPE=$2
Q=$3
SHORT=$4

INPUT_DIR=$ROOT/projects/$PROJECT/$TYPE/output/filter_data
REFGENOME=$ROOT/projects/$PROJECT/ref_genome/4/output/final_assembly/AntFio1.fasta
OUT_DIR=$ROOT/projects/$PROJECT/$TYPE/output/reads_to_ref
mkdir -p $OUT_DIR

# Align reads to reference genome and produce SAM
minimap2 -ax splice -uf --secondary=no -k14 $REFGENOME $INPUT_DIR/$PROJECT'_chopper_'$SHORT'kbp_q'$Q.fastq.gz -t 128\
    | samtools sort -o $OUT_DIR/reads_$SHORT'kbp_q'$Q'_to_refgenome_sorted'.bam

# Index the BAM
samtools index $OUT_DIR/reads_$SHORT'kbp_q'$Q'_to_refgenome_sorted'.bam

conda deactivate

sbatch /nfs/scratch/pierotmi/scripts/phd_blackcoral/rna_pipeline/4_bambu.sh black_coral rna

