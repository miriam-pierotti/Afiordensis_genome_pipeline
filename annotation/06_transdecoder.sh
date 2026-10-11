#!/bin/bash
#SBATCH -a 1-4
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --partition=gpu
#SBATCH --time=0-5:00:00
#SBATCH --job-name=transdecoder_array
#SBATCH -o /nfs/scratch/pierotmi/stdout/transdecoder/transdecoder_array_%A_%a.out
#SBATCH -e /nfs/scratch/pierotmi/stdout/transdecoder/transdecoder_array_%A_%a.err
#SBATCH --mail-type=END,FAIL

ROOT=/nfs/scratch/$USER
PROJECT=black_coral
N=${SLURM_ARRAY_TASK_ID}

INPUT=$ROOT/projects/$PROJECT/rna/output/gffread
OUTDIR=$ROOT/projects/$PROJECT/rna/output/transdecoder
mkdir -p $OUTDIR

# Get the N-th .fa file based on SLURM array index
fa_file=$(ls $INPUT/*.fa | sed -n "${N}p")

# Extract just the filename without path and extension
base_name=$(basename "$fa_file" .fa)

echo "Processing: $fa_file"

cd $OUTDIR

source activate transdecoder
TransDecoder.LongOrfs -t "$fa_file" #--complete_orfs_only
TransDecoder.Predict -t "$fa_file" --no_refine_starts
conda deactivate

INPUT_GFF=$OUTDIR/${base_name}.transdecoder.gff3

# source activate agat
#     agat_sp_keep_longest_isoform.pl --gff "$INPUT_GFF" -o "${base_name}.final.gff"
# conda deactivate