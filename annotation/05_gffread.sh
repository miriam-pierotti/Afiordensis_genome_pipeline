#!/bin/bash
#SBATCH -a 1-4
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --partition=gpu
#SBATCH --time=0-5:00:00
#SBATCH --job-name=gffread_array
#SBATCH -o /nfs/scratch/pierotmi/stdout/gffread/gffread_array_%A_%a.out
#SBATCH -e /nfs/scratch/pierotmi/stdout/gffread/gffread_array_%A_%a.err
#SBATCH --mail-type=END,FAIL

ROOT=/nfs/scratch/$USER
PROJECT=black_coral
N=${SLURM_ARRAY_TASK_ID}

REFGENOME=$ROOT/projects/$PROJECT/ref_genome/4/output/final_assembly/AntFio1.fasta
INPUT=$ROOT/projects/$PROJECT/rna/output/bambu
OUTDIR=$ROOT/projects/$PROJECT/rna/output/gffread
mkdir -p $OUTDIR

# List GTF files and check if any exist
gtf_files=($INPUT/*.gtf)
if [ ${#gtf_files[@]} -eq 0 ] || [ ! -f "${gtf_files[$((N-1))]}" ]; then
    echo "Error: No GTF files found or invalid array index $N. Skipping task."
    exit 1
fi

# Get the N-th .gtf file (array is 0-indexed in bash)
gtf_file=${gtf_files[$((N-1))]}

# Extract just the filename without path and extension
base_name=$(basename "$gtf_file" .gtf)

# Create output filename with bambu prefix
output_file="$OUTDIR/bambu.${base_name}.fa"

echo "Processing: $gtf_file -> $output_file"

source activate gffread

gffread "$gtf_file" -g "$REFGENOME" -w "$output_file"
# Compress the output file with gzip (keep original)
echo "Compressing: $output_file -> $output_file.gz"
gzip -k "$output_file"

conda deactivate