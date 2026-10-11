#!/bin/bash
#SBATCH --cpus-per-task=8
#SBATCH --mem=256G
#SBATCH --partition=gpu
#SBATCH --time=1-0:00
#SBATCH --job-name=bambu
#SBATCH -o /nfs/scratch/pierotmi/stdout/bambu.%j.out
#SBATCH -e /nfs/scratch/pierotmi/stdout/bambu.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=matteo.collina@vuw.ac.nz

ROOT=/nfs/scratch/$USER
PROJECT=black_coral
TYPE=rna

OUT_DIR=$ROOT/projects/$PROJECT/$TYPE/output/bambu
mkdir -p $OUT_DIR

source activate r

# Force correct R library path
export R_LIBS_USER=/nfs/home/pierotmi/.conda/envs/r/lib/R/library
# Run R
Rscript /nfs/scratch/pierotmi/scripts/phd_blackcoral/rna_pipeline/bambu.r

conda deactivate