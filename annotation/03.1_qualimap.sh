#!/bin/bash
#SBATCH --cpus-per-task=32
#SBATCH --mem=512G
#SBATCH --partition=bigmem
#SBATCH --time=1-0:00
#SBATCH --job-name=qualimap
#SBATCH -o /nfs/scratch/pierotmi/stdout/qualimap.%j.out
#SBATCH -e /nfs/scratch/pierotmi/stdout/qualimap.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=matteo.collina@vuw.ac.nz




ROOT=/nfs/scratch/$USER
PROJECT=$1
TYPE=$2

source activate qualimap
qualimap bamqc -bam /nfs/scratch/pierotmi/projects/black_coral/rna/output/reads_to_ref/reads_1kbp_q10_to_refgenome_sorted.bam -outdir /nfs/scratch/pierotmi/projects/black_coral/rna/output/reads_to_ref/qualimap --java-mem-size=128G -outformat HTML
conda deactivate