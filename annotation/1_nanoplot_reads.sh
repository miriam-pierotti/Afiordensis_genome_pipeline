#!/bin/bash
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --partition=quicktest
#SBATCH --time=0-5:00:00
#SBATCH --job-name=nanoplot_pychopper
#SBATCH -o /nfs/scratch/pierotmi/stdout/nanoplot_pychopper.%j.out
#SBATCH -e /nfs/scratch/pierotmi/stdout/nanoplot_pychopper.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=miriam.pierotti@vuw.ac.nz

ROOT=/nfs/scratch/$USER
PROJECT=$1
TYPE=$2

INPUT=$ROOT/projects/$PROJECT/$TYPE/data
OUTPUT=$ROOT/projects/$PROJECT/$TYPE/data

########## RUN NANOPLOT ON ALL FASTQ FILES ##########
echo "Launching Nanoplot"
source activate nanoplot
    NanoPlot -t 16 --fastq $INPUT/$PROJECT'_'$TYPE'_mergedreads'.fastq.gz -o $OUTPUT/$PROJECT'_'$TYPE'_mergedreads_nanoplot' 
echo "Nanoplot completed"
conda deactivate