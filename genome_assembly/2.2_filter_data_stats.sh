#!/bin/bash
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --partition=quicktest
#SBATCH --time=0-5:00:00
#SBATCH --job-name=filter_data_stats
#SBATCH -o /nfs/scratch/pierotmi/stdout/filter_data_stats.%j.out
#SBATCH -e /nfs/scratch/pierotmi/stdout/filter_data_stats.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=miriam.pierotti@vuw.ac.nz

CONFIG_FILE=$1
source "$CONFIG_FILE"

INPUT=$ROOT/projects/$PROJECT/$TYPE/$VERSION/output/filter_data
OUTPUT=$ROOT/projects/$PROJECT/$TYPE/$VERSION/output/filter_data

########## RUN NANOPLOT ON ALL FASTQ FILES ##########
echo "Launching Nanoplot"
source activate nanoplot
    NanoPlot -t 10 --fastq $INPUT/$PROJECT'_uncorrected_'$SHORT'kbp_q'$Q'.fastq.gz' -o $OUTPUT/$PROJECT'_uncorrected_'$SHORT'kbp_q'$Q'_nanoplot'
echo "Nanoplot completed"
conda deactivate