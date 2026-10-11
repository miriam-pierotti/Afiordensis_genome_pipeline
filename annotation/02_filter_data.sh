#!/bin/bash
#SBATCH --ntasks=10
#SBATCH --mem=30G
#SBATCH --partition=parallel
#SBATCH --time=4-0:00
#SBATCH --job-name=filter_data
#SBATCH -o /nfs/scratch/pierotmi/stdout/filter_data.%j.out
#SBATCH -e /nfs/scratch/pierotmi/stdout/filter_data.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=miriam.pierotti@vuw.ac.nz

#variables
ROOT=/nfs/scratch/$USER
PROJECT=$1
TYPE=$2
Q=$3
SHORT=$4


#paths
INPUT=$ROOT/projects/$PROJECT/$TYPE/data/$PROJECT'_rna_mergedreads'.fastq.gz
OUTPUT=$ROOT/projects/$PROJECT/$TYPE/output/filter_data

mkdir -p $OUTPUT

########## FILTERING FOR HIGH QUALITY READS ##########
source activate chopper

module load GCC/8.3.0 
module load SAMtools/1.10 
module load HTSlib/1.10.2

echo "Modules loaded"

###run chopper to filter low quality reads and short reads
echo "Filtering at $SHORT kbp"
chopper -q $Q  -t 1 -l $SHORT'000' -i $INPUT | bgzip -c > $OUTPUT/$PROJECT'_chopper_'$SHORT'kbp_q'$Q'.fastq.gz'
echo "Filtering at $SHORT kbp done."
echo "Quality and Lenght filtering completed"
conda deactivate

source activate nanoplot
    NanoPlot -t 16 --fastq $OUTPUT/$PROJECT'_chopper_'$SHORT'kbp_q'$Q'.fastq.gz' -o $OUTPUT/$PROJECT'_chopper_'$SHORT'kbp_q'$Q'_nanoplot' 
echo "Nanoplot completed"
conda deactivate
