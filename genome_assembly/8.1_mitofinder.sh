#!/bin/bash
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --partition=gpu
#SBATCH --time=0-3:00
#SBATCH --job-name=mitofinder
#SBATCH -o /nfs/scratch/pierotmi/stdout/mitofinder.%j.out
#SBATCH -e /nfs/scratch/pierotmi/stdout/mitofinder.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=miriam.pierotti@vuw.ac.nz

# Paths
CONFIG_FILE=$1
source "$CONFIG_FILE"

source activate mitofinder 

module load GCC/10.2.0
module load OpenMPI/4.0.5


INPUT_DATA=$ROOT/projects/$PROJECT/$TYPE/$VERSION/output/racon_2
OUTPUT_DIR=$ROOT/projects/$PROJECT/$TYPE/$VERSION/output/mitofinder
INPUT_REF=$ROOT/data/M_japonica_dataset/ncbi_dataset/data

mkdir -p $OUTPUT_DIR

cd $INPUT_DATA

mitofinder -a $INPUT_DATA/polished_iter2.fasta -j mtDNA_black_coral_polished -r $INPUT_REF/sequence.gb -o 5 -p 32 

conda deactivate

#mv black_coral_mito_corrected.fasta black_coral_mito_corrected black_coral_mito_corrected_MitoFinder.log $OUTPUT_DIR

####### SBATCH CLEANUP ########
#cd $ROOT/scripts/phd_blackcoral/herro_pipeline_v2
#sbtach x_clean_up.sh $PROJECT $TYPE $VERSION $Q $SHORT $ULTRALONG 