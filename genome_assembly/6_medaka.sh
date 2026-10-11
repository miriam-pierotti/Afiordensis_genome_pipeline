#!/bin/bash
#SBATCH --cpus-per-task=50
#SBATCH --mem=1000G
#SBATCH --partition=bigmem
#SBATCH --time=10-0:00
#SBATCH --job-name=Medaka
#SBATCH -o /nfs/scratch/pierotmi/stdout/Medaka.%j.out
#SBATCH -e /nfs/scratch/pierotmi/stdout/Medaka.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=miriam.pierotti@vuw.ac.nz

module load Miniconda3/22.11.1-1

#parameters
CONFIG_FILE=$1
source "$CONFIG_FILE"

#PATHS
READ_DIR=$ROOT/projects/$PROJECT/$TYPE/raw_data/fastq
OUTPUT_DIR=$ROOT/projects/$PROJECT/$TYPE/$VERSION/

#activate conda environment
source activate medaka

#run medaka on R10.4.1 reads (PromethION)			
medaka_consensus 	-i $READ_DIR/$PROJECT.fastq.gz											\
					-d $OUTPUT_DIR/output/purgedups/$PROJECT'_purged.fa'					\
					-o $OUTPUT_DIR/output/medaka											\
					-t 10					 												\
					-m r1041_e82_400bps_sup_v4.3.0

#-m needs to be changed for new cell models

conda deactivate

cd $OUTPUT_DIR/output/medaka/
mv consensus.fasta $PROJECT'_consensus.fasta'

#SBATCH CIRCLEPLOT PURGEDUPS
# cd $ROOT/scripts/phd_blackcoral/herro_pipeline_v2/
# sbatch 6.1_medaka_stats.sh $CONFIG_FILE

# cd $ROOT/scripts/phd_blackcoral/herro_pipeline_v2/
# sbatch 7_FCS_GX.sh $CONFIG_FILE

# cd $ROOT/scripts/phd_blackcoral/herro_pipeline_v2/
# sbatch 8_mitohifi.sh $CONFIG_FILE

