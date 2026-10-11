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
CONFIG_FILE=$1
source $CONFIG_FILE

#paths
INPUT=$ROOT/projects/$PROJECT/$TYPE/raw_data/fastq/$PROJECT.fastq.gz #true_duplex comes from the last assembly. We used the simplex at the end.
OUTPUT=$ROOT/projects/$PROJECT/$TYPE/$VERSION/output/filter_data

mkdir -p $OUTPUT

########## FILTERING FOR HIGH QUALITY READS ##########
source activate chopper

module load GCC/8.3.0 
module load SAMtools/1.10 
module load HTSlib/1.10.2

echo "Modules loaded"

###run chopper to filter low quality reads and short reads
echo "Filtering at '$SHORT'kbp"
bgzip -cd $INPUT | chopper -q $Q  -t 10 -l $SHORT'000' | bgzip -c > $OUTPUT/$PROJECT'_uncorrected_'$SHORT'kbp_q'$Q'.fastq.gz'
echo "Filtering at '$SHORT'kbp done."
echo "Quality and Lenght filtering completed"
conda deactivate

########## LAUNCH NANOPLOT ##########
echo "Launching Nanoplot"
cd $ROOT/scripts/phd_blackcoral/
sbatch 2.2_filter_data_stats.sh $CONFIG_FILE

########## FILTERING FOR ULTRALONG READS ##########
source activate chopper

echo "Loading modules for UL filtering"
module load GCC/8.3.0 
module load SAMtools/1.10 
module load HTSlib/1.10.2
echo "Modules loaded for UL filtering"

echo "Filtering at '$ULTRALONG'kbp"
bgzip -cd $INPUT | chopper -t 15 -l $ULTRALONG'000' | bgzip -c > $OUTPUT/$PROJECT'_ONT_simplex_'$ULTRALONG'kbp.fastq.gz'
conda deactivate

echo "Filtering at '$ULTRALONG'kbp done."
echo "Quality and Lenght filtering for UL reads completed"

######### LAUNCH HIFIASM  #########
# echo "Launching Dorado Correct"
# cd $ROOT/scripts/phd_blackcoral/herro_pipeline_v2
# sbatch 4_hifiasm.sh $CONFIG_FILE
# echo "Dorado Correct launched"


echo "DONE"

