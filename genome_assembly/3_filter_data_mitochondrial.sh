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


CONFIG_FILE=$1
source "$CONFIG_FILE"

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
echo "Filtering at '$MITO_SHORT'kbp"
bgzip -cd $INPUT | chopper -q $MITO_Q  -t 1 -l $MITO_SHORT'000' | bgzip -c > $OUTPUT/$PROJECT'_mitochondrial_'$MITO_SHORT'kbp_q'$MITO_Q'.fastq.gz'
echo "Filtering at '$MITO_SHORT'kbp done."
echo "Quality and Lenght filtering completed"
conda deactivate

