#!/bin/bash
#SBATCH --cpus-per-task=32
#SBATCH --mem=320G
#SBATCH --partition=gpu
#SBATCH --time=1-0:00
#SBATCH --job-name=Hifiasm
#SBATCH -o /nfs/scratch/pierotmi/stdout/Hifiasm.%j.out
#SBATCH -e /nfs/scratch/pierotmi/stdout/Hifiasm.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=miriam.pierotti@vuw.ac.nz


# variables
CONFIG_FILE=$1
source "$CONFIG_FILE"

READ_DIR=$ROOT/projects/$PROJECT/$TYPE/$VERSION
LONG_DIR=$READ_DIR/output/filter_data
OUT_DIR=$READ_DIR/output/hifiasm_$SHORT'k_'$ULTRALONG'k'

#set HiFiasm path
HIFIASM=$ROOT/bin/hifiasm

mkdir -p $OUT_DIR

#set modules
module load GCC/10.2.0
module load CUDA/11.1.1
module load OpenMPI/4.0.5

echo "Reading from correct directory"
INPUT_DIR=$READ_DIR/output/filter_data

########## RUN THE ASSEMBLY ##########

cd $HIFIASM

./hifiasm -o $OUT_DIR/$PROJECT.asm -t 32 --ul $LONG_DIR/$PROJECT'_ONT_simplex_'$ULTRALONG'kbp.fastq.gz' $INPUT_DIR/$PROJECT'_uncorrected_8kbp_q10.fastq.gz' --ont

echo "Assembly done."


########## CONVERT AND REORDER CONTIGS ##########

source activate seqkit  #Activate this conda env as it has all the necessary tools installed in.

#convert gfa to fasta
awk '/^S/{print ">"$2;print $3}' $OUT_DIR/$PROJECT.asm.bp.p_ctg.gfa > $OUT_DIR/$PROJECT.asm.bp.p_ctg.fa
echo "Convertion from gfa to fasta done."

#reorder and rename contigs
echo "Reordering contigs by length."
seqkit sort --by-length --reverse $OUT_DIR/$PROJECT.asm.bp.p_ctg.fa | seqkit replace --pattern '.+' --replacement 'Contig{nr}' > $OUT_DIR/$PROJECT.asm.bp.p_ctg.reordered.fasta
echo "Sort by lenght done."

conda deactivate

########## CREATE AN INDEXED FASTA ##########

source activate chopper #Activate this conda env as it has SAMtools working - no need for chopper itself

module load GCC/8.3.0 
module load SAMtools/1.10 
module load HTSlib/1.10.2

samtools faidx $OUT_DIR/$PROJECT.asm.bp.p_ctg.reordered.fasta
echo "Index done."
conda deactivate


# #SBATCH
# cd $ROOT/scripts/phd_blackcoral/herro_pipeline_v2
# sbatch 5_purgedups.sh $CONFIG_FILE 
# echo "SBATCH PURGEDUPS DONE."


# #SBATCH HIFIASM STATS
# cd $ROOT/scripts/phd_blackcoral/herro_pipeline_v2
# sbatch 4.1_hifiasm_stats.sh $CONFIG_FILE
# echo "SBATCH HIFIASM STATS DONE."

# echo "DONE"