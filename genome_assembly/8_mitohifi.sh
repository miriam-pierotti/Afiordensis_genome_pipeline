#!/bin/bash
#SBATCH --cpus-per-task=64
#SBATCH --mem=320G
#SBATCH --partition=gpu
#SBATCH --time=1-00:00
#SBATCH --job-name=mitohifi
#SBATCH -o /nfs/scratch/pierotmi/stdout/mitohifi.%j.out
#SBATCH -e /nfs/scratch/pierotmi/stdout/mitohifi.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=miriam.pierotti@vuw.ac.nz

# Paths
CONFIG_FILE=$1
source "$CONFIG_FILE"

module load GCC/10.2.0
module load OpenMPI/4.0.5
module load Singularity/3.10.2
module load Java


MITOHIFI=$ROOT/bin/MitoHiFi/src
INPUT_DATA=$ROOT/projects/$PROJECT/$TYPE/raw_data/bam
#REF_GENOME=$ROOT/projects/$PROJECT/$TYPE/$VERSION/output/hifiasm_$SHORT'k_'$ULTRALONG'k'
REF_GENOME=$ROOT/projects/$PROJECT/$TYPE/$VERSION/output/ragtag_scaffold
OUTPUT_DIR=$ROOT/projects/$PROJECT/$TYPE/$VERSION/output/mitohifi
INPUT_REF=$ROOT/data/P_pennacea_dataset/ncbi_dataset/data

rm -rf $OUTPUT_DIR

mkdir -p $OUTPUT_DIR

cd $OUTPUT_DIR #move to the output directory as results of mitohifi as stored in the current directory


########### SINGULARITY APPROACH ############
# source activate mitofinder 

module load GCC/10.2.0
module load OpenMPI/4.0.5
module load Singularity/3.10.2
module load Java

singularity exec --bind $MITOHIFI,$INPUT_DATA,$OUTPUT_DIR,$INPUT_REF,$REF_GENOME --home $OUTPUT_DIR \
  docker://ghcr.io/marcelauliano/mitohifi:master mitohifi.py \
  -f $INPUT_REF/mitochondrial.fasta \
  -g $INPUT_REF/sequence.gb \
  -t 64 \
  -o 5 \
  -c $REF_GENOME'/ragtag.scaffold.fasta' \
  -a animal 
  # -p 90 \
  # --circular-size 18000 \
  # --max-read-len 2

# conda deactivate
########### END SINGULARITY APPROACH ###########

########### CONDA APPROACH ############

# source activate mitohifi 

# module load GCC/10.2.0
# module load OpenMPI/4.0.5
# module load Singularity/3.10.2
# module load Java

# python $MITOHIFI/mitohifi.py \
#   -f $INPUT_REF/mitochondrial.fasta \
#   -g $INPUT_REF/sequence.gb \
#   -t 64 \
#   -o 5 \
#   -r $INPUT_DATA/$PROJECT'_mitochondrial_'$MITO_SHORT'kbp_q'$MITO_Q'.fastq' \
#   -a animal \
#   -p 80

# conda deactivate

####### SBATCH ########
# cd $ROOT/scripts/phd_blackcoral/herro_pipeline_v2
# sbatch 8.1_mito_minimap.sh $CONFIG_FILE   -c $REF_GENOME/$PROJECT.asm.bp.p_ctg.fa \ -r $INPUT_DATA/$PROJECT'_mitochondrial_'$MITO_SHORT'kbp_q'$MITO_Q'.fastq' \ -r $INPUT_DATA/$PROJECT'_uncorrected_'$SHORT'kbp_q'$Q'.fastq' \

