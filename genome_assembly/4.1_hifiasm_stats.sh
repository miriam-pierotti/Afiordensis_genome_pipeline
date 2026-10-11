#!/bin/bash
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --partition=quicktest
#SBATCH --time=0-5:00:00
#SBATCH --job-name=hifiasm_stats
#SBATCH -o /nfs/scratch/pierotmi/stdout/hifiasm_stats.%j.out
#SBATCH -e /nfs/scratch/pierotmi/stdout/hifiasm_stats.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=miriam.pierotti@vuw.ac.nz

#set parameters
module load GCC/10.2.0
module load CUDA/11.1.1
module load OpenMPI/4.0.5

# Paths
CONFIG_FILE=$1
source "$CONFIG_FILE"

READ_DIR=$ROOT/projects/$PROJECT/$TYPE/$VERSION
OUT_DIR=$READ_DIR/output/hifiasm_$SHORT'k_'$ULTRALONG'k'


mkdir -p $OUT_DIR

########## CONVERT AND REORDER CONTIGS ##########

source activate seqkit  #Activate this conda env as it has all the necessary tools installed in.

#convert gfa to fasta
awk '/^S/{print ">"$2;print $3}' $OUT_DIR/$PROJECT.asm.bp.hap1.p_ctg.gfa > $OUT_DIR/$PROJECT.asm.bp.hap1.p_ctg.fa
awk '/^S/{print ">"$2;print $3}' $OUT_DIR/$PROJECT.asm.bp.hap2.p_ctg.gfa > $OUT_DIR/$PROJECT.asm.bp.hap2.p_ctg.fa
echo "Convertion from gfa to fasta done."

seqkit stats -a $OUT_DIR/$PROJECT.asm.bp.hap1.p_ctg.fa > $OUT_DIR/$PROJECT.asm.bp.hap1.p_ctg.tsv
seqkit stats -a $OUT_DIR/$PROJECT.asm.bp.hap2.p_ctg.fa > $OUT_DIR/$PROJECT.asm.bp.hap2.p_ctg.tsv
echo "Stats calculated."

echo "DONE"

conda deactivate

########## CALCULATE STATS ##########
source activate seqkit  #Activate this conda env as it has all the necessary tools installed in.

gfastats -f $OUT_DIR/$PROJECT.asm.bp.p_ctg.reordered.fasta              > $OUT_DIR/$PROJECT.asm.bp.p_ctg.reordered.summary.tsv
gfastats -f $OUT_DIR/$PROJECT.asm.bp.p_ctg.reordered.fasta --seq-report > $OUT_DIR/$PROJECT.asm.bp.p_ctg.reordered.seq_summary.tsv

echo "Seqkit stats done."

conda deactivate


####### SBATCH CIRCLEPLOT ########
cd $ROOT/scripts/phd_blackcoral/herro_pipeline_v2/circeplot_busco
sbatch circle_plot.sh $CONFIG_FILE
echo "SBATCH CIRCLEPLOT DONE."


######## CALCULATE QUAST STATS #######

QUAST_DIR=$ROOT/projects/$PROJECT/$TYPE/$VERSION/output/hifiasm_stats/quast
mkdir -p $QUAST_DIR
 
source activate quast
quast $OUT_DIR/$PROJECT.asm.bp.p_ctg.reordered.fasta -o $QUAST_DIR
conda deactivate