#!/bin/bash
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=5G
#SBATCH --partition=quicktest
#SBATCH --time=0-1:00
#SBATCH --job-name=medaka_stats
#SBATCH -o /nfs/scratch/pierotmi/stdout/medaka_stats.%j.out
#SBATCH -e /nfs/scratch/pierotmi/stdout/medaka_stats.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=miriam.pierotti@vuw.ac.nz

#parameters
CONFIG_FILE=$1
source "$CONFIG_FILE"

#directory where your purge haplotigs assembly.fasta is
ASSEMBLY_DIR=$ROOT/projects/$PROJECT/$TYPE/$VERSION/output/medaka

#output directory where output is written to
OUT_DIR=$ROOT/projects/$PROJECT/$TYPE/$VERSION/output/medaka_stats/quast
mkdir -p $OUT_DIR

source activate quast

quast $ASSEMBLY_DIR/$PROJECT'_consensus.fasta' -o $OUT_DIR

echo "QUAST COMPLETED"

conda deactivate

# #SBATCH CIRCLEPLOT MEDAKA
# cd $ROOT/scripts/phd_blackcoral/herro_pipeline_v2/circleplot_busco
# sbatch circle_plot_medaka.sh $CONFIG_FILE

