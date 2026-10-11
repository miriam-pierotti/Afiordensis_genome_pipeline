#!/bin/bash
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=32G
#SBATCH --partition=quicktest
#SBATCH --time=0-1:00
#SBATCH --job-name=busco_summaries
#SBATCH -o /nfs/scratch/pierotmi/stdout/busco_summaries.%j.out
#SBATCH -e /nfs/scratch/pierotmi/stdout/busco_summaries.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=miriam.pierotti@vuw.ac.nz

# Paths
CONFIG_FILE=$1
source "$CONFIG_FILE"

HIFIASM_DIR=$ROOT/projects/$PROJECT/$TYPE/$VERSION/output/hifiasm_stats
MEDAKA_DIR=$ROOT/projects/$PROJECT/$TYPE/$VERSION/output/medaka_stats
PURGEDUPS_DIR=$ROOT/projects/$PROJECT/$TYPE/$VERSION/output/purgedups_stats
RAGTAG_DIR=$ROOT/projects/$PROJECT/$TYPE/$VERSION/output/ragtag_stats
OUT_DIR=$ROOT/projects/$PROJECT/$TYPE/$VERSION/output/busco_summaries

HIFIASM_RENAME=HI
MEDAKA_RENAME=ME
PURGEDUPS_RENAME=PD
RAGTAG_RENAME=RT

mkdir -p $OUT_DIR

source activate busco

#module load R

cp $HIFIASM_DIR/busco_metazoa_odb10/short_summary.specific.metazoa_odb10.busco.txt $OUT_DIR/short_summary.specific.metazoa_odb10.$HIFIASM_RENAME'_odb10'.txt
cp $HIFIASM_DIR/busco_metazoa_odb12/short_summary.specific.metazoa_odb12.busco.txt $OUT_DIR/short_summary.specific.metazoa_odb12.$HIFIASM_RENAME'_odb12'.txt
echo "Hifiasm BUSCO summaries copied to $OUT_DIR"
cp $MEDAKA_DIR/busco_metazoa_odb10/short_summary.specific.metazoa_odb10.busco.txt $OUT_DIR/short_summary.specific.metazoa_odb10.$MEDAKA_RENAME'_odb10'.txt
cp $MEDAKA_DIR/busco_metazoa_odb12/short_summary.specific.metazoa_odb12.busco.txt $OUT_DIR/short_summary.specific.metazoa_odb12.$MEDAKA_RENAME'_odb12'.txt
echo "Medaka BUSCO summaries copied to $OUT_DIR"
cp $PURGEDUPS_DIR/busco_metazoa_odb10/short_summary.specific.metazoa_odb10.busco.txt $OUT_DIR/short_summary.specific.metazoa_odb10.purgedups_odb10.txt
cp $PURGEDUPS_DIR/busco_metazoa_odb12/short_summary.specific.metazoa_odb12.busco.txt $OUT_DIR/short_summary.specific.metazoa_odb12.purgedups_odb12.txt
echo "Purgedups BUSCO summaries copied to $OUT_DIR"
cp $RAGTAG_DIR/busco_metazoa_odb10/short_summary.specific.metazoa_odb10.busco.txt $OUT_DIR/short_summary.specific.metazoa_odb10.ragtag_odb10.txt
cp $RAGTAG_DIR/busco_metazoa_odb12/short_summary.specific.metazoa_odb12.busco.txt $OUT_DIR/short_summary.specific.metazoa_odb12.ragtag_odb12.txt
echo "Ragtag BUSCO summaries copied to $OUT_DIR"
cd $ROOT/scripts/phd_blackcoral/herro_pipeline_v2/circleplot_busco

python3 generate_plot.py -wd $OUT_DIR 

conda deactivate
