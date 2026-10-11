#!/bin/bash
#SBATCH --ntasks=32
#SBATCH --mem=100G
#SBATCH --partition=bigmem
#SBATCH --time=10-0:00
#SBATCH --job-name=purgedups
#SBATCH -o /nfs/scratch/pierotmi/stdout/purgedups.%j.out
#SBATCH -e /nfs/scratch/pierotmi/stdout/purgedups.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=miriam.pierotti@vuw.ac.nz

#variables
CONFIG_FILE=$1
source "$CONFIG_FILE"

#paths
INPUT=$ROOT/projects/$PROJECT/$TYPE/raw_data/fastq/$PROJECT'.fastq.gz'
HIFIASM_ASS=$ROOT/projects/$PROJECT/$TYPE/$VERSION/output/hifiasm_$SHORT'k_'$ULTRALONG'k'/$PROJECT.asm.bp.p_ctg.reordered.fasta
OUTPUT=$ROOT/projects/$PROJECT/$TYPE/$VERSION/output/purgedups

mkdir -p $OUTPUT


module load GCC/8.3.0
module load SAMtools/1.10 
module load HTSlib/1.10.2

# Load purge_haplotigs conda environment because it has minimap2 installed
source activate purge_haplotigs

minimap2 -xasm20 -t 30 $HIFIASM_ASS $INPUT | gzip -c - > $OUTPUT/aligned.paf.gz

# Calculate base-level read coverage and produce necessary files 
cd $ROOT/bin/purge_dups

./bin/pbcstat $OUTPUT/aligned.paf.gz
./bin/calcuts PB.stat > cutoffs 2>$OUTPUT/calcuts.log

./bin/split_fa $HIFIASM_ASS > $OUTPUT/hifiasm_ass.split

minimap2 -xasm5 -DP $OUTPUT/hifiasm_ass.split $OUTPUT/hifiasm_ass.split | gzip -c - > $OUTPUT/hifiasm_ass.split.self.paf.gz

./bin/purge_dups -2 -T cutoffs -c PB.base.cov $OUTPUT/hifiasm_ass.split.self.paf.gz > dups.bed 2> purge_dups.log

./bin/get_seqs -e dups.bed $HIFIASM_ASS > $OUTPUT/purged.fa

conda deactivate

# Move Purge_Dups files to the correct folder.
cd $ROOT/bin/purge_dups
mv cutoffs dups.bed hap.fa PB.base.cov PB.cov.wig PB.stat purge_dups.log purged.fa $ROOT/projects/$PROJECT/$TYPE/$VERSION/output/purgedups/

echo "Purgedups done."

# Rename purged.fa in $PROJECT_purged.fa
cd $ROOT/projects/$PROJECT/$TYPE/$VERSION/output/purgedups/
mv purged.fa $PROJECT'_purged.fa'
echo "Renaming done"


# ############ SBATCH PURGEDUPS STATS ###############
# cd $ROOT/scripts/phd_blackcoral/herro_pipeline_v2
# sbatch 5.1_purgedups_stats.sh $CONFIG_FILE
# echo "Purgedups stats launched"

# ############ SBATCH MEDAKA ###############
# cd $ROOT/scripts/phd_blackcoral/herro_pipeline_v2
# sbatch 6_medaka.sh $CONFIG_FILE
# echo "Medaka launched"

echo "DONE"