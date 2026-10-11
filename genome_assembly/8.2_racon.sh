#!/bin/bash
#SBATCH --ntasks=16
#SBATCH --mem=256G
#SBATCH --partition=bigmem
#SBATCH --time=0-6:00
#SBATCH --job-name=racon_mitogenome
#SBATCH -o /nfs/scratch/pierotmi/stdout/racon.%j.out
#SBATCH -e /nfs/scratch/pierotmi/stdout/racon.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=miriam.pierotti@vuw.ac.nz

#variables
CONFIG_FILE=$1
source "$CONFIG_FILE"

#paths
REF=$ROOT/projects/$PROJECT/$TYPE/$VERSION/output/filter_data
OUT_DIR=$ROOT/projects/$PROJECT/$TYPE/$VERSION/output/racon
REF_DIR=$ROOT/projects/$PROJECT/$TYPE/$VERSION/output/mitohifi
RACON=$ROOT/projects/$PROJECT/$TYPE/$VERSION/output/racon

rm -rf $OUT_DIR
mkdir -p $OUT_DIR

rm -rf $RACON
mkdir -p $RACON

# concatenate refs
cat $ROOT/projects/$PROJECT/$TYPE/$VERSION/output/medaka/$PROJECT'_consensus.fasta' $ROOT/projects/$PROJECT/$TYPE/$VERSION/output/mitohifi/final_mitogenome.fasta > $OUT_DIR/combined_genomes.fasta

source activate minimap2
#map the fastq reads (q20.1kb plus) to this concatenated reference
minimap2 -ax lr:hq -t 16 $OUT_DIR/combined_genomes.fasta $REF/$PROJECT'_mitochondrial_'$MITO_SHORT'kbp_q'$MITO_Q'.fastq' > $OUT_DIR/nuc_mito.sam
conda deactivate

module load GCC/8.3.0 
module load SAMtools/1.10 
module load HTSlib/1.10.2

#sort and index to allow pulling out of mito contig only
samtools sort -@ 8 -o $OUT_DIR/nuc_mito.sorted.bam $OUT_DIR/nuc_mito.sam
samtools index $OUT_DIR/nuc_mito.sorted.bam

#Extract mito contig mappings only (because this is the only contig you want to polish)
samtools view -bh $OUT_DIR/nuc_mito.sorted.bam ptg000001l_rc_rotated > $OUT_DIR/mito_after_nucmapping.bam
#decompress to sam for racon
samtools view -h -o $OUT_DIR/mito_after_nucmapping.sam $OUT_DIR/mito_after_nucmapping.bam

source activate racon
#Run racon with all reads against the mapping and unpolished mito reference
racon $REF/$PROJECT'_mitochondrial_'$MITO_SHORT'kbp_q'$MITO_Q'.fastq' $OUT_DIR/mito_after_nucmapping.sam $ROOT/projects/$PROJECT/$TYPE/$VERSION/output/mitohifi/final_mitogenome.fasta > $OUT_DIR/$PROJECT'_mito_corrected.fasta'
conda deactivate
