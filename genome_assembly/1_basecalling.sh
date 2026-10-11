#!/bin/bash
#SBATCH --cpus-per-task=64
#SBATCH --mem=128G
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --time=10:00:00
#SBATCH --job-name=basecalling
#SBATCH -o /nfs/scratch/pierotmi/stdout/basecalling.%j.out
#SBATCH -e /nfs/scratch/pierotmi/stdout/basecalling.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=miriam.pierotti@vuw.ac.nz

module load GCC/12.3.0
module load SAMtools/1.19.2

CONFIG_FILE=$1
source "$CONFIG_FILE"


INDIR=$ROOT/data/$PROJECT/pod5
OUTDIR=$ROOT/projects/$PROJECT/ref_genome/raw_data/bam
DORADO=$ROOT/bin/dorado-0.8.0-linux-x64/bin

mkdir -p $OUTDIR

#BASECALLING
#echo "Starting Basecalling"
$DORADO/dorado basecaller sup $INDIR/ > $OUTDIR/$PROJECT'.bam'
$DORADO/dorado duplex sup $INDIR/ > $OUTDIR/$PROJECT'_duplex.bam'
#echo "Basecalling completed"

#LOAD LIBRARIES FOR SAMTOOLS

echo "Loading modules"

module load GCC/12.3.0
module load SAMtools/1.19.2

echo "Modules loaded"

#EXTRACT DUPLEX
samtools view -b -h -d dx:1 $OUTDIR/$PROJECT'_duplex.bam' > $OUTDIR/$PROJECT'_trueduplex.bam'
samtools view -b -h -d dx:0 $OUTDIR/$PROJECT'_duplex.bam' > $OUTDIR/$PROJECT'_simplexnoduplexoffsprings.bam'
samtools view -b -h -d dx:-1 $OUTDIR/$PROJECT'_duplex.bam' > $OUTDIR/$PROJECT'_simplexduplexoffsprings.bam'
echo "Duplex extraction completed"

#CONVERT TO FASTQ
samtools fastq $OUTDIR/$PROJECT'_trueduplex.bam' > $OUTDIR/$PROJECT'_trueduplex.fastq'
samtools fastq $OUTDIR/$PROJECT'_simplexnoduplexoffsprings.bam' > $OUTDIR/$PROJECT'_simplexnoduplexoffsprings.fastq'
samtools fastq $OUTDIR/$PROJECT'_simplexduplexoffsprings.bam' > $OUTDIR/$PROJECT'_simplexduplexoffsprings.fastq'
samtools fastq $OUTDIR/$PROJECT'.bam' > $OUTDIR/$PROJECT'.fastq'
echo "Conversion to fastq completed"

#CALCULATE STATS
source activate nanoplot
    NanoPlot -t 64 --fastq $OUTDIR/$PROJECT.fastq -o $OUTDIR/$PROJECT'_nanoplot'
conda deactivate

#GZIP
gzip $OUTDIR/$PROJECT'_trueduplex.fastq'
gzip $OUTDIR/$PROJECT'_simplexnoduplexoffsprings.fastq'
gzip $OUTDIR/$PROJECT'_simplexduplexoffsprings.fastq'
gzip $OUTDIR/$PROJECT'.fastq'
echo "Gzip completed"

echo "Launching filter data"

### SBATCH FILTER DATA ###
# cd $ROOT/scripts/phd_blackcoral/genome_pipeline
# sbatch 2_filter_data.sh $CONFIG_FILE
# sbatch 2_filter_data_mitochondrial.sh $CONFIG_FILE
# echo "DONE"