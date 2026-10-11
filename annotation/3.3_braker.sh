#!/bin/bash
#SBATCH --cpus-per-task=50
#SBATCH --mem=256G
#SBATCH --partition=bigmem
#SBATCH --time=0-15:00
#SBATCH --job-name=braker
#SBATCH -o /nfs/scratch/pierotmi/stdout/braker.%j.out
#SBATCH -e /nfs/scratch/pierotmi/stdout/braker.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=miriam.pierotti@vuw.ac.nz

#variables
CONFIG_FILE=$1
source "$CONFIG_FILE"

#load modules
module load GCC/10.2.0
module load OpenMPI/4.0.5
module load Singularity/3.10.2

OUT_DIR=$ROOT/projects/$PROJECT/$TYPE/$VERSION/output/braker
PROTEIN=$ROOT/Braker_libs/Metazoa.fa
GENOME=$ROOT/projects/$PROJECT/$TYPE/$VERSION/output/repeat_masker/results/AntFio1.fasta.masked

export BRAKER_SIF=$ROOT/Braker_libs/braker3.sif

rm -rf $OUT_DIR
mkdir -p $OUT_DIR
cd $OUT_DIR

#run braker
singularity exec --bind $OUT_DIR,$PROTEIN,$GENOME $BRAKER_SIF 'braker.pl' --genome=$GENOME --prot_seq=$PROTEIN --threads=32 --gff3 
