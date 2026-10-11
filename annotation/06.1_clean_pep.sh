#!/bin/bash
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --partition=gpu
#SBATCH --time=0-5:00:00
#SBATCH --job-name=clean_files
#SBATCH -o /nfs/scratch/pierotmi/stdout/clean_files_%j.out
#SBATCH -e /nfs/scratch/pierotmi/stdout/clean_files_%j.err
#SBATCH --mail-type=END,FAIL

ROOT=/nfs/scratch/$USER
PROJECT=black_coral

INPUT=$ROOT/projects/$PROJECT/rna/output/transdecoder
OUTDIR=$ROOT/projects/$PROJECT/rna/output/clean_transdecoder
mkdir -p $OUTDIR

# Clean the GFF3 file: remove empty lines
sed '/^$/d' $INPUT/bambu.extended_annotations.fa.transdecoder.gff3 > $OUTDIR/bambu.extended_annotations.fa.transdecoder.final_clean.gff3

# Clean the PEP file: remove '*' characters and empty lines
sed 's/\*//g' $INPUT/bambu.extended_annotations.fa.transdecoder.pep | sed '/^$/d' > $OUTDIR/bambu.extended_annotations.fa.transdecoder_clean.pep

