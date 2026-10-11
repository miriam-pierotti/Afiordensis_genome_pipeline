#!/bin/bash
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --partition=gpu
#SBATCH --time=0-5:00:00
#SBATCH --job-name=eggnog
#SBATCH -o /nfs/scratch/pierotmi/stdout/eggnog.%j.out
#SBATCH -e /nfs/scratch/pierotmi/stdout/eggnog.%j.err
#SBATCH --mail-type=END,FAIL


#load modules
module purge

#set variables
ROOT=/nfs/scratch/$USER
PROJECT=black_coral

#define directories
dir=$ROOT/projects/$PROJECT/rna/output/clean_transdecoder
#proteins=$dir/bambu.novelTranscripts.fa.transdecoder_clean.pep
proteins=$dir/bambu.extended_annotations.fa.transdecoder_clean.pep
#gff3=$dir/bambu.novelTranscripts.fa.transdecoder_clean.gff3
gff3=$ROOT/projects/$PROJECT/rna/output/clean_transdecoder/bambu.extended_annotations.fa.transdecoder.final_clean.gff3
out_dir=$ROOT/projects/$PROJECT/rna/output/eggnog
mkdir -p $out_dir

source activate eggnog
export EGGNOG_DATA_DIR=/nfs/scratch/pierotmi/bin/eggnog-mapper/data

## run
# emapper.py \
# 	--data_dir /nfs/scratch/pierotmi/bin/eggnog-mapper/data \
# 	-i $proteins \
# 	-o $out_dir \
# 	--decorate_gff $gff3 \
# 	--tax_scope auto \
# 	--cpu 12

emapper.py -i $proteins --itype proteins -o $out_dir/eggnog --cpu 16 --decorate_gff $gff3 --tax_scope auto --override

conda deactivate