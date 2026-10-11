#!/bin/bash
#SBATCH --cpus-per-task=12
#SBATCH --mem=32G
#SBATCH --partition=quicktest
#SBATCH --time=0-3:00:00
#SBATCH --job-name=interpro
#SBATCH -o /nfs/scratch/pierotmi/stdout/interpro%j.out
#SBATCH -e /nfs/scratch/pierotmi/stdout/interpro.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=miriam.pierotti@vuw.ac.nz

#load modules
module load Java/11.0.20

#set variables
ROOT=/nfs/scratch/pierotmi/projects/black_coral/rna

#define directories
proteins=$ROOT/output/clean_transdecoder/bambu.extended_annotations.fa.transdecoder_clean.pep

out_dir=$ROOT/output/interpro/
out_file=$out_dir/interpro
mkdir -p $out_dir

cd /nfs/scratch/pierotmi/bin/interproscan

#run interpro
./interproscan.sh \
	-appl pfam \
	--disable-precalc \
	-f JSON \
	--goterms \
	--iprlookup \
	--pathways \
	--seqtype p \
	--cpu 22 \
	--input $proteins \
	--outfile $out_file.json

./interproscan.sh \
	-appl pfam \
	--disable-precalc \
	-f TSV \
	--goterms \
	--iprlookup \
	--pathways \
	--seqtype p \
	--cpu 22 \
	--input $proteins \
	--outfile $out_file.tsv

./interproscan.sh \
	-appl pfam \
	--disable-precalc \
	-f GFF3 \
	--goterms \
	--iprlookup \
	--pathways \
	--seqtype p \
	--cpu 22 \
	--input $proteins \
	--outfile $out_file.gff3



# -appl allows you to choose which analyses to run. If not specified all analyses will be run
# --disable-precalc disables use of precalculated match lookup service. All match calculations will be run locally. Not entirely sure what this means but can increase accuracy when working with novel sequences (at expense of increased compute time/resource use).
# --goterms switches on lookup of corresponding GO annotation
# --iprlookup Also include lookup of corresponding InterPro annotation in the TSV and GFF3 output formats.
# --pathways switch on lookup of corresponding Pathway annotation


