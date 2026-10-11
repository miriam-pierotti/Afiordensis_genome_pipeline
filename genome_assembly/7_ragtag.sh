#!/bin/bash
#SBATCH --cpus-per-task=50
#SBATCH --mem=32G
#SBATCH --partition=bigmem
#SBATCH --time=0-15:00
#SBATCH --job-name=ragtag
#SBATCH -o /nfs/scratch/pierotmi/stdout/ragtag.%j.out
#SBATCH -e /nfs/scratch/pierotmi/stdout/ragtag.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=miriam.pierotti@vuw.ac.nz           

#variables
CONFIG_FILE=$1
source "$CONFIG_FILE"

#OUT_DIR_CORRECT=$ROOT/projects/$PROJECT/$TYPE/$VERSION/output/ragtag_correct
OUT_DIR_SCAFFOLD=$ROOT/projects/$PROJECT/$TYPE/$VERSION/output/ragtag_scaffold
REFERENCE=$ROOT/data/P_pennacea_dataset/ncbi_dataset/data/GCA_964036295.1_jaPluPenn1.1_genomic.fna
QUERY=$ROOT/projects/$PROJECT/$TYPE/$VERSION/output/medaka/$PROJECT'_consensus.fasta'

#Prepare output directories
rm -rf $OUT_DIR_SCAFFOLD
mkdir -p $OUT_DIR_SCAFFOLD 
echo "Output directories prepared"

#Activate conda environment
source activate ragtag
echo "Conda Ragtag environment activated"

#Run RagTag scaffold
#Homology-based assembly scaffolding: Order and orient sequences in 'query.fa' by comparing them to sequences in 'reference.fa'
#https://github.com/malonge/RagTag/wiki/scaffold
echo "Running RagTag scaffold"
ragtag.py scaffold $REFERENCE $QUERY -t 50 -o $OUT_DIR_SCAFFOLD
echo "RagTag scaffold completed"

conda deactivate
echo "Conda environment deactivated"
echo "RagTag pipeline completed successfully"