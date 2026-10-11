#!/bin/bash
#SBATCH --cpus-per-task=50
#SBATCH --mem=500G
#SBATCH --partition=gpu
#SBATCH --time=1-0:00
#SBATCH --job-name=NCBI_FCS_GX
#SBATCH -o /nfs/scratch/pierotmi/stdout/NCBI_FCS_GX.%j.out
#SBATCH -e /nfs/scratch/pierotmi/stdout/NCBI_FCS_GX.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=miriam.pierotti@vuw.ac.nz

# FOR MORE INFO ON THIS TOOL, PLEASE VISIT THE FOLLOWING LINK:
# https://github.com/ncbi/fcs/wiki/FCS-GX-quickstart

CONFIG_FILE=$1
source "$CONFIG_FILE"


#set environmental variables
LOCAL_DB=$ROOT/bin/NCBI-FCS-GX
INPUT_DIR=$ROOT/projects/black_coral/ref_genome/4/output/medaka
OUT_DIR=$ROOT/projects/black_coral/ref_genome/4/output/fcs_gx
SOURCE_DB_MANIFEST=https://ncbi-fcs-gx.s3.amazonaws.com/gxdb/latest/all.manifest

#load modules
module load GCC/10.2.0
module load OpenMPI/4.0.5
module load Singularity/3.10.2


#move to the tool directory
cd $ROOT/bin/NCBI-FCS-GX

export FCS_DEFAULT_IMAGE=$ROOT/bin/NCBI-FCS-GX/fcs-gx.sif

#copy database to tmp
python3 fcs.py db get --mft $LOCAL_DB/gxdb/all.manifest --dir /tmp/pierotmi/ncbi/gxdb

#check the database
python3 fcs.py db check --mft $SOURCE_DB_MANIFEST --dir $LOCAL_DB/gxdb
python3 fcs.py db check --mft $LOCAL_DB"/gxdb/all.manifest" --dir /tmp/pierotmi/ncbi/gxdb

#set environmental variable for the database on the RAM
GXDB_LOC=/tmp/pierotmi/ncbi

#create the output directory
mkdir -p $OUT_DIR

#run the screen genome - 86253 is the tax id for A.Fiordensis
python3 ./fcs.py screen genome --fasta $INPUT_DIR/$PROJECT'_consensus.fasta' --out-dir $OUT_DIR --gx-db $GXDB_LOC/gxdb --tax-id 86253 