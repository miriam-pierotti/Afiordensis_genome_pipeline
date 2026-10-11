#!/bin/bash
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=32G
#SBATCH --partition=quicktest
#SBATCH --time=0-1:00
#SBATCH --job-name=cleanup
#SBATCH -o /nfs/scratch/pierotmi/stdout/cleanup.%j.out
#SBATCH -e /nfs/scratch/pierotmi/stdout/cleanup.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=miriam.pierotti@vuw.ac.nz

#parameters
ROOT=/nfs/scratch/pierotmi
PROJECT=$1 #black_coral
TYPE=$2 #this variable could be ref_genome or population
VERSION=$3 #this variable could be 1, 2, 3 etc.
Q=$4 #quality score
SHORT=$5 
ULTRALONG=$6

# Set the to the directory you want to check
IN_DIR=$ROOT/projects/$PROJECT/$TYPE/$VERSION/output
ASSEMBLY_DIR=$ROOT/projects/$PROJECT/$TYPE/$VERSION/assemblies
STDOUT_DIR=$ROOT/stdout


# Count the number of directories in the specified path and add 1
NUMBER=$(( $(find "$ASSEMBLY_DIR" -mindepth 1 -maxdepth 1 -type d | wc -l) + 1 ))

# Print the folder count
echo "Number of folders: $NUMBER"

OUT_DIR=$ROOT/projects/$PROJECT/$TYPE/$VERSION/assemblies/$NUMBER'_assembly_q'$Q'_'$SHORT'kbp_'$ULTRALONG'kbp'
LOG_DIR=$ROOT/projects/$PROJECT/$TYPE/$VERSION/assemblies/$NUMBER'_assembly_q'$Q'_'$SHORT'kbp_'$ULTRALONG'kbp'/logs

mkdir -p $OUT_DIR
mkdir -p $LOG_DIR

mv $IN_DIR/* $OUT_DIR
mv $STDOUT_DIR/* $LOG_DIR

echo "Files moved"

echo "PIPELINE COMPLETED!"
