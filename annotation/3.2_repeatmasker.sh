#!/bin/bash
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --partition=bigmem
#SBATCH --time=1-0:00
#SBATCH --job-name=Repeat_masker
#SBATCH -o /nfs/scratch/pierotmi/stdout/Repeat_masker.%j.out
#SBATCH -e /nfs/scratch/pierotmi/stdout/Repeat_masker.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=miriam.pierotti@vuw.ac.nz

module load GCC/10.2.0
module load OpenMPI/4.0.5
module load Singularity/3.10.2
module load h5py

#Create paths
ROOT=/nfs/scratch/pierotmi
PROJECT=$1
TYPE=$2
VERSION=$3
TYPEOFASSEMBLY=$4
ASSEMBLY=$PROJECT
READ_DIR=$ROOT/projects/$PROJECT/$TYPE/$VERSION/output
RAGTAG_DIR=$READ_DIR/final_assembly
REPEAT_DIR=$READ_DIR/repeat_masker

#export PATH to libraries created in previous section
export LIBDIR=$ROOT/Repeat_masker_libs/Libraries

#make an output dir for repeat_masker processes
mkdir -p $REPEAT_DIR


mkdir -p $REPEAT_DIR/libraries

cd $ROOT/Repeat_masker_libs

# 1. build database 
singularity exec --bind $ROOT,$REPEAT_DIR,$RAGTAG_DIR,$READ_DIR dfam-tetools-latest.sif \
    BuildDatabase -name "$ROOT/projects/$PROJECT/$TYPE/$VERSION/output/repeat_masker/libraries/${PROJECT}_db" \
                  "$ROOT/projects/$PROJECT/$TYPE/$VERSION/output/final_assembly/AntFio1.fasta" 

# #2. Identify TEs
singularity exec --bind $ROOT,$REPEAT_DIR,$RAGTAG_DIR,$READ_DIR dfam-tetools-latest.sif \
    RepeatModeler -threads 32 \
                  -LTRStruct \
                  -database "$ROOT/projects/$PROJECT/$TYPE/$VERSION/output/repeat_masker/libraries/${ASSEMBLY}_db"
#                 -engine ncbi #remove because it is not working with the latest version of RepeatModeler

# #3. check if taxonomic group is available
singularity exec --bind $ROOT,$REPEAT_DIR,$RAGTAG_DIR,$READ_DIR dfam-tetools-latest.sif \
    famdb.py -i "$ROOT/Repeat_masker_libs/Libraries/famdb" lineage -ad "Hexacorallia"
# #extract repeats
singularity exec --bind $ROOT,$REPEAT_DIR,$RAGTAG_DIR,$READ_DIR dfam-tetools-latest.sif \
    famdb.py -i "$ROOT/Repeat_masker_libs/Libraries/famdb" \
             families -ad --add-reverse-complement "Hexacorallia" \
    > "$REPEAT_DIR/Hexacorallia_library.fa"

# #4. merge dfam and repeatmodeler databases
cat $REPEAT_DIR/libraries/${ASSEMBLY}_db-families.fa $REPEAT_DIR/Hexacorallia_library.fa > $REPEAT_DIR/libraries/combined_library.fa

# #5. run repeatmasker
singularity exec --bind $ROOT,$REPEAT_DIR,$RAGTAG_DIR,$READ_DIR dfam-tetools-latest.sif \
    RepeatMasker -pa 32 \
                 -dir "$ROOT/projects/$PROJECT/$TYPE/$VERSION/output/repeat_masker/results" \
                 -xsmall \
                 -lib "$ROOT/projects/$PROJECT/$TYPE/$VERSION/output/repeat_masker/libraries/combined_library.fa" \
                 "$ROOT/projects/$PROJECT/$TYPE/$VERSION/output/final_assembly/AntFio1.fasta"

# cp $REPEAT_DIR/results/assembly.fasta.masked $REPEAT_DIR/assembly.softmasked.ordered.fasta

#6. convert softmasked to hardmasked
#sed 's/(acgt)/N/g' $REPEAT_DIR/assembly.softmasked.ordered.fasta > $REPEAT_DIR/assembly.hardmasked.ordered.fasta

# SBATCH BRAKER
# sbatch /nfs/scratch/pierotmi/scripts/phd_blackcoral/genome_pipeline/10.1_braker.sh /nfs/scratch/pierotmi/projects/black_coral/ref_genome/4/output/config.sh