#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --partition=quicktest
#SBATCH --time=0-5:00:00
#SBATCH --job-name=launcher
#SBATCH -o /nfs/scratch/pierotmi/stdout/launcher.%j.out
#SBATCH -e /nfs/scratch/pierotmi/stdout/launcher.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=miriam.pierotti@vuw.ac.nz

############################################################
########## SET ALL THE VARIABLES FOR THE PIPELINE ##########
############################################################

#REFERENCE GENOME PARAMETERS
ROOT=/nfs/scratch/$USER
PROJECT=black_coral #set the name of the project. ES: black_coral
TYPE=ref_genome #set the type of the assembly. ES: ref_genome
VERSION=4 #set the version of the assembly. Use this number for convience or for running multiple versions at the same time.
Q=10 #set the quality threshold ES: 10 is equal to q10
SHORT=8 #set the size of the short reads in kbp ES: 1 is equal to 1000
ULTRALONG=10 #set the size of the ultralong reads in kbp ES: 1 is equal to 1000
BASECALLING=false #set to true if you want to run the basecalling step

#MITOCHONDRIAL PARAMETERS
MITO_Q=20
MITO_SHORT=1

############################################################

########## CLEAN DIRECTORIES BEFORE STARTING ##########
OUT_FILE='launcher.'${SLURM_JOB_ID}'.out'
ERR_FILE='launcher.'${SLURM_JOB_ID}'.err'

find $ROOT/stdout -type f ! -name $OUT_FILE ! -name $ERR_FILE -delete
rm -r $ROOT/projects/$PROJECT/$TYPE/$VERSION/output/*
echo "Directories cleaned"

########## SAVE THE CONFIG FILE TO SOURCE IN THE PIPELINE ##########
# Define config file path
CONFIG_FILE="$ROOT/projects/$PROJECT/$TYPE/$VERSION/output/config.sh"

# Save variables in shell script format
cat <<EOL > "$CONFIG_FILE"
#!/bin/bash
ROOT=$ROOT
PROJECT=$PROJECT
TYPE=$TYPE
VERSION=$VERSION
Q=$Q
SHORT=$SHORT
ULTRALONG=$ULTRALONG
BASECALLING=$BASECALLING
MITO_Q=$MITO_Q
MITO_SHORT=$MITO_SHORT
EOL

chmod +x "$CONFIG_FILE"  # Ensure it's executable

echo "Configuration saved to $CONFIG_FILE"

########## SAVE THE CONFIG FILE AS TXT FOR REFERENCE ##########
# Define output file
OUTPUT_FILE="$ROOT/projects/$PROJECT/$TYPE/$VERSION/output/config.txt"

# Save variables to a text file
echo "CONFIGURATION" > "$OUTPUT_FILE"
echo "Root: $ROOT" >> "$OUTPUT_FILE"
echo "Project: $PROJECT" >> "$OUTPUT_FILE"
echo "Type: $TYPE" >> "$OUTPUT_FILE"
echo "Version: $VERSION" >> "$OUTPUT_FILE"
echo "Q: $Q" >> "$OUTPUT_FILE"
echo "Short Read Lenghts: $SHORT" >> "$OUTPUT_FILE"
echo "Ultralong Read Lenghts: $ULTRALONG" >> "$OUTPUT_FILE"
echo "Basecalling: $BASECALLING" >> "$OUTPUT_FILE"
echo "Mito Q: $MITO_Q" >> "$OUTPUT_FILE"
echo "Mito Short Read Lenghts: $MITO_SHORT" >> "$OUTPUT_FILE"
echo "" >> "$OUTPUT_FILE"

# Show saved content
cat "$OUTPUT_FILE"

########## SBATCHING THE PIPELINE ##########

# if [ $BASECALLING == "true" ]; then
#     sbatch 1_basecalling.sh $CONFIG_FILE
# else
#     sbatch 2_filter_data.sh $CONFIG_FILE
#     sbatch 2_filter_data_mitochondrial.sh $CONFIG_FILE
# fi
