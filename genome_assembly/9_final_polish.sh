#!/bin/bash
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=16G
#SBATCH --partition=quicktest
#SBATCH --time=0-1:00
#SBATCH --job-name=final_polish
#SBATCH -o /nfs/scratch/pierotmi/stdout/final_polish.%j.out
#SBATCH -e /nfs/scratch/pierotmi/stdout/final_polish.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=miriam.pierotti@vuw.ac.nz

# Load configuration parameters from the config file
CONFIG_FILE=$1
source "$CONFIG_FILE"

# Path to the curated mitogenome FASTA file
MITO=/nfs/scratch/pierotmi/projects/black_coral/ref_genome/4/mitogenome_annabel/Antfio_mito_curated.fasta

# Directory containing the ragtag scaffolded assembly
ASSEMBLY_DIR=$ROOT/projects/$PROJECT/$TYPE/$VERSION/output/ragtag_scaffold

# Output directory for final assembly files
OUT_DIR=$ROOT/projects/$PROJECT/$TYPE/$VERSION/output/final_assembly
mkdir -p $OUT_DIR
mkdir -p $OUT_DIR/temp

# Activate conda environment with seqkit
source activate seqkit

# === STEP 1: Sort and rename scaffolds ===

# Sort scaffolds by length in descending order (-l = by length, -r = reverse/descending)
seqkit sort -l -r $ASSEMBLY_DIR/ragtag.scaffold.fasta > $OUT_DIR/temp/ragtag.scaffold.sorted.fasta

# Convert sorted FASTA to tab format (header\tsequence)
seqkit fx2tab $OUT_DIR/temp/ragtag.scaffold.sorted.fasta > $OUT_DIR/temp/ragtag.scaffold.sorted.tab

# Extract only the headers (first column) from the tab file
cut -f1 $OUT_DIR/temp/ragtag.scaffold.sorted.tab > $OUT_DIR/temp/ragtag.scaffold.sorted_headers.txt

# Generate new sequential scaffold names (AntFio1_01, AntFio1_02, etc.)
# Using awk to create numbered names with zero-padding (%02d = 2 digits with leading zeros)
awk -F"\t" 'BEGIN{OFS="\t"; i=1}
{
    printf("AntFio1_%02d\n", i);
    i++;
}' $OUT_DIR/temp/ragtag.scaffold.sorted_headers.txt > $OUT_DIR/temp/ragtag.scaffold.sorted_new_headers.txt

# Create a TSV file mapping original scaffold names to new names for future reference
paste $OUT_DIR/temp/ragtag.scaffold.sorted_headers.txt $OUT_DIR/temp/ragtag.scaffold.sorted_new_headers.txt > $OUT_DIR/AntFio1_scaffold_name_correspondence.tsv

# Extract sequences (second column) from the sorted tab file
cut -f2 $OUT_DIR/temp/ragtag.scaffold.sorted.tab > $OUT_DIR/temp/ragtag.scaffold.sorted_seqs.txt

# Combine new headers with sequences
paste $OUT_DIR/temp/ragtag.scaffold.sorted_new_headers.txt $OUT_DIR/temp/ragtag.scaffold.sorted_seqs.txt > $OUT_DIR/temp/ragtag.scaffold.sorted_edited.tab

# Convert back to FASTA format with new headers
seqkit tab2fx $OUT_DIR/temp/ragtag.scaffold.sorted_edited.tab > $OUT_DIR/AntFio1_genome.fasta

# Generate assembly statistics for the nuclear genome
seqkit stats $OUT_DIR/AntFio1_genome.fasta > $OUT_DIR/AntFio1_genome.stats.txt

# === STEP 2: Process mitogenome ===

# Convert mitogenome to tab format to check/modify header
seqkit fx2tab $MITO > $OUT_DIR/temp/mito.tab

# Extract the mitogenome header
MITO_HEADER=$(cut -f1 $OUT_DIR/temp/mito.tab)

# Check if mitogenome has the correct header name
if [ "$MITO_HEADER" != "AntFio_mt" ]; then
    # If header is incorrect, extract sequence
    cut -f2 $OUT_DIR/temp/mito.tab > $OUT_DIR/temp/mito_seq.txt
    
    # Create file with correct header
    echo -e "AntFio_mt" > $OUT_DIR/temp/mito_header.txt
    
    # Combine correct header with sequence
    paste $OUT_DIR/temp/mito_header.txt $OUT_DIR/temp/mito_seq.txt > $OUT_DIR/temp/mito_corrected.tab
    
    # Convert back to FASTA format with corrected header
    seqkit tab2fx $OUT_DIR/temp/mito_corrected.tab > $OUT_DIR/temp/mito_final.fasta
else
    # If header is already correct, just copy the file
    cp $MITO $OUT_DIR/temp/mito_final.fasta
fi

# Copy the processed mitogenome to output directory as a separate file
cp $OUT_DIR/temp/mito_final.fasta $OUT_DIR/AntFio1_mito.fasta

# === STEP 3: Create final assembly with mitogenome ===

# Concatenate nuclear genome and mitogenome into final assembly file
cat $OUT_DIR/AntFio1_genome.fasta $OUT_DIR/temp/mito_final.fasta > $OUT_DIR/AntFio1.fasta

# Generate statistics for the complete assembly (nuclear + mitochondrial)
seqkit stats $OUT_DIR/AntFio1.fasta > $OUT_DIR/AntFio1.stats.txt

# Deactivate conda environment
conda deactivate
