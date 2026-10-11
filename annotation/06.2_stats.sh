#!/bin/bash
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=10G
#SBATCH --partition=quicktest
#SBATCH --time=0-5:00
#SBATCH --job-name=BUSCO_final
#SBATCH -o /nfs/scratch/pierotmi/stdout/BUSCO_final.%j.out
#SBATCH -e /nfs/scratch/pierotmi/stdout/BUSCO_final.%j.err
#SBATCH --mail-type=END,FAIL

#directory where your medaka consensus.fasta is 
ASSEMBLY_DIR=/nfs/scratch/pierotmi/projects/black_coral/rna/output/clean_transdecoder

#output directory where output is written to
BUSCO_DIR=/nfs/scratch/pierotmi/projects/black_coral/rna/output/final_stats
mkdir -p $BUSCO_DIR

#deactivate all module otherwise busco will not work in the conda environment (to be checked)
module purge

#activate the conda environment
source activate busco

#RUN BUSCO FOR ODB10
busco			-i $ASSEMBLY_DIR'/bambu.extended_annotations.fa.transdecoder_clean.pep' 		\
				-o busco										\
				-l metazoa_odb10		 						\
				-m proteins 				 						\
				--cpu 32										\
				-f									

#move the folder from the script file to the output folder
mkdir -p $BUSCO_DIR/busco_metazoa_odb10
mv -f busco $BUSCO_DIR/busco_metazoa_odb10

echo "BUSCO analysis for odb10 completed."

#RUN BUSCO FOR ODB12
busco			-i $ASSEMBLY_DIR'/bambu.extended_annotations.fa.transdecoder_clean.pep' 		\
				-o busco										\
				-l metazoa_odb12		 						\
				-m proteins 				 						\
				--cpu 32										\
				-f									

#move the folder from the script file to the output folder
mkdir -p $BUSCO_DIR/busco_metazoa_odb12
mv -f busco $BUSCO_DIR/busco_metazoa_odb12

echo "BUSCO analysis for odb12 completed."

conda deactivate

echo "BUSCO analyses completed successfully."
echo "BUSCO results are located in: $BUSCO_DIR"

#####################################################
#####################################################

echo "Start Agat for final statistics"
OUTPUT=$BUSCO_DIR/agat
rm -rf $OUTPUT
REF=$ROOT/projects/$PROJECT/ref_genome/4/output/final_assembly/AntFio1_genome.fasta

module load GCC/11.2.0
module load OpenMPI/4.1.1
module load R/4.2.0

rm -rf $OUTPUT.txt

echo "Launching Agat"
source activate agat
    agat_sp_statistics.pl -gff $ASSEMBLY_DIR/bambu.extended_annotations.fa.transdecoder.final_clean.gff -o $OUTPUT -f $REF -d -p
echo "Agat completed"
conda deactivate