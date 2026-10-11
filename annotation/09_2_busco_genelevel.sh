#!/bin/bash
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=10G
#SBATCH --partition=quicktest
#SBATCH --time=0-5:00
#SBATCH --job-name=BUSCO_genelevel
#SBATCH -o /nfs/scratch/pierotmi/stdout/BUSCO_genelevel.%j.out
#SBATCH -e /nfs/scratch/pierotmi/stdout/BUSCO_genelevel.%j.err
#SBATCH --mail-type=END,FAIL

# Step 09_2: gene-level BUSCO on the longest protein per gene of the annotation
# submitted to ENA, written by 09_1_annotation_summary.py.
# Uses the shared, already downloaded lineages (offline) so the dataset
# versions match the other BUSCO runs of the project.

ASSEMBLY=/nfs/scratch/pierotmi/projects/black_coral/rna/output/final_stats/ena_genelevel/ena_gene_longest_protein.pep
BUSCO_DIR=/nfs/scratch/pierotmi/projects/black_coral/rna/output/final_stats/ena_genelevel
LINEAGES=/nfs/scratch/pierotmi/busco_downloads
mkdir -p $BUSCO_DIR
cd $BUSCO_DIR

module purge
source activate busco

for odb in metazoa_odb10 metazoa_odb12; do
	busco	-i $ASSEMBLY				\
			-o busco						\
			-l $LINEAGES/lineages/$odb		\
			-m proteins 					\
			--offline						\
			--cpu 32						\
			-f
	rm -rf $BUSCO_DIR/busco_$odb
	mkdir -p $BUSCO_DIR/busco_$odb
	mv -f busco $BUSCO_DIR/busco_$odb
	echo "BUSCO analysis for $odb completed."
done

conda deactivate

grep -H "C:" $BUSCO_DIR/busco_*/busco/short_summary.specific*.txt
echo "Gene-level BUSCO results are located in: $BUSCO_DIR"
