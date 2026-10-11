#!/bin/bash
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --partition=quicktest
#SBATCH --time=0-5:00:00
#SBATCH --job-name=ena_embl
#SBATCH -o /nfs/scratch/pierotmi/stdout/ena_embl.%j.out
#SBATCH -e /nfs/scratch/pierotmi/stdout/ena_embl.%j.err
#SBATCH --mail-type=END,FAIL

# Converts the nuclear GFF (08_1) and the mitogenome GFF (08_2) to EMBL with
# EMBLmyGFF3, separately (nuclear: linear, translation table 1; mitogenome:
# circular, table 4), post-fixes what EMBLmyGFF3 cannot write (08_3_fix_embl.py),
# and writes the nuclear, mito and combined annotations (GFF3 + EMBL) next to the
# assembly FASTAs in final_assembly. The combined EMBL, gzipped, is the flat file
# of the ENA genome submission, with the chromosome list that tells ENA the
# circular entry is the mitochondrion.
# ENA flat files may hold several entries, each with its own topology on its ID
# line (https://ena-docs.readthedocs.io/en/latest/submit/fileprep/flat-file-example.html).
#
# Requires conda environment: emblmygff3_py310 (EMBLmyGFF3 2.4, Biopython 1.88).
#   The env "emblmygff3" (EMBLmyGFF3 2, Biopython 1.70) crashes: bcbio-gff needs Biopython >= 1.80.

set -eo pipefail

#set variables
ROOT=/nfs/scratch/$USER
PROJECT=black_coral
SCRIPTS=$ROOT/scripts/phd_blackcoral/rna_pipeline

#ENA metadata
SPECIES="Antipathella fiordensis" # must match the organism registered with ENA
LOCUS_TAG=ANTFIO                  # locus tag prefix registered with the ENA study
PROJECT_ID=PRJEB127395            # ENA study accession (PRJEB...), NOT the Webin account ID
TAXONOMY=INV                      # ENA taxonomic division: invertebrates

#define directories
dir=$ROOT/projects/$PROJECT/rna/output
assembly=$ROOT/projects/$PROJECT/ref_genome/4/output/final_assembly
nuclear_gff=$dir/functional_annotation/AntFio1_complete_ena.gff3   # from 08_1
nuclear_fa=$assembly/AntFio1_genome.fasta
mito_gff=$dir/functional_annotation/AntFio1_mito_ena.gff3          # from 08_2
mito_fa=$assembly/AntFio1_mito.fasta
combined_fa=$assembly/AntFio1.fasta                                 # nuclear + mito

out_dir=$dir/ena_submission
submission=$ROOT/projects/$PROJECT/ena_submission
FLATFILE=$submission/afiordensis_1.embl.gz            # FLATFILE in manifest_genome.txt
CHROMLIST=$submission/chromosome_list.txt.gz          # CHROMOSOME_LIST in manifest_genome.txt

for f in "$nuclear_gff" "$nuclear_fa" "$mito_gff" "$mito_fa" "$combined_fa"; do
	if [ ! -s "$f" ]; then
		echo "ERROR: missing or empty file: $f" >&2
		exit 1
	fi
done
mkdir -p $out_dir/nuclear $out_dir/mito

#EMBLmyGFF3 reads translation_gff_feature_to_embl_feature.json from the working
#directory: exon features are left out (with isoforms they are duplicated or abut,
#and the mRNA/CDS join() already holds them); the mitogenome also leaves out mRNA
#(the gene -> mRNA level is only there so EMBLmyGFF3 joins CDS pieces across introns)
write_json() {
	cat > $1/translation_gff_feature_to_embl_feature.json << EOF
{
 "five_prime_UTR": {"target": "5'UTR"},
 "three_prime_UTR": {"target": "3'UTR"},
 "exon": {"remove": true}$2
}
EOF
}
write_json $out_dir/nuclear ""
write_json $out_dir/mito ',
 "mRNA": {"remove": true}'

module purge
source activate emblmygff3_py310

#####################################################
# Step 1: nuclear genome (linear, translation table 1)
#####################################################

echo "=== Step 1: nuclear genome -> EMBL ==="
cd $out_dir/nuclear
EMBLmyGFF3 "$nuclear_gff" "$nuclear_fa" \
	--data_class STD --topology linear --molecule_type "genomic DNA" \
	--transl_table 1 --species "$SPECIES" --taxonomy $TAXONOMY \
	--locus_tag $LOCUS_TAG --project_id $PROJECT_ID --no_progress \
	-o raw.embl < /dev/null > emblmygff3.log 2>&1
python3 $SCRIPTS/08_3_fix_embl.py raw.embl "$nuclear_gff" fixed.embl

#####################################################
# Step 2: mitogenome (circular, translation table 4)
#####################################################

echo "=== Step 2: mitogenome -> EMBL ==="
#locus tags continue after the nuclear ones so they do not collide
last_locus=$(grep -o "${LOCUS_TAG}_LOCUS[0-9]*" $out_dir/nuclear/fixed.embl | sed 's/.*LOCUS//' | sort -n | tail -1)
cd $out_dir/mito
EMBLmyGFF3 "$mito_gff" "$mito_fa" \
	--data_class STD --topology circular --molecule_type "genomic DNA" \
	--transl_table 4 --species "$SPECIES" --taxonomy $TAXONOMY \
	--locus_tag $LOCUS_TAG --locus_numbering_start $((last_locus + 1)) \
	--project_id $PROJECT_ID --no_progress \
	-o raw.embl < /dev/null > emblmygff3.log 2>&1
python3 $SCRIPTS/08_3_fix_embl.py raw.embl "$mito_gff" fixed.embl
conda deactivate

#####################################################
# Step 3: annotation files next to the assembly FASTAs, flat file for ENA
#####################################################

#final_assembly gets, for each FASTA, its annotation in GFF3 and EMBL:
#  AntFio1_genome.fasta  AntFio1_genome.gff3  AntFio1_genome.embl   nuclear
#  AntFio1_mito.fasta    AntFio1_mito.gff3    AntFio1_mito.embl     mitogenome
#  AntFio1.fasta         AntFio1.gff3         AntFio1.embl          both
#AntFio1.embl, gzipped, is the flat file submitted to ENA.
echo "=== Step 3: annotation files in $assembly ==="
cp $nuclear_gff $assembly/AntFio1_genome.gff3
cp $mito_gff $assembly/AntFio1_mito.gff3
( echo "##gff-version 3"; grep -hv '^#' $nuclear_gff $mito_gff ) > $assembly/AntFio1.gff3
cp $out_dir/nuclear/fixed.embl $assembly/AntFio1_genome.embl
cp $out_dir/mito/fixed.embl $assembly/AntFio1_mito.embl
cat $assembly/AntFio1_genome.embl $assembly/AntFio1_mito.embl > $assembly/AntFio1.embl
gzip -c $assembly/AntFio1.embl > $FLATFILE
#ENA takes the organelle from the chromosome list ("chromosome" includes organelles:
#https://ena-docs.readthedocs.io/en/latest/submit/assembly/genome.html); the nuclear
#sequences are not listed, so they stay unplaced scaffolds
mito_name=$(grep -m1 '^>' $mito_fa | cut -c2- | cut -d' ' -f1)
printf "%s\tMT\tCircular-Chromosome\tMitochondrion\n" $mito_name | gzip > $CHROMLIST

#checks
#AntFio1.fasta = AntFio1_genome.fasta + AntFio1_mito.fasta, same names, order and sequences
seqs() { awk '/^>/ {if (s) print n "\t" s; n = substr($1, 2); s = ""; next} {s = s toupper($0)} END {print n "\t" s}' "$@" | md5sum; }
[ "$(seqs $combined_fa)" = "$(seqs $nuclear_fa $mito_fa)" ] \
	|| { echo "ERROR: $combined_fa is not $nuclear_fa + $mito_fa" >&2; exit 1; }
#one EMBL entry per sequence, in the FASTA order
diff <(grep '^>' $combined_fa | cut -c2- | cut -d' ' -f1) <(grep '^AC \* _' $assembly/AntFio1.embl | cut -c7-) > /dev/null \
	|| { echo "ERROR: EMBL entries do not match the sequences of $combined_fa" >&2; exit 1; }
#the same features in the combined GFF as in the two GFFs
[ "$(grep -vc '^#' $assembly/AntFio1.gff3)" -eq $(( $(grep -vc '^#' $nuclear_gff) + $(grep -vc '^#' $mito_gff) )) ] \
	|| { echo "ERROR: AntFio1.gff3 does not hold all features" >&2; exit 1; }
cmp -s <(zcat $FLATFILE) $assembly/AntFio1.embl || { echo "ERROR: $FLATFILE differs from AntFio1.embl" >&2; exit 1; }
awk '/^FT/ && length > 80 {bad++} END {if (bad) {print "ERROR: " bad " FT lines > 80 characters"; exit 1}}' $assembly/AntFio1.embl
if grep -q '/anticodon="' $assembly/AntFio1.embl; then echo "ERROR: quoted /anticodon left" >&2; exit 1; fi
echo "sequences in AntFio1.fasta: $(grep -c '^>' $combined_fa), entries in AntFio1.embl: $(grep -c '^ID ' $assembly/AntFio1.embl)"
echo "topology: $(grep '^ID ' $assembly/AntFio1.embl | cut -d';' -f3 | sort | uniq -c | tr -s ' ' | tr '\n' ',')"
echo "features:"; grep -o '^FT   [^ ]\+' $assembly/AntFio1.embl | sort | uniq -c
echo "EMBLmyGFF3 logs: $out_dir/nuclear/emblmygff3.log, $out_dir/mito/emblmygff3.log"
echo "=== Done ==="
ls -l $assembly/AntFio1_genome.{gff3,embl} $assembly/AntFio1_mito.{gff3,embl} $assembly/AntFio1.{gff3,embl}
echo "Flat file for ENA: $FLATFILE"
echo "Chromosome list: $CHROMLIST"
echo "Next: 08_4_check_embl.sh, then validate with webin-cli (ena_submission/submit.sh, -validate)"
