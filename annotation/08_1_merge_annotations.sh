#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --partition=quicktest
#SBATCH --time=0-3:00:00
#SBATCH --job-name=merge_annot
#SBATCH -o /nfs/scratch/pierotmi/stdout/merge_annot.%j.out
#SBATCH -e /nfs/scratch/pierotmi/stdout/merge_annot.%j.err
#SBATCH --mail-type=END,FAIL

# Merges all nuclear annotations into one GFF3 in genome coordinates, ready for
# EMBL conversion (08_3_ena_embl.sh):
#   protein-coding genes  TransDecoder (06) + eggNOG-mapper (07) + InterProScan (07.1)
#   rRNA genes            barrnap (07.2)
#   tRNA genes            tRNAscan-SE high-confidence set (07.3)
#
# Steps:
#   0. Check input files
#   1. Check that protein IDs match between the proteins, GFF, eggNOG and InterProScan
#   2. Project the TransDecoder ORFs from transcript to genome coordinates
#   3. Protein-coding genes: product, gene, EC_number, Dbxref, partial ends (08_1_protein_coding.py)
#   4. Add rRNA and tRNA genes, sort, report ncRNA-CDS overlaps (08_1_ncrna.py)
#   5. Check that no protein-coding feature was lost
#
# Requires conda environments: transdecoder, eggnog (for python3)

set -eo pipefail

#set variables
ROOT=/nfs/scratch/$USER
PROJECT=black_coral
SCRIPTS=$ROOT/scripts/phd_blackcoral/rna_pipeline

#options
ADD_GENE_NAMES=true               # add eggNOG Preferred_name as /gene (set false to skip)
# rRNA hits to leave out:
# - every hit in this region: the inverted (- strand) rDNA units at the start of the
#   AntFio1_12 array, a likely assembly fold-back artifact (projects/black_coral/rdna_check/README.md)
RRNA_EXCLUDE_REGION="AntFio1_12:6432234-6454703"
# - these loci (comma-separated IDs). AntFio1_12.rrna1 is a 28S hit ~1.5 kb longer than
#   the others, with a reverse transcriptase ORF (g27875.t1) inside it on the same strand:
#   an rDNA copy interrupted by an insertion (R2-like?), not an intact 28S.
EXCLUDE_NCRNA="AntFio1_12.rrna1"

#define directories
dir=$ROOT/projects/$PROJECT/rna/output
assembly=$ROOT/projects/$PROJECT/ref_genome/4/output/final_assembly
proteins=$dir/clean_transdecoder/bambu.extended_annotations.fa.transdecoder_clean.pep
gff3=$dir/clean_transdecoder/bambu.extended_annotations.fa.transdecoder.final_clean.gff3
transcripts_fa=$dir/gffread/bambu.extended_annotations.fa
bambu_gtf=$dir/bambu/extended_annotations.gtf
eggnog_ann=$dir/eggnog/eggnog.emapper.annotations
interpro_tsv=$dir/interpro/interpro.tsv
genome_fa=$assembly/AntFio1_genome.fasta
#GO term names, used to spot eggNOG "descriptions" that are GO terms, not protein names
#(downloaded 2026-10-02 from http://purl.obolibrary.org/obo/go/go-basic.obo)
go_obo=$dir/go/go-basic.obo
barrnap_gff=$dir/barrnap/AntFio1_genome_rRNA.gff
trna_out=$dir/trnascan-se/AntFio1_trna.out
trna_hc_out=$dir/trnascan-se/AntFio1_trna_hc.out
trna_gff=$dir/trnascan-se/AntFio1_trna.gff
trna_ss=$dir/trnascan-se/AntFio1_trna_hc.ss

out_dir=$dir/functional_annotation
mkdir -p $out_dir

orf_type_gff=$out_dir/transdecoder.stranded_only.gff3   # transcript coordinates, stranded ORFs only
genome_gff=$out_dir/AntFio1_transdecoder.genome.gff3     # genome coordinates
coding_gff=$out_dir/AntFio1_protein_coding.gff3
final_gff=$out_dir/AntFio1_complete_ena.gff3
report=$out_dir/product_name_report.tsv
overlap_report=$out_dir/ncrna_cds_overlaps.tsv

#####################################################
# Step 0: check input files
#####################################################

echo "=== Step 0: checking input files ==="
for f in "$proteins" "$gff3" "$transcripts_fa" "$bambu_gtf" "$eggnog_ann" "$interpro_tsv" \
	"$genome_fa" "$go_obo" "$barrnap_gff" "$trna_out" "$trna_hc_out" "$trna_gff" "$trna_ss"; do
	if [ ! -s "$f" ]; then
		echo "ERROR: missing or empty file: $f" >&2
		exit 1
	fi
	echo "OK: $f"
done

#####################################################
# Step 1: check that the IDs match
#####################################################

echo "=== Step 1: checking ID consistency ==="

grep '^>' "$proteins" | sed 's/^>//; s/[[:space:]].*//' | sort -u > $out_dir/ids_proteins.txt
awk -F'\t' '$3=="mRNA" {match($9, /ID=[^;]+/); print substr($9, RSTART+3, RLENGTH-3)}' "$gff3" \
	| sort -u > $out_dir/ids_gff_mrna.txt
cut -f1 "$interpro_tsv" | sort -u > $out_dir/ids_interpro.txt
grep -v '^#' "$eggnog_ann" | cut -f1 | sort -u > $out_dir/ids_eggnog.txt

echo "Proteins in .pep:                 $(wc -l < $out_dir/ids_proteins.txt)"
echo "mRNA features in GFF:             $(wc -l < $out_dir/ids_gff_mrna.txt)"
for f in proteins interpro eggnog; do
	n=$(wc -l < $out_dir/ids_$f.txt)
	missing=$(comm -23 $out_dir/ids_$f.txt $out_dir/ids_gff_mrna.txt | wc -l)
	echo "IDs in $f: $n, not an mRNA in the GFF: $missing"
	if [ "$missing" -gt 0 ]; then
		echo "ERROR: $f IDs do not match the mRNA IDs in $gff3" >&2
		exit 1
	fi
done

#####################################################
# Step 2: project the ORFs onto the genome
#####################################################

echo "=== Step 2: projecting TransDecoder ORFs to genome coordinates ==="

module purge
source activate transdecoder
td_util=$CONDA_PREFIX/opt/transdecoder/util
cd $out_dir
$td_util/gtf_to_alignment_gff3.pl "$bambu_gtf" > $out_dir/bambu_transcripts.alignment.gff3
#bambu leaves some transcripts unstranded (strand "."): their introns mix GT-AG on
#both strands (sense/antisense chimeras), and the projection puts any ORF on them
#on the minus strand, giving CDS with internal stops. Leave their ORFs out.
awk -F'\t' '$3=="transcript" && $7=="." {
	match($9, /transcript_id "[^"]*"/); print substr($9, RSTART+15, RLENGTH-16) }' \
	"$bambu_gtf" | sort -u > $out_dir/unstranded_transcripts.txt
awk -F'\t' 'NR==FNR {u[$1]; next} /^#/ || !($1 in u)' \
	$out_dir/unstranded_transcripts.txt "$gff3" > $orf_type_gff
echo "Unstranded bambu transcripts: $(wc -l < $out_dir/unstranded_transcripts.txt)" \
	"(ORFs left out: $(awk -F'\t' 'NR==FNR {u[$1]; next} ($1 in u) && $3=="mRNA"' \
	$out_dir/unstranded_transcripts.txt "$gff3" | wc -l))"
#ORFs on the opposite strand of a spliced transcript cannot be placed and are skipped
$td_util/cdna_alignment_orf_to_genome_orf.pl \
	$orf_type_gff \
	$out_dir/bambu_transcripts.alignment.gff3 \
	"$transcripts_fa" \
	> "$genome_gff" 2> $out_dir/genome_projection.log
conda deactivate
echo "ORFs projected: $(awk -F'\t' '$3=="mRNA"' $genome_gff | wc -l)" \
	"(skipped, antisense on a spliced transcript: $(grep -c 'Skipping entry' $out_dir/genome_projection.log || true))"

#####################################################
# Step 3: protein-coding genes
#####################################################

echo "=== Step 3: product, gene, EC_number, Dbxref and partial ends on the protein-coding genes ==="

module purge
source activate eggnog
python3 $SCRIPTS/08_1_protein_coding.py \
	"$genome_gff" "$orf_type_gff" "$eggnog_ann" "$interpro_tsv" "$go_obo" "$genome_fa" \
	"$coding_gff" "$report" "$ADD_GENE_NAMES"

#####################################################
# Step 4: rRNA and tRNA genes
#####################################################

echo "=== Step 4: adding rRNA and tRNA genes ==="

python3 $SCRIPTS/08_1_ncrna.py \
	"$coding_gff" "$barrnap_gff" "$trna_out" "$trna_hc_out" "$trna_gff" "$trna_ss" \
	"$RRNA_EXCLUDE_REGION" "$EXCLUDE_NCRNA" "$final_gff" "$overlap_report"
conda deactivate

#####################################################
# Step 5: check that no protein-coding feature was lost
#####################################################

echo "=== Step 5: feature counts ==="

count() { awk -F'\t' -v t="$2" '$3==t' "$1" | wc -l; }
for t in gene mRNA exon CDS five_prime_UTR three_prime_UTR; do
	a=$(count $genome_gff $t); b=$(count $final_gff $t)
	if [ $t = gene ]; then b=$((b - $(count $final_gff rRNA) - $(count $final_gff tRNA))); fi
	if [ $t = exon ]; then b=$((b - $(awk -F'\t' '$3=="exon" && $9 ~ /\.(rrna|trna)/' $final_gff | wc -l))); fi
	echo "$t: projected $a, final $b"
	if [ "$a" -ne "$b" ]; then
		echo "ERROR: number of protein-coding $t features changed" >&2
		exit 1
	fi
done
echo "rRNA: $(count $final_gff rRNA), tRNA: $(count $final_gff tRNA)"

echo "=== Done ==="
echo "Complete GFF (protein-coding + rRNA + tRNA): $final_gff"
echo "Review $report and $overlap_report, then run 08_2_mito_gff.py and 08_3_ena_embl.sh"
