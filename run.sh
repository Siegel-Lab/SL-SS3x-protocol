#!/usr/bin/env bash
#SBATCH --job-name=run_SL_SS3x_pipeline	
#SBATCH --time=2-00:00:00		# job time limit
#SBATCH --cpus-per-task=16              # number of threads
#SBATCH --mem 128G                      #memory pool for all cores 
#SBATCH -o slurm.%u.%j.out              #STDOUT 
#SBATCH -e slurm.%u.%j.err              #STDERR

# make script fail on first error and treat unset vars as error
set -euo pipefail
IFS=$'\n\t'

###########################################################
# Example SL_SS3xpress library mapping and counting pipeline
###########################################################

# This module aims to generate count matrices for single-cell sequencing data obtained with SL-SS3xpress pipeline. From downloading the data to generating plots.
 

# Handle help option
if [[ "${1:-}" == "--help" || "${1:-}" == "-h" ]]; then
    echo "
Usage: ./run.sh [options]

This script runs the SL_SS3xpress pipeline end-to-end.

Options:
  --help       Show this help message and exit

Notes:
  • Designed for SLURM (with #SBATCH headers), but runs locally on macOS/Linux.
  • Requires conda environment: SL_SS3xpress_env
"
    exit 0
fi

#conda env create -f SL_SS3xpress_env.yaml

# activate conda env
if command -v conda >/dev/null 2>&1; then
    eval "$(conda shell.bash hook)"
    conda activate SL_SS3xpress_env
else
    echo "ERROR: conda not found in PATH. Please load conda module or add to PATH." >&2
    exit 1
fi

# check commands
for cmd in conda wget cutadapt seqtk STAR python awk sed grep paste gzip gunzip; do
    command -v "$cmd" >/dev/null 2>&1 || { echo "ERROR: $cmd not found" >&2; exit 1; }
done

main(){
	set_variables
	create_folders
	build_metadata
	make_barcodes_only_list
	merge_bc_umi_reads
	filter_out_artefact_reads
	make_fastq_random_subsamples
	build_ERCC_hybrid_genome_and_annotation
	star_index_genome
	star_solo
	ix_hopping_filtering

}

###############################################################################
# Setting of global variables
###############################################################################

set_variables(){
        INPUT_FOLDER=input
        OUTPUT_FOLDER=output
        SCRIPTS_FOLDER=bin
	INPUT_GENOME_FOLDER="${INPUT_FOLDER}"/genome_and_annotation
	INPUT_READS_FOLDER="${INPUT_FOLDER}"/reads

	CONFIG_FILE=SS3x_example_config_file.txt
	ALL_TAGS_FASTA=TAGs.fasta
	ASSIGNMENTS_FILE=SS3x_example_assignments.tsv
	METADATA_FILE=SS3x_example_desc.txt

	GENOME_FILE=Tb427v11.fasta
	SL_ERCC_FASTA_FILE=SL_ERCC.fasta
	ANNOTATION_FILE=Tb427v11.gff3
	SL_ERCC_ANNOTATION_FILE=SL_ERCC.gff3

	READS_FILE=example_R1.fastq.gz
	TAG_UMI_READS_FILE=example_R2.fastq.gz
	I7_READS_FILE=example_i7.fastq.gz
	I5_READS_FILE=example_i5.fastq.gz
	I7_I5_TAG_UMI_READS_FILE=i7_i5_tag_UMI_R2.fastq.gz

	BC_FILE=BCs.txt
	TAG_FASTA_FILE="${METADATA_FILE%desc.txt}used_TAGs.fasta"

	SL_ERCC_GENOME_FOLDER="${OUTPUT_FOLDER}"/SL_ERCC_genome
	STAR_INDEX_FOLDER="${OUTPUT_FOLDER}"/star_index
	STAR_MAPPING_FOLDER="${OUTPUT_FOLDER}"/star_mapping

	CPUS="${SLURM_CPUS_PER_TASK:-${SLURM_CPUS_ON_NODE:-4}}"

}


###############################################################################
# Prelude
###############################################################################

create_folders(){
mkdir -p \
	"$OUTPUT_FOLDER" \
	"$INPUT_READS_FOLDER" \
	"$SL_ERCC_GENOME_FOLDER" \
	"$STAR_INDEX_FOLDER" \
	"$STAR_MAPPING_FOLDER"
	
}


build_metadata(){


python "${SCRIPTS_FOLDER}"/SL_SS3x_metadata_table_generation.py \
	--config "${INPUT_FOLDER}/${CONFIG_FILE}" \
	--assignments "${INPUT_FOLDER}/${ASSIGNMENTS_FILE}" \
	--tag_fasta "${INPUT_FOLDER}/${ALL_TAGS_FASTA}" \
	--output "${INPUT_FOLDER}/${METADATA_FILE}" \
	--tag_output "${INPUT_FOLDER}/${TAG_FASTA_FILE}" \
	--index_sets_folder "${INPUT_FOLDER}"



}


make_barcodes_only_list(){

col="BC_TAG"; 

# Check if BC_TAG column exists on the metadata file
if ! awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="BC_TAG") {print i; exit}}' "${INPUT_FOLDER}/${METADATA_FILE}" > /dev/null; then
  echo "ERROR: BC_TAG column not found in ${METADATA_FILE}" >&2; exit 1
fi

# Extract the barcode-tags from BC_TAG column
awk -F'\t' -v colname="$col" 'NR==1 {for (i=1; i<=NF; i++) if ($i == colname) col=i; next} {print $col}' "${INPUT_FOLDER}/${METADATA_FILE}" > "$OUTPUT_FOLDER/${BC_FILE}"


}


merge_bc_umi_reads(){
# set local variables

local i7="${INPUT_READS_FOLDER}/${I7_READS_FILE}"
local i5="${INPUT_READS_FOLDER}/${I5_READS_FILE}"
local tagumi="${INPUT_READS_FOLDER}/${TAG_UMI_READS_FILE}"
local out="${OUTPUT_FOLDER}/${I7_I5_TAG_UMI_READS_FILE}"

# Ensure inputs exist
for f in "$i7" "$i5" "$tagumi"
do
  [[ -s "$f" ]] || { echo "Missing input: $f" >&2; exit 1; }
done

# paste the reads containing i7 , i5 and the one containing the TAG and UMI into one read

paste <(gunzip -c "$i7") \
	<(gunzip -c "$i5") \
	<( gunzip -c "$tagumi") \
| awk -F'\t' '{ if(NR%4==1||NR%4==3) {print $1} else {print $1 $2 $3} }' \
| gzip -c > "$out"


}


filter_out_artefact_reads(){
# Filter out reads that contain the TAG sequence (or its reverse complement) in the cDNA read

cutadapt -j "$CPUS" --discard-trimmed -O 10 -b file:"${INPUT_FOLDER}/${TAG_FASTA_FILE}" \
	-o "${OUTPUT_FOLDER}/${READS_FILE%_R1.fastq.gz}_filtered_R1.fastq.gz" \
	-p "${OUTPUT_FOLDER}/${I7_I5_TAG_UMI_READS_FILE%_R2.fastq.gz}_filtered_R2.fastq.gz" \
	"${INPUT_READS_FOLDER}/${READS_FILE}" \
	"${OUTPUT_FOLDER}/${I7_I5_TAG_UMI_READS_FILE}" \
	&> "${OUTPUT_FOLDER}/SS3x_example_cutadapt_summary.txt"

}

make_fastq_random_subsamples(){

local subsamples=(100000)

for SUBSAMPLE in "${subsamples[@]}"
do
        # to get an average of this amounts of reads per cell, considering that there are 369 cells(15 wells are no-cell controls ), we need to subsample so many reads:

	SUBSAMPLED_READS=$(( SUBSAMPLE * 369 ))
        SUBSAMPLE_FOLDER="${OUTPUT_FOLDER}/SS3x_example_${SUBSAMPLE%000}K_reads_subsample"

        mkdir -p "${SUBSAMPLE_FOLDER}"

	seqtk sample -s100 \
	"${OUTPUT_FOLDER}/${I7_I5_TAG_UMI_READS_FILE%_R2.fastq.gz}_filtered_R2.fastq.gz" \
	"${SUBSAMPLED_READS}" \
	| gzip -c \
	> "${SUBSAMPLE_FOLDER}/${I7_I5_TAG_UMI_READS_FILE%_R2.fastq.gz}_filtered_${SUBSAMPLE%000}K_R2.fastq.gz"

	seqtk sample -s100 \
	"${OUTPUT_FOLDER}/${READS_FILE%_R1.fastq.gz}_filtered_R1.fastq.gz" \
	"${SUBSAMPLED_READS}" \
	| gzip -c \
	> "${SUBSAMPLE_FOLDER}/${READS_FILE%_R1.fastq.gz}_filtered_${SUBSAMPLE%000}K_R1.fastq.gz"

done

}


build_ERCC_hybrid_genome_and_annotation(){

# I need to add the sequences of the 10 SL_ERCC constructs to the genome file
#SL_ERCC.fasta
cat "${INPUT_GENOME_FOLDER}/${GENOME_FILE}" \
	"${INPUT_GENOME_FOLDER}/${SL_ERCC_FASTA_FILE}" \
	> "${SL_ERCC_GENOME_FOLDER}/${GENOME_FILE%.fasta}_SL_ERCC.fasta"


# Make hybrid annotation file:

HYBRID_GFF="${ANNOTATION_FILE%.gff3}_SL_ERCC.gff3"

# Create header

echo "##gff-version 3" > "${SL_ERCC_GENOME_FOLDER}/${HYBRID_GFF}" 

# Append main annotations without headers

grep -v '^##' "${INPUT_GENOME_FOLDER}/${ANNOTATION_FILE}" >> "${SL_ERCC_GENOME_FOLDER}/${HYBRID_GFF}"

# Append SL_ERCC features without headers

grep -v '^##' "${INPUT_GENOME_FOLDER}/${SL_ERCC_ANNOTATION_FILE}" >> "${SL_ERCC_GENOME_FOLDER}/${HYBRID_GFF}"

# change features "pseudogene" to "gene" and "pseudogenic_transcript" to "mRNA" --> to be able to use "mRNA" as feature and quantify also those genes sometimes "wrongly" annotated as pseudogenes.

HYBRID_MRNA_GFF="${HYBRID_GFF%_SL_ERCC.gff3}_mRNA_SL_ERCC.gff3"

sed $'s/\tpseudogene\t/\tgene\t/g; s/\tpseudogenic_transcript\t/\tmRNA\t/' "${SL_ERCC_GENOME_FOLDER}/${HYBRID_GFF}" > "${SL_ERCC_GENOME_FOLDER}/${HYBRID_MRNA_GFF}"


}

star_index_genome(){

STAR --runThreadN "$CPUS" --runMode genomeGenerate --genomeSAindexNbases 11 --genomeDir "${STAR_INDEX_FOLDER}" --genomeFastaFiles "${SL_ERCC_GENOME_FOLDER}/${GENOME_FILE%.fasta}_SL_ERCC.fasta"

}

star_solo(){

HYBRID_GFF="${ANNOTATION_FILE%.gff3}_SL_ERCC.gff3"

HYBRID_MRNA_GFF="${HYBRID_GFF%_SL_ERCC.gff3}_mRNA_SL_ERCC.gff3"

OUTPUT_PREFIX=SS3x_example

local subsamples=(100000)

for SUBSAMPLE in "${subsamples[@]}"
do
	SUBSAMPLE_FOLDER="${OUTPUT_FOLDER}/SS3x_example_${SUBSAMPLE%000}K_reads_subsample"

	OUTPUT_PREFIX="SS3x_example_${SUBSAMPLE%000}K_subsample"

	STAR --runThreadN "$CPUS" --genomeDir "${STAR_INDEX_FOLDER}" \
		--readFilesIn "${SUBSAMPLE_FOLDER}/${READS_FILE%_R1.fastq.gz}_filtered_${SUBSAMPLE%000}K_R1.fastq.gz" "${SUBSAMPLE_FOLDER}/${I7_I5_TAG_UMI_READS_FILE%_R2.fastq.gz}_filtered_${SUBSAMPLE%000}K_R2.fastq.gz" \
		--readFilesCommand gunzip -c --soloType CB_UMI_Simple --soloCBwhitelist "${OUTPUT_FOLDER}/${BC_FILE}" \
		--soloCBstart 1 --soloCBlen 27 --soloUMIstart 28 --soloUMIlen 8 \
		--sjdbGTFfile "${SL_ERCC_GENOME_FOLDER}/${HYBRID_MRNA_GFF}" \
		--outFileNamePrefix "${STAR_MAPPING_FOLDER}/${OUTPUT_PREFIX}" \
		--outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 10000000000 --outBAMsortingBinsN 200 \
		--outSAMattributes NH HI nM AS CB UB GX GN --outSAMattrRGline ID:rg1 SM:sm1 LB:lib1 PL:illumina PU:unit1 --outSAMmapqUnique 60 \
		--sjdbGTFfeatureExon mRNA --sjdbGTFtagExonParentTranscript ID --sjdbGTFtagExonParentGene Parent  \
		--soloCBmatchWLtype 1MM --soloUMIdedup 1MM_All --soloUMIfiltering MultiGeneUMI --outFilterMultimapNmax 40 \
		2> "${STAR_MAPPING_FOLDER}/star.${OUTPUT_PREFIX}.err"

done

}

ix_hopping_filtering(){

local subsamples=(100000)

for SUBSAMPLE in "${subsamples[@]}"
do
	BAM_PREFIX="SS3x_example_${SUBSAMPLE%000}K_subsample"
	BAM="${BAM_PREFIX}Aligned.sortedByCoord.out.bam"
	OUTPUT_PREFIX="SS3x_example_${SUBSAMPLE%000}K"

	python "${SCRIPTS_FOLDER}/compute_neg_matrix.py" "${STAR_MAPPING_FOLDER}/${BAM}" "${STAR_MAPPING_FOLDER}/${OUTPUT_PREFIX}"

done

}

start_time=$(date +%s)

main

end_time=$(date +%s)
runtime=$(( end_time - start_time ))

printf "\n✅ Pipeline completed in %02dh:%02dm:%02ds\n" \
  $((runtime/3600)) $(((runtime%3600)/60)) $((runtime%60))

