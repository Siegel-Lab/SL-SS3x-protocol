#!/bin/bash
#SBATCH --ntasks=16                     # number of threads
#SBATCH --mem 128G                      #memory pool for all cores 
#SBATCH -o slurm.%u.%j.out              #STDOUT 
#SBATCH -e slurm.%u.%j.err              #STDERR

#######################################
# Example library mapping and counting
#######################################

# This module aims to generate count matrices for single-cell sequencing data obtained with SL-SS3xpress pipeline. From downloading the data to generating plots.
 

#to run ./run.sh
#
###Requirements##



#conda env create -f SL_SS3xpress_env.yaml

source /work/project/ladsie_003/miniconda3/etc/profile.d/conda.sh # how do i do this for the tutorial?
conda activate SL_SS3xpress_env 

main(){
	set_variables
	create_folders

#	download_reads

	build_metadata

#	make_barcodes_only_list

##	link_genomes

#	merge_bc_umi_reads
#	filter_out_artefact_reads
#	make_fastq_random_subsamples
#	build_ERCC_hybrid_genome_and_annotation
#	star_index_genome
#	star_solo
#	ix_hopping_filtering

}

###############################################################################
# Setting of global variables
###############################################################################

set_variables(){
        INPUT_FOLDER=input
        OUTPUT_FOLDER=output
        SCRIPTS_FOLDER=bin
	INPUT_GENOME_FOLDER=${INPUT_FOLDER}/genome_and_annotation
	INPUT_READS_FOLDER=${INPUT_FOLDER}/reads

	CONFIG_FILE=SS3x_example_config_file.txt
	ALL_TAGS_FASTA=TAGs.fasta
	ASSIGNMENTS_FILE=SS3x_example_assignments.tsv
	METADATA_FILE=SS3x_example_desc.txt
#	METADATA_FILE=test.txt

#	BC_ASSOCIATION_FILE=SS3x_003_desc.txt
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
	TAG_FASTA_FILE=${METADATA_FILE%desc.txt}used_TAGs.fasta
#	TAG_FASTA_FILE=SS3x_example_TAGs_list.fasta
	SL_ERCC_GENOME_FOLDER=${OUTPUT_FOLDER}/SL_ERCC_genome
	STAR_INDEX_FOLDER=${OUTPUT_FOLDER}/star_index
	STAR_MAPPING_FOLDER=${OUTPUT_FOLDER}/star_mapping
}


###############################################################################
# Prelude
###############################################################################

create_folders(){
        mkdir -p \
        $OUTPUT_FOLDER \
	$INPUT_READS_FOLDER \
	$SL_ERCC_GENOME_FOLDER \
	$STAR_INDEX_FOLDER \
	$STAR_MAPPING_FOLDER
	
}


## download the data or put a small enough dataset?
download_reads(){

wget -O ${INPUT_READS_FOLDER}/${READS_FILE} ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR127/ERR12711090/SS3x_003_filtered_lane1_R1.fixed.fastq.gz
wget -O ${INPUT_READS_FOLDER}/${TAG_UMI_READS_FILE} ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR127/ERR12711090/SS3x_003_filtered_lane1_R4.fixed.fastq.gz

wget -O ${INPUT_READS_FOLDER}/${I7_READS_FILE} ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR127/ERR12711091/SS3x_003_filtered_lane1_R2.fixed.fastq.gz
wget -O ${INPUT_READS_FOLDER}/${I5_READS_FILE} ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR127/ERR12711091/SS3x_003_filtered_lane1_R3.fixed.fastq.gz

}

build_metadata(){


python ${SCRIPTS_FOLDER}/SL_SS3x_metadata_table_generation.py \
	--config ${INPUT_FOLDER}/${CONFIG_FILE} \
	--assignments ${INPUT_FOLDER}/${ASSIGNMENTS_FILE} \
	--tag_fasta ${INPUT_FOLDER}/${ALL_TAGS_FASTA} \
	--output ${INPUT_FOLDER}/${METADATA_FILE} \
	--tag_output ${INPUT_FOLDER}/${TAG_FASTA_FILE} \
	--index_sets_folder ${INPUT_FOLDER}



}


make_barcodes_only_list(){

col="BC_TAG"; 

awk -F'\t' -v colname="$col" 'NR==1 {for (i=1; i<=NF; i++) if ($i == colname) col=i; next} {print $col}' ${INPUT_FOLDER}/${METADATA_FILE} > $OUTPUT_FOLDER/${BC_FILE}


}


merge_bc_umi_reads(){
# paste the reads containing i7 , i5 and the one containing the TAG and UMI into one read

paste <(zcat ${INPUT_READS_FOLDER}/$I7_READS_FILE) <(zcat ${INPUT_READS_FOLDER}/$I5_READS_FILE) <(zcat ${INPUT_READS_FOLDER}/$TAG_UMI_READS_FILE) | awk -F '\t' '{ if(NR%4==1||NR%4==3) {print $1} else {print $1 $2 $3} }' | gzip > $OUTPUT_FOLDER/$I7_I5_TAG_UMI_READS_FILE


#jid1=$(sbatch --export=I7_READS_FILE=${INPUT_READS_FOLDER}/${I7_READS_FILE},I5_READS_FILE=${INPUT_READS_FOLDER}/${I5_READS_FILE},TAG_UMI_READS_FILE=${INPUT_READS_FOLDER}/${TAG_UMI_READS_FILE},OUTPUT_FOLDER=${OUTPUT_FOLDER},I7_I5_TAG_UMI_READS_FILE=${I7_I5_TAG_UMI_READS_FILE} ${SCRIPTS_FOLDER}/paste_ss3_bc_fastq_files.sh)
#jid1=$(echo $jid1 | sed "s/Submitted batch job //")


}


filter_out_artefact_reads(){
# Filter out reads that contain the TAG sequence (or its reverse complement) in the cDNA read

cutadapt -j 16 --discard-trimmed -O 10 -b file:${INPUT_FOLDER}/${TAG_FASTA_FILE} \
	-o ${OUTPUT_FOLDER}/${READS_FILE%_R1.fastq.gz}_filtered_R1.fastq.gz \
	-p ${OUTPUT_FOLDER}/${I7_I5_TAG_UMI_READS_FILE%_R2.fastq.gz}_filtered_R2.fastq.gz \
	${INPUT_READS_FOLDER}/${READS_FILE} \
	${OUTPUT_FOLDER}/${I7_I5_TAG_UMI_READS_FILE} \
	> ${OUTPUT_FOLDER}/SS3x_example_cutadapt_summary.txt

#jid2=$(sbatch --job-name=cutadapt --wrap="cutadapt -j 16 --discard-trimmed -O 10 -b file:${INPUT_FOLDER}/${TAG_FASTA_FILE} -o ${OUTPUT_FOLDER}/${READS_FILE%_R1.fastq.gz}_filtered_R1.fastq.gz -p ${OUTPUT_FOLDER}/${I7_I5_TAG_UMI_READS_FILE%_R2.fastq.gz}_filtered_R2.fastq.gz ${INPUT_READS_FOLDER}/${READS_FILE} ${OUTPUT_FOLDER}/${I7_I5_TAG_UMI_READS_FILE} > ${OUTPUT_FOLDER}/SS3x_003_cutadapt_summary.txt")
#jid2=$(echo $jid2 | sed "s/Submitted batch job //")

#== Read fate breakdown ==
#Pairs discarded as trimmed:          7,634,246 (1.4%)
#Pairs written (passing filters):   527,092,681 (98.6%)


}

make_fastq_random_subsamples(){


for SUBSAMPLE in 100000
do
        # to get an average of this amounts of reads per cell, considering that there are 369 cells(15 wells are no-cell controls ), we need to subsample so many reads:
        let "SUBSAMPLED_READS=$SUBSAMPLE * 369"
        SUBSAMPLE_FOLDER=${OUTPUT_FOLDER}/SS3x_example_${SUBSAMPLE%000}K_reads_subsample
        mkdir -p ${SUBSAMPLE_FOLDER}

	seqtk sample -s100 ${OUTPUT_FOLDER}/${I7_I5_TAG_UMI_READS_FILE%_R2.fastq.gz}_filtered_R2.fastq.gz ${SUBSAMPLED_READS} | gzip > ${SUBSAMPLE_FOLDER}/${I7_I5_TAG_UMI_READS_FILE%_R2.fastq.gz}_filtered_${SUBSAMPLE%000}K_R2.fastq.gz
	seqtk sample -s100 ${OUTPUT_FOLDER}/${READS_FILE%_R1.fastq.gz}_filtered_R1.fastq.gz ${SUBSAMPLED_READS} | gzip > ${SUBSAMPLE_FOLDER}/${READS_FILE%_R1.fastq.gz}_filtered_${SUBSAMPLE%000}K_R1.fastq.gz

done



#for SUBSAMPLE in 100000 1000000
#do
#        # to get an average of this amounts of reads per cell, considering 369 cells we need to subsample so many reads:
#        let "SUBSAMPLED_READS=$SUBSAMPLE * 369"
#        SUBSAMPLE_FOLDER=${OUTPUT_FOLDER}/SS3x_003_${SUBSAMPLE%000}K_reads_subsample
#        mkdir -p ${SUBSAMPLE_FOLDER}
#
#        jid3=$(sbatch --job-name=ss3x_sub --mem 250G --wrap="seqtk sample -s100 ${OUTPUT_FOLDER}/${I7_I5_TAG_UMI_READS_FILE%_R2.fastq.gz}_filtered_R2.fastq.gz ${SUBSAMPLED_READS} | gzip > ${SUBSAMPLE_FOLDER}/${I7_I5_TAG_UMI_READS_FILE%_R2.fastq.gz}_filtered_${SUBSAMPLE%000}K_R2.fastq.gz")
#        jid3=$(echo $jid3 | sed "s/Submitted batch job //")
#
#        jid4=$(sbatch --job-name=ss3x_sub --mem 250G --wrap="seqtk sample -s100 ${OUTPUT_FOLDER}/${READS_FILE%_R1.fastq.gz}_filtered_R1.fastq.gz ${SUBSAMPLED_READS} | gzip > ${SUBSAMPLE_FOLDER}/${READS_FILE%_R1.fastq.gz}_filtered_${SUBSAMPLE%000}K_R1.fastq.gz")
#        jid4=$(echo $jid4 | sed "s/Submitted batch job //")
#
#done

}


build_ERCC_hybrid_genome_and_annotation(){

# I need to add the sequences of the 10 SL_ERCC constructs to the genome file
#SL_ERCC.fasta
cat ${INPUT_GENOME_FOLDER}/${GENOME_FILE} ${INPUT_GENOME_FOLDER}/${SL_ERCC_FASTA_FILE} > ${SL_ERCC_GENOME_FOLDER}/${GENOME_FILE%.fasta}_SL_ERCC.fasta


# Make hybrid annotation file:

HYBRID_GFF=${ANNOTATION_FILE%.gff3}_SL_ERCC.gff3

grep -v '^##' ${INPUT_GENOME_FOLDER}/${SL_ERCC_ANNOTATION_FILE} | cat ${INPUT_GENOME_FOLDER}/${ANNOTATION_FILE} - > ${SL_ERCC_GENOME_FOLDER}/${HYBRID_GFF}

# change features "pseudogene" to "gene" and "pseudogenic_transcript" to "mRNA" --> to be able to use "mRNA" as feature and quantify also those genes sometimes "wrongly" annotated as pseudogenes.
HYBRID_MRNA_GFF=${HYBRID_GFF%_SL_ERCC.gff3}_mRNA_SL_ERCC.gff3

sed $'s/\tpseudogene\t/\tgene\t/' ${SL_ERCC_GENOME_FOLDER}/${HYBRID_GFF} | sed $'s/\tpseudogenic_transcript\t/\tmRNA\t/' - > ${SL_ERCC_GENOME_FOLDER}/${HYBRID_MRNA_GFF}

}


star_index_genome(){

STAR --runThreadN 16 --runMode genomeGenerate --genomeDir ${STAR_INDEX_FOLDER} --genomeFastaFiles ${SL_ERCC_GENOME_FOLDER}/${GENOME_FILE%.fasta}_SL_ERCC.fasta

#jid8=$(sbatch --job-name=StarIx --wrap="STAR --runThreadN 16 --runMode genomeGenerate --genomeDir ${STAR_INDEX_FOLDER} --genomeFastaFiles ${SL_ERCC_GENOME_FOLDER}/${GENOME_FILE%.fasta}_SL_ERCC.fasta")
#jid8=$(echo $jid8 | sed "s/Submitted batch job //")

}

star_solo(){
HYBRID_GFF=${ANNOTATION_FILE%.gff3}_SL_ERCC.gff3

HYBRID_MRNA_GFF=${HYBRID_GFF%_SL_ERCC.gff3}_mRNA_SL_ERCC.gff3

OUTPUT_PREFIX=SS3x_example

for SUBSAMPLE in 100000
do
	SUBSAMPLE_FOLDER=${OUTPUT_FOLDER}/SS3x_example_${SUBSAMPLE%000}K_reads_subsample

	OUTPUT_PREFIX=SS3x_example_${SUBSAMPLE%000}K_subsample

	STAR --runThreadN 16 --genomeDir ${STAR_INDEX_FOLDER} \
		--readFilesIn ${SUBSAMPLE_FOLDER}/${READS_FILE%_R1.fastq.gz}_filtered_${SUBSAMPLE%000}K_R1.fastq.gz ${SUBSAMPLE_FOLDER}/${I7_I5_TAG_UMI_READS_FILE%_R2.fastq.gz}_filtered_${SUBSAMPLE%000}K_R2.fastq.gz \
		--readFilesCommand zcat --soloType CB_UMI_Simple --soloCBwhitelist ${OUTPUT_FOLDER}/${BC_FILE} \
		--soloCBstart 1 --soloCBlen 27 --soloUMIstart 28 --soloUMIlen 8 \
		--sjdbGTFfile ${SL_ERCC_GENOME_FOLDER}/${HYBRID_MRNA_GFF} \
		--outFileNamePrefix ${STAR_MAPPING_FOLDER}/${OUTPUT_PREFIX} \
		--outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 10000000000 --outBAMsortingBinsN 200 \
		--outSAMattributes NH HI nM AS CB UB GX GN --outSAMattrRGline ID:rg1 SM:sm1 LB:lib1 PL:illumina PU:unit1 --outSAMmapqUnique 60 \
		--sjdbGTFfeatureExon mRNA --sjdbGTFtagExonParentTranscript ID --sjdbGTFtagExonParentGene Parent  \
		--soloCBmatchWLtype 1MM --soloUMIdedup 1MM_All --soloUMIfiltering MultiGeneUMI --outFilterMultimapNmax 40

# "--limitBAMsortRAM 10000000000" AND "--outBAMsortingBinsN 200" are to avoid it crushing in the cluster because to much RAM usage during the sorting of the BAM.

## IT TOOK 16hs (mostly sorting the bam, counting finished after ~22min)


## DOUBLE-CHECK THIS:
# I don't think i need the "--outSAMattrRGline ID:rg1 SM:sm1 LB:lib1 PL:illumina PU:unit1" nor the --outSAMmapqUnique 60

done

}

ix_hopping_filtering(){

for SUBSAMPLE in 100000
do
	BAM_PREFIX=SS3x_example_${SUBSAMPLE%000}K_subsample
	BAM=${BAM_PREFIX}Aligned.sortedByCoord.out.bam
	OUTPUT_PREFIX=SS3x_example_${SUBSAMPLE%000}K

	python ${SCRIPTS_FOLDER}/compute_neg_matrix.py ${STAR_MAPPING_FOLDER}/${BAM} ${STAR_MAPPING_FOLDER}/${OUTPUT_PREFIX}

done


#for SUBSAMPLE in 100000 1000000
#do
#	BAM_PREFIX=SS3x_003_${SUBSAMPLE%000}K_subsample
#	BAM=${BAM_PREFIX}Aligned.sortedByCoord.out.bam
#	OUTPUT_PREFIX=SS3x_003_${SUBSAMPLE%000}K
#	jid12=$(sbatch --job-name=ix_hop_filt --wrap="python ${SCRIPTS_FOLDER}/compute_neg_matrix.py ${STAR_MAPPING_FOLDER}/${BAM} ${STAR_MAPPING_FOLDER}/${OUTPUT_PREFIX}")
#	jid11=$(echo $jid11 | sed "s/Submitted batch job //")
#done


}


main
