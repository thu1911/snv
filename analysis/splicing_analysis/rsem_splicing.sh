#!/bin/bash
#this script is for generate index and mapping with star software
#take 2 or 3 args
#1st is gtf
#2nd is fasta
#3rd is fastq_name without _1.fastq or _2.fastq

####################
path_to_gtf=/data8t/mtx/useful_data/realdata/hg19/gencode.v19.annotation.gtf
path_to_fasta=/data8t/mtx/useful_data/realdata/hg19/hg19.fa
path_to_fastq=/data8t/mtx/disk_out/BACKUP/scSNV/GSE57872/fastq
path_to_bam=/data8t/mtx/disk_out/BACKUP/scSNV/GSE57872/bam
path_to_resm_ref=/data8t/mtx/scSNV/just_for_test/code/analysis/splicing_analysis/data/rsem_ref/rsem_hg19_ref
path_to_star_index=/data8t/mtx/scSNV/just_for_test/code/analysis/splicing_analysis/data/star_index
####################


path_to_filename=/data8t/mtx/scSNV/just_for_test/code/analysis/splicing_analysis/data/TAF7_140698400/chr5__140698400__no__filename.txt
path_to_fastq_1=/data8t/mtx/scSNV/just_for_test/code/analysis/splicing_analysis/data/TAF7_140698400/TAF7_140698400_no_1.fastq
path_to_fastq_2=/data8t/mtx/scSNV/just_for_test/code/analysis/splicing_analysis/data/TAF7_140698400/TAF7_140698400_no_2.fastq
path_to_result=/data8t/mtx/scSNV/just_for_test/code/analysis/splicing_analysis/data/TAF7_140698400/rsem/
path_to_extracted_sam=/data8t/mtx/scSNV/just_for_test/code/analysis/splicing_analysis/data/TAF7_140698400/TAF7_140698400_extracted_no.sam
path_to_extracted_bam=/data8t/mtx/scSNV/just_for_test/code/analysis/splicing_analysis/data/TAF7_140698400/TAF7_140698400_extracted_no.bam
path_to_extracted_bam_sorted=/data8t/mtx/scSNV/just_for_test/code/analysis/splicing_analysis/data/TAF7_140698400/TAF7_140698400_extracted_no.sorted.bam


# mapping path
path_to_mapping_result=/data8t/mtx/scSNV/just_for_test/code/analysis/splicing_analysis/data/TAF7_140698400/mapping
path_to_mapping=${path_to_mapping_result}/TAF7_140698400_no

mkdir -p ${path_to_mapping_result}
mkdir -p ${path_to_result}

if false; then
while IFS= read -r line
do
  # chr5:140697000-140702000
  input_bam=${path_to_bam}/${line}Aligned.sortedByCoord.out.bam
  samtools view ${input_bam} chr5:140697000-140702000  >> ${path_to_extracted_sam}
  #echo ${input_fastq_1}
  #echo ${input_fastq_2}
  echo ${input_bam}
done < $path_to_filename
fi


#if false; then
# sort by name
samtools view -bS ${path_to_extracted_sam} > ${path_to_extracted_bam}

samtools sort -n ${path_to_extracted_bam} >> ${path_to_extracted_bam_sorted}

# bam to fastq
bedtools bamtofastq -i ${path_to_extracted_bam_sorted} -fq $path_to_fastq_1 -fq2 $path_to_fastq_2

# mapping

STAR --runThreadN 10 --genomeDir $path_to_star_index \
        --readFilesIn $path_to_fastq_1 $path_to_fastq_2 \
        --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix ${path_to_mapping} \
        --quantMode TranscriptomeSAM 

# quantification
rsem_input=${path_to_mapping}Aligned.toTranscriptome.out.bam
path_to_rsem=${path_to_result}/TAF7_140698400_no_rsem
rsem-calculate-expression --bam --paired-end -p 16 ${rsem_input} ${path_to_resm_ref} ${path_to_rsem}
#fi
#fi