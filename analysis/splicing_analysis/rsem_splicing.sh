#!/bin/bash
#this script is for generate index and mapping with star software
#take 2 or 3 args
#1st is gtf
#2nd is fasta
#3rd is fastq_name without _1.fastq or _2.fastq

path_to_gtf=/data8t/mtx/useful_data/realdata/hg19/gencode.v19.annotation.gtf
path_to_fasta=/data8t/mtx/useful_data/realdata/hg19/hg19.fa
path_to_resm_ref=/data8t/mtx/scSNV/just_for_test/code/analysis/splicing_analysis/data/rsem_ref/rsem_hg19_ref
path_to_result=/data8t/mtx/scSNV/just_for_test/code/analysis/splicing_analysis/data/ZC3H11A/ZC3H11A
path_to_extracted_sam=/data8t/mtx/scSNV/just_for_test/code/analysis/splicing_analysis/data/ZC3H11A/ZC3H11A_extracted.sam
path_to_extracted_sam_no_dup=/data8t/mtx/scSNV/just_for_test/code/analysis/splicing_analysis/data/ZC3H11A/ZC3H11A_extracted_no_dup.sam
path_to_fastq_1=/data8t/mtx/scSNV/just_for_test/code/analysis/splicing_analysis/data/ZC3H11A/ZC3H11A_1.fastq
path_to_fastq_2=/data8t/mtx/scSNV/just_for_test/code/analysis/splicing_analysis/data/ZC3H11A/ZC3H11A_2.fastq

if false; then

# build index for resm
rsem-prepare-reference -p 16 --gtf ${path_to_gtf} --star ${path_to_fasta} ${path_to_resm_ref}

# extract reads from target

for reads in $(cat $path_to_extracted_sam | awk  '{print $1}' | sort | uniq -c |gawk '$1==2{print $2}');
do
    grep ${reads} ${path_to_extracted_sam} >> ${path_to_extracted_sam_no_dup}
done
cat ${path_to_extracted_sam_no_dup} | grep -v ^@ | awk 'NR%2==1 {print "@"$1"\n"$10"\n+\n"$11}' > ${path_to_fastq_1}
cat ${path_to_extracted_sam_no_dup} | grep -v ^@ | awk 'NR%2==0 {print "@"$1"\n"$10"\n+\n"$11}' > ${path_to_fastq_2}

fi

# quantification
rsem-calculate-expression --paired-end --star --paired-end -p 20 --output-genome-bam ${path_to_fastq_1} ${path_to_fastq_2} ${path_to_resm_ref} ${path_to_result}








star_index(){
    # for generating index for star
    # 1st parameter is gtf
    # 2nd parameter is fasta
    mkdir -p "$path_to_star_index"
    STAR --runThreadN 30  --runMode genomeGenerate \
        --genomeDir $path_to_star_index \
        --genomeFastaFiles $2 \
        --sjdbGTFfile $1 \
        --sjdbOverhang 100  
}


star_mapping(){
    # one arg
    # we need filename (without location) without _1.fastq or _1.fastq as input 
    # for example, "ESS12345"
    filename=$1
    fastq1="$path_to_fastq""$filename"_1.fastq
    fastq2="$path_to_fastq""$filename"_2.fastq
    output="$path_to_bam""$filename"
    start_STAR=`date +%s`
    # the threads of star mapping can be set as 20
    STAR --runThreadN 10 --genomeDir "$path_to_star_index" \
        --readFilesIn $fastq1 $fastq2 \
        --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix $output \
        --quantMode TranscriptomeSAM 
     samtools index -@ 10 ""$output"Aligned.sortedByCoord.out.bam"

    # this part is for quantmode, we need to sort and index it
     samtools sort -@ 10 ""$output"Aligned.toTranscriptome.out.bam" \
         -T ""$output"tmp" \
         -o ""$output"Aligned.toTranscriptome.sortedByCoord.out.bam"
     samtools index -@ 10 ""$output"Aligned.toTranscriptome.sortedByCoord.out.bam"
    stop_STAR=`date +%s`
    time_file="$path_to_time_stats"time_STAR.csv
    # check if time file exist. if not, add header
    if [ ! -f $time_file ]; then
        echo 'filename,time' >> $time_file
    fi
    echo $filename","$((stop_STAR-start_STAR)) >> $time_file 
}


samtools_index_bam(){
    # one arg
    # we need BAM filename (without locaton) 
    filename=$1
    samtools index -@ 10 $filename
}
    
samtools_sort_bam(){
    # 3 args
    # arg1: input filename
    # arg2: prefix of tmp file
    # arg3: output 
    # we need BAM filename (without location)
    filename=$1
    samtools sort -@ 20 $filename
}
