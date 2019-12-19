#!/bin/bash
#this script is for generate index and mapping with star software
#take 2 or 3 args
#1st is gtf
#2nd is fasta
#3rd is fastq_name without _1.fastq or _2.fastq

path_to_star_index=../data/star_index/
path_to_bam=../data/bam/
path_to_time_stats=../data/time_stats/
path_to_fastq=../data/fastq/
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
    STAR --runThreadN 20 --genomeDir "$path_to_star_index" \
        --readFilesIn $fastq1 $fastq2 \
        --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix $output \
        --quantMode TranscriptomeSAM 
     samtools index -@ 20 ""$output"Aligned.sortedByCoord.out.bam"

    # this part is for quantmode, we need to sort and index it
     samtools sort -@ 20 ""$output"Aligned.toTranscriptome.out.bam" \
         -T ""$output"tmp" \
         -o ""$output"Aligned.toTranscriptome.sortedByCoord.out.bam"
     samtools index -@ 20 ""$output"Aligned.toTranscriptome.sortedByCoord.out.bam"
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
    samtools index -@ 20 $filename
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
"$@"
