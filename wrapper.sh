#!/bin/bash
ref_fasta="/data8t/mtx/useful_data/realdata/hg19/hg19.fa"
ref_gtf="/data8t/mtx/useful_data/realdata/hg19/gencode.v19.annotation.gtf"
picard="/home/mtx/software/picard.jar"
gatk="/home/mtx/software/gatk-4.1.4.1/gatk"
ref_snp="/data8t/mtx/useful_data/realdata/hg19/dbsnp_138.hg19.vcf"

# path setting
path_to_fastq="../data/fastq/"
if [ ! -d $path_to_fastq ]; then
    mkdir -p $path_to_fastq
fi
path_to_bam="../data/bam/"
if [ ! -d $path_to_bam ]; then
    mkdir -p $path_to_bam
fi
path_to_time_stats="../data/time_stats/"
if [ ! -d $path_to_time_stats ]; then
    mkdir -p $path_to_time_stats
fi
path_to_star_index="../data/star_index/"

if false; then
# build index for star
./mapping.sh star_index "$ref_gtf" "$ref_fasta"
fi

# cellwise-operation
for i in "$path_to_fastq"*_1.fastq;
do
    filename=`echo $i |awk -F/ '{print $NF}' |  awk 'gsub("_1.fastq","")'` 
    ./cell_level_analysis.sh $ref_gtf $ref_fasta $filename $picard $gatk $ref_snp 
done


# generate statistic
./make_matrix.sh star_statistic
./make_matrix.sh featurecounts_statistic
