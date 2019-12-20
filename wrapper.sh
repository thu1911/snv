#!/bin/bash
# one argument, organism, human or mouse

# set reference 
organism=$1
if [ "$#" -ne 1 ]; then
    echo "Please check your argument. It should be "human" or "mouse""
    exit 1
fi
if [ $organism == 'mouse' ]; then
    ref_fasta="/data8t/mtx/useful_data/realdata/mm10/mm10.fa"
    ref_gtf="/data8t/mtx/useful_data/realdata/mm10/gencode.vM23.annotation.gtf"
    ref_snp="/data8t/mtx/useful_data/realdata/mm10/mgp.v3.snps.rsIDdbSNPv137.mm10.sorted.vcf"
    ref_indel="/data8t/mtx/useful_data/realdata/mm10/mgp.v3.indels.rsIDdbSNPv137.mm10.sorted.vcf"
elif [ $organism == 'human' ]; then
    ref_fasta="/data8t/mtx/useful_data/realdata/hg19/hg19.fa"
    ref_gtf="/data8t/mtx/useful_data/realdata/hg19/gencode.v19.annotation.gtf"
    ref_snp="/data8t/mtx/useful_data/realdata/hg19/dbsnp_138.hg19.vcf"
    ref_indel="/data8t/mtx/useful_data/realdata/hg19/1000G_phase1.indels.hg19.sites.vcf"
else
    echo "Please check your argument. It should be "human" or "mouse""
    exit 1
fi
picard="/home/mtx/software/picard.jar"
gatk="/home/mtx/software/gatk-4.1.4.1/gatk"

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

# build index for star
./mapping.sh star_index "$ref_gtf" "$ref_fasta"

# cellwise-operation
for i in "$path_to_fastq"*_1.fastq;
do
    filename=`echo $i |awk -F/ '{print $NF}' |  awk 'gsub("_1.fastq","")'` 
    ./cell_level_analysis.sh $ref_gtf $ref_fasta $filename $picard $gatk $ref_snp $ref_indel
done


# generate statistic
./make_matrix.sh star_statistic
./make_matrix.sh featurecounts_statistic
