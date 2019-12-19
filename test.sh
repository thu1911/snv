#!/bin/bash
ref_fasta="/data8t/mtx/useful_data/realdata/hg19/hg19.fa"
ref_gtf="/data8t/mtx/useful_data/realdata/hg19/gencode.v19.annotation.gtf"
picard="/home/mtx/software/picard.jar"
gatk="/home/mtx/software/gatk-4.1.4.1/gatk"
ref_snp1="/data8t/mtx/useful_data/realdata/hg19/dbsnp_138.hg19.vcf"
ref_snp2="/data8t/mtx/useful_data/realdata/hg19/1000G_phase1.snps.high_confidence.hg19.sites.vcf"
ref_indel1="/data8t/mtx/useful_data/realdata/hg19/1000G_phase1.indels.hg19.sites.vcf"
ref_indel2="/data8t/mtx/useful_data/realdata/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf"

for i in ../data/fastq/*_1.fastq;
do
    filename=`echo $i |awk -F/ '{print $NF}' |  awk 'gsub("_1.fastq","")'`
    ./SNV_calling.sh $ref_gtf $ref_fasta $filename $picard \
        $gatk $ref_snp1 $ref_snp2 $ref_indel1 $ref_indel2
done
