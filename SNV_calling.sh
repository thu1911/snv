#!/bin/bash
# this code is for calling snv
# at least 2 parameters
# 1st arg, gtf file
# 2nd arg, fasta file
# 3rd arg, input file name without location
# 4th arg, picard.jar location
# 5th arg, gatk location
# 6th arg(optional), snp reference 1
# 7th arg(optional), snp reference 2
# 8th arg(optional), indel reference 1
# 9th arg(optional), indel reference 2


# argument setting
ref_gtf=$1
ref_fasta=$2
filename=$3
picard=$4
gatk=$5
ref_snp1=${6:SNP_reference_without_given}
ref_snp2=${7:SNP_reference_without_given}
ref_indel1=${8:Indel_reference_without_given}
ref_indel2=${9:Indel_reference_without_given}

# folder path
path_to_bam="../data/bam/"
path_to_snv="../data/snv/"
path_to_cell_level_snv="$path_to_snv"'cell_level_snv/'
if [ ! -d $path_to_cell_level_snv ]; then
    mkdir -p $path_to_cell_level_snv
fi
path_to_time_stats="../data/time_stats/"
if [ ! -d $path_to_time_stats ]; then
    mkdir -p $path_to_time_stats
fi
gatk_time_stats="$path_to_time_stats"'time_GATK.csv'
if [ ! -f $gatk_time_stats ]; then
    echo "Filename, Time" >> $gatk_time_stats
fi


# For GATK, it requires that reference folder contains .fa .fai .dict
# check if .fa, .fai, .dict are in reference folder
path_to_ref_fasta=${ref_fasta%/*}
ref_fasta_fai="$ref_fasta".fai
ref_fasta_dict=`echo $ref_fasta | sed -e 's/.fa/.dict/g'`
if [ ! -f $ref_fasta_fai ]; then
    echo "genome.fa.fai doesn't exist, this code just build a new one."
    samtools faidx "$ref_fasta"
fi
if [ ! -f $ref_fasta_dict ]; then
    echo "genome.dict doesn't exist, this code just build a new one."
    java -jar $picard CreateSequenceDictionary R=$ref_fasta
fi

# Start real work
start_GATK=`date +%s`
# add read group
in_bam="$path_to_bam""$filename"Aligned.sortedByCoord.out.bam
addrg_bam="$path_to_cell_level_snv""$filename"'_sort_rg.bam'
myRGID="$filename"'_RGID'
myRGLB=$filename
myRGPU=$filename
myRGSM=$filename
java -jar $picard AddOrReplaceReadGroups \
    I=$in_bam \
    O=$addrg_bam \
    RGID=$myRGID \
    RGLB=$myRGLB \
    RGPL=illumina \
    RGPU=$myRGPU \
    RGSM=$myRGSM
# markduplicates
dedup_bam="$path_to_cell_level_snv""$filename"'_sort_rg_dedup.bam'
metrics="$path_to_cell_level_snv""$filename"'_metrics.txt'
java -jar $picard MarkDuplicates \
    INPUT=$addrg_bam \
    OUTPUT=$dedup_bam \
    METRICS_FILE=$metrics
# split N CIGAR for RNA
split_bam="$path_to_cell_level_snv""$filename"'_sort_rg_dedup_split.bam'
$gatk SplitNCigarReads -R $ref_fasta -I $dedup_bam -O $split_bam 
# Base Quality Recalibration
recal_table="$path_to_cell_level_snv""$filename"'_recal_table.txt'
$gatk BaseRecalibrator -R $ref_fasta -I $split_bam -O $recal_table \
    --known-sites $ref_snp1 --known-sites $ref_snp2 \
    --known-sites $ref_indel1 --known-sites $ref_indel2 \
    -nct 20
# PrintReads
recal_reads_bam="$path_to_cell_level_snv""$filename"'_sort_rg_dedup_recal.bam'
$gatk PrintReads -R $ref_fasta -I $split_bam -O $recal_reads_bam -nct 20
# ApplyBQSR
gatk_bam="$path_to_cell_level_snv""$filename"'_gatk.bam'
$gatk ApplyBQSR -I $recal_reads_bam -O $gatk_bam --bqsr-recal-file $recal_table
# HaplotypeCaller
raw_variants="$path_to_cell_level_snv""$filename"'_raw_variants.vcf'
$gatk HaplotypeCaller -R $ref_fasta -I $gatk_bam  --dontUseSoftClippedBases \
    -stand_call_conf 20 -O $raw_variants
# VariantFiltration
filtered_variants="$path_to_cell_level_snv""$filename"'_filtered.vcf'
$gatk VariantFiltration -R $ref_fasta -V $raw_variants --filter-name FS -filter "FS > 30.0" --filter-name QD -filter "QD < 2.0" -O $filtered_variants
# Filter out PASS item
filtered_pass_variants="$path_to_cell_level_snv""$filename"'_filtered_pass.vcf'
cat $filtered_variants | grep -e "#\|PASS"  > $filtered_pass_variants
stop_GATK=`date +%s`
echo $filename","$((stop_GATK-start_GATK)) >> $gatk_time_stats

