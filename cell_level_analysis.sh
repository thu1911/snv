#!/bin/bash
#for cell level analysis
#3 parameter: 
#1st parameter: gtf
#2nd parameter: fasta
#3rd parameter: cell name without location
ref_gtf=$1
ref_fasta=$2
filename=$3
picard=$4
gatk=$5
ref_snp=${6:SNP_reference_without_given}
# mapping
./mapping.sh star_mapping $filename

# quantify
./quantify.sh featurecounts $ref_gtf $filename 
./SNV_calling.sh $ref_gtf $ref_fasta $filename $picard $gatk $ref_snp

