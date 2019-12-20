#!/bin/bash
ref_snp="/data8t/mtx/useful_data/realdata/mm10/mgp.v3.snps.rsIDdbSNPv137.vcf"
ref_indel="/data8t/mtx/useful_data/realdata/mm10/mgp.v3.indels.rsIDdbSNPv137.vcf"
picard="/home/mtx/software/picard.jar"
gatk="/home/mtx/software/gatk-4.1.4.1/gatk"
### snp
# convert from GRCm38 to mm10 format
ref_snp_mm10=`echo "${ref_snp/.vcf/.mm10.vcf}" `
cat $ref_snp | awk '{if($1 ~ "^#") {gsub("contig=<ID=", "contig=<ID=chr"); \
    gsub("contig=<ID=chrMT", "contig=<ID=chrM"); print} else {gsub("^MT", "M"); \
    print "chr"$0}}' > $ref_snp_mm10
# sort vcf
ref_snp_mm10_sorted=`echo "${ref_snp/.vcf/.mm10.sorted.vcf}" `
java -jar $picard SortVcf I=$ref_snp_mm10 O=$ref_snp_mm10_sorted

### indel
# convert from GRCm38 to mm10 format
ref_indel_mm10=`echo "${ref_indel/.vcf/.mm10.vcf}" `
cat $ref_indel | awk '{if($1 ~ "^#") {gsub("contig=<ID=", "contig=<ID=chr"); \
    gsub("contig=<ID=chrMT", "contig=<ID=chrM"); print} else {gsub("^MT", "M"); \
    print "chr"$0}}' > $ref_indel_mm10
# sort vcf
ref_indel_mm10_sorted=`echo "${ref_indel/.vcf/.mm10.sorted.vcf}" `
java -jar $picard SortVcf I=$ref_indel_mm10 O=$ref_indel_mm10_sorted

