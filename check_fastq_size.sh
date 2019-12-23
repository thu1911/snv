#!/bin/bash
# 1 argument: fastq folder
path_to_fastq=$1
# check if "/" is given
if [[ $path_to_fastq != */ ]]; then
   path_to_fastq="$path_to_fastq"'/'
fi
for i in "$path_to_fastq"*_1.fastq;
do
    file1=$i
    file2=`echo $i |awk 'gsub("_1.fastq","_2.fastq")'`
    size1=$(ls -la $file1 | awk '{ print $5}')
    size2=$(ls -la $file2 | awk '{ print $5}')
    echo $size1
    if [[ $size1 != $size2 ]]; then
        echo $file1
    fi
done

for i in "$path_to_fastq"*_2.fastq;
do
    file1=$i
    file2=`echo $i |awk 'gsub("_2.fastq","_1.fastq")'`
    size1=$(ls -la $file1 | awk '{ print $5}')
    size2=$(ls -la $file2 | awk '{ print $5}')
    echo $size1
    if [[ $size1 != $size2 ]]; then
        echo $file1
    fi
done
echo 'Finish Job'
