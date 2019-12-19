#!/bin/bash
location="/data8t/mtx/scSNV/kim_data/data/fastq/"
for i in `find $location -maxdepth 1 -type d`;
do
    if [ -d $i ];then
        echo $i | awk -F/ '{print $NF}'
    fi
done

