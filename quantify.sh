#!/bin/bash
featurecounts(){
    # this is a cell-level operation
    # 2 parameter, 1st parameter is gtf, 2nd parameter is bamfile name without location
    ref_gtf=$1
    bamfilename=$2
    path_to_time_stats="../data/time_stats/"
    path_to_bam_file="../data/bam/"
    subfix='Aligned.sortedByCoord.out.bam'
    bamfile="$path_to_bam_file""$bamfilename""$subfix"
    threads=30
    path_to_quantification_stats="../data/quantification_stats/cell_level_files/"
    # check if path of output of featurecounts exist
    if [ ! -d "$path_to_quantification_stats" ]; then
        mkdir -p $path_to_quantification_stats
    fi
    # check if path of time stats of featurecounts exist
    if [ ! -d "$path_to_time_stats" ]; then
        echo "time folder doesn't exist, buid it first."
        mkdir -p $path_to_time_stats
    fi
    output="$path_to_quantification_stats""$bamfilename"'.featurecounts'
    start_featurecounts=`date +%s`
    featureCounts -T $threads -p -a $ref_gtf -o $output $bamfile 
    stop_featurecounts=`date +%s`
    time_file="$path_to_time_stats"time_featurecounts.csv
    # check if time file of feature counts exist, if not, add header
    if [ ! -f "$time_file" ]; then
        echo 'filename,time' >> $time_file
    fi
    echo $bamfilename","$((stop_featurecounts-start_featurecounts)) >> $time_file
}

"$@"
