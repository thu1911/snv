#!/bin/bash
star_statistic(){
    # this function is for calculate the metrics in log.final.out
    path_to_output='../data/mapping_stat/'
    mkdir -p  $path_to_output
    outputfile="$path_to_output"'mapping_stat.csv'
    if [ -f $outputfile ]; then
        echo "rewrite a new CSV file!"
        rm $outputfile
    fi
    echo "Filename,StartTime,FinishTime,ReadLength,NumReads,Unique,Unique_Percent,Multi,Multi_Percent" >> $outputfile
    path_to_bam="../data/bam/"
    for i in "$path_to_bam"*Log.final.out;
    do
        Filename=`echo $i | awk -F/ '{print $NF}'|awk '{gsub("Log.final.out", "");print}'`
        StartTime=`cat $i | grep "Started mapping on" | awk -F"|" '{print $NF}'`
        FinishTime=`cat $i | grep "Finished on" | awk -F"|" '{print $NF}'`
        ReadLength=`cat $i | grep "Number of input reads" | awk -F"|" '{print $NF}'`
        NumReads=`cat $i | grep "Average input read length" | awk -F"|" '{print $NF}'` 
        Unique=`cat $i | grep "Uniquely mapped reads number" | awk -F"|" '{print $NF}'` 
        Unique_Percent=`cat $i | grep "Uniquely mapped reads %" | awk -F"|" '{print $NF}'` 
        Multi=`cat $i | grep "Number of reads mapped to multiple loci" | awk -F"|" '{print $NF}'` 
        Multi_Percent=`cat $i | grep "% of reads mapped to multiple loci" | awk -F"|" '{print $NF}'` 
        echo $Filename","$StartTime","$FinishTime","$ReadLength","$NumReads","$Unique","$Unique_Percent","$Multi","\
            $Multi_Percent>> $outputfile
    done
}
featurecounts_statistic(){
    # this function is for generating gene expression matrix and summary matrix
    path_to_input="../data/quantification_stats/cell_level_files/"
    path_to_output="../data/quantification_stats/gene_expression_matrix/"
    outputfile_gene_expression="$path_to_output"gene_expression_matrix.csv
    outputfile_featurecounts_summary="$path_to_output"featurecounts_summary_matrix.csv
    outputfile_gtf_parse="$path_to_output"feature_parse_matrix.csv
    # to check if input folder exist
    if [ ! -d "$path_to_input" ]; then
        echo 'Folder '"$path_to_output"' does not exist! Can not make gene expression profile.'
        exit 1
    fi
    # to check if output folder exist
    if [ ! -d "$path_to_output" ]; then
        mkdir -p $path_to_output
    fi
    # parse gtf files
    # I need to remove header(first line) and read counts(last column)
    for i in "$path_to_input"*.featurecounts ;
    do
        NumColumn=`cat $i | tail -n +2 | awk '{print NF}' | sort -nu | head -n 1`
        cat $i | tail -n +2 |tr '\t' ',' |cut --complement -f $NumColumn -d, > "$outputfile_gtf_parse"
        break
    done
    # generate gene expression matrix
    # use python to do process
    # 1st parameter is input folder
    # 2nd parameter is output file name with location
    python ./generate_gene_expression_matrix.py "$path_to_input" "$outputfile_gene_expression"
    python ./normalization_gene_expression_matrix.py "$path_to_output"
    # feature summary matrix
    if [ ! -f "$outputfile_featurecounts_summary" ]; then
        echo 'no feature summary CSV. Build a new one!'
        echo 'Filename,Assigned,Unassigned_Unmapped,Unassigned_Read_Type,Unassigned_Singleton,Unassigned_MappingQuality,Unassigned_Chimera,Unassigned_FragmentLength,Unassigned_Duplicate,Unassigned_MultiMapping,Unassigned_Secondary,Unassigned_NonSplit,Unassigned_NoFeatures,Unassigned_Overlapping_Length,Unassigned_Ambiguity' >> "$outputfile_featurecounts_summary"
    fi
    for i in "$path_to_input"*.featurecounts.summary ;
    do
        Filename=`echo $i | awk -F/ '{print $NF}'|awk '{gsub(".featurecounts.summary", "");print}'`
        Assigned=`cat $i | grep "Assigned" | awk '{print $NF}'`
        Unassigned_Unmapped=`cat $i | grep "Unassigned_Unmapped" | awk '{print $NF}'`
        Unassigned_Read_Type=`cat $i | grep "Unassigned_Read_Type" | awk '{print $NF}'`
        Unassigned_Singleton=`cat $i | grep "Unassigned_Singleton" | awk '{print $NF}'`
        Unassigned_MappingQuality=`cat $i | grep "Unassigned_MappingQuality" | awk '{print $NF}'`
        Unassigned_Chimera=`cat $i | grep "Unassigned_Chimera" | awk '{print $NF}'`
        Unassigned_FragmentLength=`cat $i | grep "Unassigned_FragmentLength" | awk '{print $NF}'`
        Unassigned_Duplicate=`cat $i | grep "Unassigned_Duplicate" | awk '{print $NF}'`
        Unassigned_MultiMapping=`cat $i | grep "Unassigned_MultiMapping" | awk '{print $NF}'`
        Unassigned_Secondary=`cat $i | grep "Unassigned_Secondary" | awk '{print $NF}'`
        Unassigned_NonSplit=`cat $i | grep "Unassigned_NonSplit" | awk '{print $NF}'`
        Unassigned_NoFeatures=`cat $i | grep "Unassigned_NoFeatures" | awk '{print $NF}'`
        Unassigned_Overlapping_Length=`cat $i | grep "Unassigned_Overlapping_Length" | awk '{print $NF}'`
        Unassigned_Ambiguity=`cat $i | grep "Unassigned_Ambiguity" | awk '{print $NF}'`
        echo $Filename","$Assigned","$Unassigned_Unmapped","$Unassigned_Read_Type","$Unassigned_Singleton","$Unassigned_MappingQuality","$Unassigned_Chimera","$Unassigned_FragmentLength","$Unassigned_Duplicate","$Unassigned_MultiMapping","$Unassigned_Secondary","$Unassigned_NonSplit","$Unassigned_NoFeatures","$Unassigned_Overlapping_Length","$Unassigned_Ambiguity >> "$outputfile_featurecounts_summary"
    done

}

"$@"
