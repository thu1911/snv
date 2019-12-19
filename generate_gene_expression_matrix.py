#!/usr/bin/python
# this function is for generate gene expression matrix 
# 2 parameters or 3 parameters
# 1st parameter, input folder, which contains cell level statistic
# 2nd parameter, folder, output filename

import os
import sys
from os import listdir
from os.path import isfile, join
import pandas as pd
import functools 

# arguments
arguments = sys.argv
input_folder = arguments[1]
output_filename = arguments[2]

suffix=".featurecounts"
onlyfiles = [f for f in listdir(input_folder) if isfile(join(input_folder, f)) if f.endswith(suffix)]
all_cell_df_list = []
for f in onlyfiles:
    cell_level_file = join(input_folder, f)
    tmp = pd.read_csv(cell_level_file,delimiter='\t',skiprows=[0], header=0)
    # rename last column
    tmp_colname = tmp.columns[-1]
    tmp_rename = tmp.rename({tmp_colname: f.split('.')[0]}, axis='columns')
    # only contain geneid and read counts
    current_cell_df =  tmp_rename.iloc[:,[0,-1]]
    del tmp, tmp_rename
    all_cell_df_list.append(current_cell_df)

gene_expression_matrix = functools.reduce(lambda x, y: pd.merge(x, y, on = 'Geneid'), all_cell_df_list)
gene_expression_matrix.to_csv(output_filename,index=False) 
