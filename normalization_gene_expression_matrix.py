import pandas as pd
import numpy as np
import os
import sys
# arguments
arguments = sys.argv
input_folder = arguments[1]

# the folder of gene expression which contains gene expression matrix and parsed data
gene_expression_folder = input_folder
GENE_MATRIX = "gene_expression_matrix.csv"
FEATURE_PARSE = "feature_parse_matrix.csv"
gene_matrix_file = os.path.join(gene_expression_folder, GENE_MATRIX)
feature_parse_file = os.path.join(gene_expression_folder, FEATURE_PARSE)
# extract gene efficient length from parsed file
feature_parse_df = pd.read_csv(feature_parse_file)
gene_length_df = feature_parse_df[['Geneid', 'Length']]
gene_length_df.set_index('Geneid', inplace=True)
# transform read counts to tpm
gene_matrix_df = pd.read_csv(gene_matrix_file, index_col='Geneid')
cnt_transformed_save_to_local(gene_matrix_df, gene_length_df, gene_expression_folder)

def cnt_transformed_save_to_local(gene_df, length_df, outputfolder):
    _ = cnt2tpm(gene_df, length_df,outputfolder)
    _ = cnt2fpkm(gene_df, length_df,outputfolder)
    _ = cnt2cpm(gene_df,outputfolder)

def cnt2tpm(gene_df, length_df, outputfolder=None):
    merged_df = pd.merge(length_df, gene_df, left_index=True,right_index=True)
    merged_df.iloc[:,1:] = merged_df.iloc[:,1:].div(merged_df['Length'], axis=0)
    tpm_df = merged_df.iloc[:,1:].div(merged_df.sum()[1:])  * (10**6)
    if outputfolder:
        outputfile = os.path.join(outputfolder, "gene_expression_matrix_tpm.csv")
        tpm_df.to_csv(outputfile)
    return tpm_df

def cnt2fpkm(gene_df, length_df,outputfolder=None):
    merged_df = pd.merge(length_df, gene_df, left_index=True,right_index=True)
    fpkm_df = merged_df.iloc[:,1:].div(merged_df['Length'], axis=0).div(merged_df.sum()[1:]) * (10**9)
    if outputfolder:
        outputfile = os.path.join(outputfolder, "gene_expression_matrix_fpkm.csv")
        fpkm_df.to_csv(outputfile)
    return fpkm_df

def cnt2cpm(gene_df,outputfolder=None):
    cpm_df = gene_df.div(gene_df.sum()) * (10**6)
    if outputfolder:
        outputfile = os.path.join(outputfolder, "gene_expression_matrix_cpm.csv")
        cpm_df.to_csv(outputfile)
    return cpm_df