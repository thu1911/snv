import os
from collections import defaultdict
import pandas as pd
import os
import sys
# arguments
arguments = sys.argv
path_to_gtf = arguments[1] #path to gtf file
path_to_output_folder = arguments[2]
path_to_output = os.path.join(path_to_output_folder,"GTF_index.csv")

gtf_dict = defaultdict(defaultdict)
gtf = open(path_to_gtf,'r')
Lines = gtf.readlines()
gene_id_list = []
chr_list = []
start_list = []
end_list = []
gene_name_list = []

for line in Lines:
    if not line.startswith("#"):
        line_type = line.split('\t')[2]
        pos_chr = line.split('\t')[0]
        pos_start = line.split('\t')[3]
        pos_end = line.split('\t')[4]
        merged_info = line.split('\t')[8]
        gene_id = merged_info.split()[merged_info.split().index('gene_id') + 1].strip(";").strip('"')
        gene_name = merged_info.split()[merged_info.split().index('gene_name') + 1].strip(";").strip('"')
        if line_type == "gene":
            if gene_id not in gene_id_list:
                gene_id_list.append(gene_id)
                chr_list.append(pos_chr)
                start_list.append(pos_start)
                end_list.append(pos_end)
                gene_name_list.append(gene_name)
gtf.close()
gtf_data = {'Geneid':gene_id_list, 'chr':chr_list, 'start':start_list, 'end':end_list, 'gene_name':gene_name_list}
gtf_df = pd.DataFrame(gtf_data)
gtf_df.set_index('Geneid', inplace=True)
gtf_df = gtf_df.astype({'start':'int64', 'end':'int64'})
gtf_df.to_csv(path_to_output)
