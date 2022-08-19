#!/usr/bin/python

import sys
import pandas as pd


cdhit_inpath = sys.argv[1]         # e.g. 'cdhit_outputs/cdhit_out_sim_095.txt.clstr'
cdhit_df_outpath = sys.argv[2]     # e.g.'cdhit_outputs/cdhit_out_sim_095.clstr.tsv'

def cluster_to_table(cluster):
    lines = cluster.split('\n')[1:]
    index = cluster.split('\n')[0].strip(' ')
    #print(f"index: {index}")

    df_rows = []

    for entry in lines:
        entry_items = entry.split(' ')
        entry_index = entry_items[0].split('\t')[0]
        entry_info = entry_items[0].split('\t')[1]
        entity = entry_items[1].strip('...')

        if entry_items[-1] == '*':
            is_centroid = True
            similarity = 100.0
        else:
            is_centroid = False
            similarity = entry_items[-1].strip('%')

        df_rows.append([index, entry_index, entity, is_centroid, similarity])

    return pd.DataFrame(df_rows, columns=['cluster_index', 'entry_index', 'entity', 'is_centroid', 'similarity'])


with open(cdhit_inpath, 'r') as f:
    raw_data = f.read()
    raw_clusters = [i.strip(' ').strip('\n') for i in raw_data.split('>Cluster')[1:]]
    raw_clusters = [i for i in raw_clusters if i != '']

    cdhit_df = pd.concat([cluster_to_table(cluster) for cluster in raw_clusters])

    cdhit_df.to_csv(cdhit_df_outpath, sep='\t', index=False)
