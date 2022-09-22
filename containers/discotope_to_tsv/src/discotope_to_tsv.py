import pandas as pd
from Bio.PDB.Polypeptide import three_to_one
import numpy as np
import sys

discotope_infile = sys.argv[1]
discotope_tsv_outfile = sys.argv[2]
protein_name = sys.argv[3]

discotope_df = pd.read_csv(discotope_infile, sep='\t', skiprows=4, names=['chain_id','residue_id','residue_name','contact_number','propensity_score','discotope_score','epitope_flag'])
discotope_df = discotope_df[pd.notna(discotope_df['residue_id'])]


prelim_epitope_residues = discotope_df[discotope_df['epitope_flag']=='<=B'].sort_values(by='residue_id').index.values
offset_array = np.append(discotope_df[discotope_df['epitope_flag']=='<=B'].sort_values(by='residue_id').index.values[1:],
         (discotope_df[discotope_df['epitope_flag']=='<=B'].sort_values(by='residue_id').index.values[-1]+1))
steps = (offset_array - prelim_epitope_residues)

counter = 0
indexes = []
tags = []

for i,v in enumerate(steps):
    try:
        if v == 1:
            indexes.append(counter)
            tags.append(original_array[i])
        elif v > 1:
            counter += 1
    except IndexError:
        pass

discotope_index_df = pd.DataFrame(data={
    'epitope_index':indexes,
    'discotope_rows':tags
})

discotope_pre_coalesce = discotope_df.reset_index().merge(discotope_index_df, left_on='index', right_on='discotope_rows', how='left')
epitope_indices = list(discotope_pre_coalesce['epitope_index'].dropna().unique())



###
#
###
def create_single_row_from_epitope_table(epitope_table):
    mean_discotope_score = [epitope_table['discotope_score'].mean()]
    max_discotope_score = [epitope_table['discotope_score'].max()]
    min_discotope_score = [epitope_table['discotope_score'].min()]
    std_discotope_score = [epitope_table['discotope_score'].std()]

    discotope_sequence = [''.join([three_to_one(i) for i in epitope_table['residue_name'].values])]
    len_discotope_sequence = [len(discotope_sequence[0])]

    out_df = pd.DataFrame(data={
        'sequence':discotope_sequence,
        'length':len_discotope_sequence,
        'mean_score':mean_discotope_score,
        'max_score':max_discotope_score,
        'min_score':min_discotope_score,
        'std_score':std_discotope_score
    })
    return out_df


out_dfs = []
for epitope_index in epitope_indices:
    out_dfs.append(create_single_row_from_epitope_table(discotope_pre_coalesce[discotope_pre_coalesce['epitope_index']==epitope_index]))

discotope_table = pd.concat(out_dfs)
discotope_table['protein_name'] = protein_name

filtered_discotope_table = discotope_table[(discotope_table['length']>8) & (discotope_table['length']<15)]

filtered_discotope_table.to_csv(discotope_tsv_outfile, sep='\t', index=False)
