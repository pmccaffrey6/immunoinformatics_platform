#!/usr/bin/python

import sys
import pandas as pd


bepipred_inpath = sys.argv[1]
bepipred_df_outpath = sys.argv[2]


with open(bepipred_inpath, 'r') as f:
    rawdat = f.read().split('input: ')[1:]

    dat = [i.split('Position\tResidue')[0] for i in rawdat]
    protname,tabledat = dat[0].split('Predicted peptides\n')

    dfs = []

    for datchunk in dat:
        protname,tabledat = datchunk.split('Predicted peptides\n')
        protname = protname.strip('\n')
        dfrows = [i.split('\t') for i in tabledat.split('\n')]
        df = pd.DataFrame(dfrows[1:], columns=dfrows[0])
        df['protein_full'] = protname
        df['protein_id'] = protname.split(' ')[0]
        df['protein_short'] = protname.split(' ')[0].split('.')[0]
        dfs.append(df)

    all_bepipred_out = pd.concat(dfs)

    all_bepipred_out.to_csv(bepipred_df_outpath, sep='\t', index=False)
