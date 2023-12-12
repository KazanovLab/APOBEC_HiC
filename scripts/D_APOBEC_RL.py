#!/usr/bin/env python
# coding: utf-8
#searching for which position are in R-loops

# In[ ]:

import sys
import os

import pandas as pd
import numpy as np

def closest_value(input_list, input_value):
    arr = np.asarray(input_list)
    i = (np.abs(arr - input_value)).argmin()
    return i

cancer = ['BLCA', 'BRCA', 'CESC', 'HNSC', 'LUAD', 'LUSC']

RL_regions = pd.read_csv('/mnt/lustre/alice/data/RL_loops_APOBEC_mutagenesis/RL_regions.csv')

snv_mnv = dict()
chrom_list_XY = list(range(1,23)) + ['X', 'Y']
chrom_list_XY = ['chr'+str(i) for i in chrom_list_XY]

for cancer_type in cancer:
    print(cancer_type, end=', ')
    df = pd.read_csv('/mnt/lustre/alice/data/RL_loops_APOBEC_mutagenesis/APOBEC_mutagenesis_cancers/'+cancer_type+'_hg38.bed', sep='\t', names=['chr', 'start', 'end', 'sample'])
    df = df[df['chr'].isin(chrom_list_XY)]
    df['inRloop'] = 0
    
    print(len(df))
    
    df = df.reset_index(drop=True)
    for row in range(0, len(df)):
        
        #checking
        if row % 100000 == 0:
            print('working')
        
        Chr = df.iloc[row]['chr']
        position = df.iloc[row]['start']
        
        RL_regions_chr = RL_regions[RL_regions["chrom"]==Chr]
        RL_regions_chr = RL_regions_chr.reset_index(drop=True)
        
        # searching of row's number with the closest to position value
        if __name__ == "__main__" :
            N = closest_value(RL_regions_chr.start.tolist(),position)
        
        # position may be in this row (N) or in previous row (N-1)
        RL_coordinates = [[RL_regions_chr.iloc[N-1]['start'], RL_regions_chr.iloc[N-1]['end']], [RL_regions_chr.iloc[N]['start'], RL_regions_chr.iloc[N]['end']]]
        
        for x in RL_coordinates:
            if position <= x[1] and position >= x[0]:
                df.at[row, 'inRloop'] = 1
                break
    snv_mnv[cancer_type] = df
    
merged_df = pd.concat(snv_mnv.values(), join='inner', keys=snv_mnv.keys())
merged_df.to_csv('/mnt/lustre/alice/data/RL_loops_APOBEC_mutagenesis/snv_mnv_RL.csv', sep='\t')



# Calculation the mutation density D_APOBEC of the APOBEC mutagenesis 

# Number of potential APOBEC targets in R-loops
N_TCN_RL = 16990012

# Number of potential APOBEC targets outside R-loops
N_TCN_out = 333779299

APOBEC_mutagenesis = pd.DataFrame(columns=['cancer','sample','APOBEC_mutagenesis_inRloops','APOBEC_mutagenesis_outsideRloops'])
for cancer_type in cancer:
    df = snv_mnv[cancer_type]
    for sample in df['sample'].unique():
        df_sample = df[df['sample']==sample]

        # number of mutated targets inside R loops
        N_APOBEC_RL = df_sample.inRloop.sum()
        
        # number of mutated targets outside R loops
        N_APOBEC_out = len(df_sample) - N_APOBEC_RL
        
        D_APOBEC_RL = N_APOBEC_RL / N_TCN_RL
        D_APOBEC_out = N_APOBEC_out / N_TCN_out
        
        APOBEC_mutagenesis.loc[len(APOBEC_mutagenesis.index)] = [cancer_type, sample, D_APOBEC_RL, D_APOBEC_out]

APOBEC_mutagenesis.to_csv('/mnt/lustre/alice/data/RL_loops_APOBEC_mutagenesis/D_APOBEC.csv', sep='\t')