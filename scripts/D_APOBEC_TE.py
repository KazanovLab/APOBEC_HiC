#!/usr/bin/env python
# coding: utf-8

# In[ ]:

import sys
import os

import pandas as pd
import numpy as np
import gzip
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio import SeqIO 
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Seq import MutableSeq

# input data
genome_path = '/mnt/lustre/alice/data/APOBEC_mutagenesis/hg38.fa.gz'
TE_path = '/mnt/lustre/alice/data/APOBEC_mutagenesis/transposons_coordinates.csv'

cancer = ['BLCA', 'BRCA', 'CESC', 'HNSC', 'LUAD', 'LUSC']
cancer_path = '/mnt/lustre/alice/data/APOBEC_mutagenesis/APOBEC_mutagenesis_cancers/'
outfile_path = '/mnt/lustre/alice/data/APOBEC_mutagenesis/D_APOBEC_TE.csv'

# number of the potential APOBEC targets in the genome
N_motifs_genome = 350769311


# read genome file
chr_list = list(range(1,23))
hg38_genome = dict()
with gzip.open(genome_path, "rt") as fasta_file:
    for sequence in SimpleFastaParser(fasta_file):
        chrom = sequence[0][3:]
        if chrom in [str(chrom) for chrom in chr_list]:
            hg38_genome[int(chrom)] = sequence[1].upper()
            print(chrom, end=', ')
        elif chrom in ['X', 'Y']:
            hg38_genome[chrom] = sequence[1].upper()
            print(chrom, end=', ')

# read file with transposon coordinates
transposons_coordinates = pd.read_csv(TE_path, sep='\t')


# count potential APOBEC targets (TC) in TE
transposons_coordinates['N_motifs'] = 0
print(len(transposons_coordinates))
for row in range(0, len(transposons_coordinates)):
    N_chr = transposons_coordinates.iloc[row]['genoName'][3:]
    if N_chr != 'X' and N_chr != 'Y':
        N_chr = int(N_chr)
    start = int(transposons_coordinates.iloc[row]['genoStart'])
    end = int(transposons_coordinates.iloc[row]['genoEnd'])
    seq = hg38_genome[N_chr][start:end]
    transposons_coordinates.at[row, 'N_motifs'] = seq.count('TC') + seq.count('GA')
N_motifs_TE = transposons_coordinates[['repClass', 'N_motifs']]
N_motifs_TE = N_motifs_TE.groupby('repClass').sum()
N_motifs_TE = N_motifs_TE.rename(columns={"N_motifs": "N_motifs_inTE"})
N_motifs_TE['N_motifs_outsideTE'] = N_motifs_genome - N_motifs_TE.N_motifs_inTE
    

# result dataframe with D_APOBEC
col_names = ['cancer','sample', 'transposon_class', 'N_motifs_inTE', 'N_motifs_outsideTE', 'N_mutated_motifs_inTE', 'N_mutated_motifs_outsideTE', 'D_APOBEC_in', 'D_APOBEC_outside']
APOBEC_mutagenesis_TE = pd.DataFrame(columns=col_names)


# count mutated targets in TE
# replace mutated position in the genome to 'N' by cancer type and sample
hg38_genome_mutated = hg38_genome.copy()
for key, value in hg38_genome_mutated.items():
    hg38_genome_mutated[key] = MutableSeq(value)

chrom_list_XY = list(range(1,23)) + ['X', 'Y']
chr_list_XY = ['chr'+str(i) for i in chrom_list_XY]

for cancer_type in cancer:
    print(cancer_type)
    df = pd.read_csv(cancer_path+cancer_type+'_hg38.bed', sep='\t', names=['chr', 'start', 'end', 'Sample'])
    df = df[df['chr'].isin(chr_list_XY)]
    
    for sample in df.Sample.unique():
        df2 = df[df['Sample']==sample]
        df2 = df2.reset_index(drop=True)
        
        mutated_sample_genome = hg38_genome_mutated.copy()
        not_mutated_motifs_count = dict()
    
        for Chr in chrom_list_XY:
            N_chr = 'chr'+str(Chr)
            df3 = df2[df2['chr']==N_chr]
            df3 = df3.reset_index(drop=True)
            mutation_list = df3.start.to_list()
            for ind in mutation_list:
                mutated_sample_genome[Chr][ind] = 'N'
                
            # count the number of NOT mutated motifs in the genome
            not_mutated_motifs_count[Chr] = pd.Series([mutated_sample_genome[Chr].count('TC'), mutated_sample_genome[Chr].count('GA')], index = ['TC', 'GA'])
        not_mutated_motifs_count = pd.concat(not_mutated_motifs_count).unstack()
        not_mutated_motifs_count['TC_GA'] = not_mutated_motifs_count['TC']+not_mutated_motifs_count['GA']
        N_mutated_motifs_genome = N_motifs_genome - not_mutated_motifs_count.TC_GA.sum()

        
        # count the number of NOT mutated motifs in TE for the sample
        transposons_coordinates['N_not_mutated_motifs'] = 0
        for row in range(0, len(transposons_coordinates)):
            N_chr = transposons_coordinates.iloc[row]['genoName'][3:]
            if N_chr != 'X' and N_chr != 'Y':
                N_chr = int(N_chr)
            start = int(transposons_coordinates.iloc[row]['genoStart'])
            end = int(transposons_coordinates.iloc[row]['genoEnd'])
            seq = mutated_sample_genome[N_chr][start:end]
            transposons_coordinates.at[row, 'N_not_mutated_motifs'] = seq.count('TC') + seq.count('GA')

        transposons_coordinates['N_mutated_motifs'] = transposons_coordinates.N_motifs - transposons_coordinates.N_not_mutated_motifs
        
        grouped_TE = transposons_coordinates[['repClass', 'N_mutated_motifs']]
        grouped_TE = grouped_TE.groupby('repClass').sum()
        grouped_TE = grouped_TE.rename(columns={'N_mutated_motifs':'N_mutated_motifs_inTE'})
        grouped_TE['N_mutated_motifs_outsideTE'] = N_mutated_motifs_genome - grouped_TE.N_mutated_motifs_inTE
        
        new_df = pd.concat([N_motifs_TE, grouped_TE], axis=1)
        new_df['D_APOBEC_in'] = new_df.N_mutated_motifs_inTE / new_df.N_motifs_inTE
        new_df['D_APOBEC_outside'] = new_df.N_mutated_motifs_outsideTE / new_df.N_motifs_outsideTE
        new_df = new_df.reset_index()
        new_df = new_df.rename(columns={'repClass':'transposon_class'})
        new_df = new_df[['transposon_class', 'N_motifs_inTE', 'N_motifs_outsideTE', 'N_mutated_motifs_inTE', 'N_mutated_motifs_outsideTE', 'D_APOBEC_in', 'D_APOBEC_outside']]
        for row in range(0, len(new_df)):
            APOBEC_mutagenesis_TE.loc[len(APOBEC_mutagenesis_TE.index)] = [cancer_type, sample]+new_df.iloc[row].to_list()


APOBEC_mutagenesis_TE.to_csv(outfile_path, sep='\t')