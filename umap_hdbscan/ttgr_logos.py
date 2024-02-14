#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import logomaker as lm
from Bio import SeqIO
from collections import Counter
import glob

# Code to generate logos

aa = "ACDEFGHIKLMNPQRSTVWY"
clustered_df = pd.read_csv("avg_foldenrichment_full_clustered.csv")

def get_aa_frequencies(sequences, aminoacids):
    ''' Returns a dataframe where the columns are aminoacids. the rows are positions in the sequences and the values are frequencies of each amino acid at each position.'''
    freq_dict = {char: [] for char in aminoacids}
    for i in range(len(sequences[0])):
        for char in aminoacids:
            count = sum(seq[i] == char for seq in sequences)
            freq_dict[char].append(count)
    freq_df = pd.DataFrame(freq_dict)
    freq_df.index = freq_df.index + 1
    freq_df.columns = [char for char in aminoacids]
    freq_df = freq_df.div(freq_df.sum(axis=1), axis=0)
    return freq_df

## For generating sequence logos for all the clusters
full = get_aa_frequencies(clustered_df['Mut_Pos'].values, aa)

fig, ax = plt.subplots(6, 4, figsize=(30, 20))
plt.subplots_adjust(wspace=0.3, hspace=0.4)
#fig.delaxes(ax[5, 3])
for cluster in range(0, 23):
    testdf = get_aa_frequencies(clustered_df.loc[clustered_df['Cluster'] ==cluster, 'Mut_Pos'].values, aa)
#enrichment_df = pd.DataFrame(np.where(full != 0, testdf / full, 0), columns=full.columns)
#enrichment_df.index = enrichment_df.index + 1
    if cluster < 24:
        row = cluster // 4
        col = cluster % 4

    logo = lm.Logo(testdf, color_scheme='chemistry', ax=ax[row, col])
    logo.style_spines(visible=False)
    logo.style_spines(['left', 'bottom'], visible=True)
    logo.ax.set_ylabel('Frequency', fontsize=20)
    ax[row, col].set_title(f'Cluster {cluster+1}', fontsize=25)
    ax[row, col].tick_params(axis='x', labelsize=15)
    ax[row, col].tick_params(axis='y', labelsize=15)
    #plt.savefig(f'cluster{cluster+1}_logo.pdf')
logo= lm.Logo(full, color_scheme='chemistry', ax=ax[5, 3])
logo.style_spines(visible=False)
logo.style_spines(['left', 'bottom'], visible=True)
logo.ax.set_ylabel('Frequency', fontsize=20)
ax[5, 3].set_title("Full Library", fontsize=25)
ax[5, 3].tick_params(axis='x', labelsize=15)
ax[5, 3].tick_params(axis='y', labelsize=15)

plt.savefig("clusterlogos_revised.pdf")
plt.show()



### For generating probability sequence logos for top clusters for each ligand using sequences which are >= 1.5 fold induction

ligands = ["4Hy", "End", "Ndes", "Tam", "Phlo", "Nar", "Quin", "Nal"]
# lig_clusters = {"4Hy": [17, 21, 6, 4, 18],
#                 "End": [17, 11, 6, 4, 21],
#                 "Ndes": [1, 2, 18, 17, 12],
#                 "Tam": [1, 15, 2, 12, 17],
#                 "Phlo": [21, 14, 12, 17, 18],
#                 "Nar": [14, 21, 13, 18, 12],
#                 "Quin": [15, 7, 8, 23, 19],
#                 "Nal": [5, 2, 8, 6]}

lig_clusters = {"4Hy": [17, 21, 6],
                "End": [17, 11, 6],
                "Ndes": [1, 2, 18],
                "Tam": [1, 15, 2],
                "Phlo": [21, 14, 12],
                "Nar": [14, 21, 13],
                "Quin": [15, 7, 8],
                "Nal": [5, 2, 8]}

clustered_df['Cluster'] = clustered_df['Cluster'].astype(float)
fig, ax = plt.subplots(1, 8, figsize=(60, 4))
plt.subplots_adjust(wspace=0.4)
for i, ligand in enumerate(ligands):
   clusters = lig_clusters[ligand]
   clusters = [x -1 for x in clusters]
   
   clustered_df[f'FLC_{ligand}'] = clustered_df[f'FLC_{ligand}'].astype(float)
   clustered_subset = clustered_df.loc[(clustered_df[f'FLC_{ligand}'] >= 1.5) & (clustered_df['Cluster'].isin(clusters)), 'Mut_Pos'].values
   print(clustered_subset)
   clustered_subset_freq = get_aa_frequencies(clustered_subset, aa)
   logo = lm.Logo(clustered_subset_freq, color_scheme='chemistry', ax=ax[i])
   logo.style_spines(visible=False)
   logo.style_spines(['left', 'bottom'], visible=True)
   logo.ax.set_ylabel('Frequency', fontsize=35)
   ax[i].set_title(ligand, fontsize=40)
   ax[i].tick_params(axis='x', labelsize=30)
   ax[i].tick_params(axis='y', labelsize=30)
plt.tight_layout()
plt.savefig("topclusters_logos_by_ligand.pdf")
plt.show()

### For generating enrichment sequence logos for top clusters for each ligand using sequences which are >= 1.5 fold induction 

def concatenate_fasta(input_files, output_file):
    with open(output_file, 'w') as output_fasta:
        for input_file in input_files:
            with open(input_file, 'r') as input_fasta:
                for line in input_fasta:
                    output_fasta.write(line)

def get_freqs(fasta):
    ''' fasta contains sequences of the variant positions where the abundance is determined by the functional score. Returns a dictionary of the frequency of each amino acid for each position'''
    seq_ls = []
    for record in SeqIO.parse(fasta, 'fasta'):
        seq_ls.append(str(record.seq))
    totallen = len(seq_ls)
    composition = list(zip(*seq_ls))
    composition_counts = {i+1: Counter(ls) for i, ls in enumerate(composition)}
    for key, inner_dict in composition_counts.items():
        for inner_key, value in inner_dict.items():
            composition_counts[key][inner_key] = value/totallen
    return composition_counts

#obtain lig_ref_weights from initial_weights.py
ref_counts = lig_ref_weights
aafreq_dict = {}
for ligand in ligands:
    files = glob.glob(f'{ligand}*averagesort.fasta')
    print(files)
    concatenate_fasta(files, f'{ligand}_combined_averagesort.fasta')
    ligandaa_dict = get_freqs(f'{ligand}_combined_averagesort.fasta')
    for positions, aa_dict in ligandaa_dict.items():
        for aa, frequency in aa_dict.items():
            ligandaa_dict[positions][aa] = frequency/ ref_counts[ligand][positions][aa]
    aafreq_dict[ligand] = ligandaa_dict

aafreq_df_dict = {}
aa_ls = [x for x in aa]
for ligand, aafreq in aafreq_dict.items():
    df = pd.DataFrame(None, index=np.arange(1, 12), columns=aa_ls)
    for position, aa_dict in ref_counts[ligand].items():
        for aminoacid, frequency in aa_dict.items():
            df.at[position, aminoacid] = np.log2(aafreq[position][aminoacid])

    #df.replace(0, np.nan, inplace=True)
    # Find the minimum value for each column (ignoring NaN)
    df_copy = df.copy()
    df_copy.replace(-np.inf, np.nan, inplace=True)
    minimumvalue = df_copy.min().min()
    # Replace amino acids that drop out with the minimum of the ligand set
    df.replace(-np.inf, minimumvalue, inplace=True)
    # Replace amino acids at each position that are not represented in the library with 0
    df.fillna(0, inplace=True)
    aafreq_df_dict[ligand] = df

fig, ax = plt.subplots(2, 4, figsize=(30, 8))
plt.subplots_adjust(wspace=0.4)
for i, (ligand, df) in enumerate(aafreq_df_dict.items()):
   row = i//4
   col = i % 4
   logo = lm.Logo(df, color_scheme='dmslogo_funcgroup', flip_below=False, ax=ax[row, col], vpad=0.15)
   logo.style_spines(visible=False)
   logo.style_spines(['left', 'bottom'], visible=True)
   logo.ax.set_ylabel('log2(Enrichment)', fontsize=35)
   logo.ax.set_xticks(np.arange(12), )
   logo.ax.set_xticklabels([None, 68, 70, 74, 78, 89, 92, 93, 96, 110, 113, 114], rotation=90)
   logo.ax.tick_params(axis='x', which='both', length=10)
   ax[row, col].set_title(ligand, fontsize=40)
   ax[row, col].tick_params(axis='x', labelsize=30)
   ax[row, col].tick_params(axis='y', labelsize=30)
plt.tight_layout()
plt.savefig("topclusters_enrichment_logos_by_ligand.pdf")
plt.show()





