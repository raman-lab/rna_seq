#Code for UMAP-HDBSCAN analysis
import umap.umap_ as umap
import hdbscan
import random
from tqdm import trange
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import statistics as stat
import os
import constants
from collections import Counter
import SeqIO

# File containing F-scores for each variant-ligand pair. File here has Cluster column added already
vl_df_full = pd.read_csv("avg_foldenrichment_full_clustered.csv")

## Create a physiochemical encoding for all variants using the 'Mut_Pos' column and the dimensinality-reduced AAIndex (Sam Gelman, Sarah A Fahlberg, Pete Heinzelman, Philip A Romero+, Anthony Gitter+. Proceedings of the National Academy of Sciences, 118:48, 2021.)
def enc_aa_index(int_seqs):
    """ encodes data in aa index properties format """
    aa_features = np.load("pca-19.npy")
    # add all zero features for stop codon
    aa_features = np.insert(aa_features, 0, np.zeros(aa_features.shape[1]), axis=0)
    aa_features_enc = aa_features[int_seqs]
    return aa_features_enc

def enc_int_seqs_from_char_seqs(char_seqs):
    seq_ints = []
    for char_seq in char_seqs:
        int_seq = [constants.C2I_MAPPING[c] for c in char_seq]
        seq_ints.append(int_seq)
    seq_ints = np.array(seq_ints)
    return seq_ints

def pc_encode(char_seqs, indexes=None):
    int_seqs = enc_int_seqs_from_char_seqs(char_seqs)
    encoded = enc_aa_index(int_seqs)
    encoded_zip = list(zip(*encoded))
    
    dfs_to_merge = []
    for columnset in encoded_zip:
        dfs_to_merge.append(pd.DataFrame(columnset))
    merged_df = pd.concat(dfs_to_merge, axis=1)
    
    if indexes is not None:
        merged_df = merged_df.set_index(pd.Index(indexes))
    print(merged_df.shape)
    return merged_df

pc_encoded = pc_encode(vl_df_full['Mut_Pos'], indexes=vl_df_full.index.tolist())
pc_encoded.to_csv("pcencoded.csv")


## Random UMAP parameter search. Taken from David Borrelli (https://towardsdatascience.com/clustering-sentence-embeddings-to-identify-intents-in-short-text-48d22d3bf02e)
def generate_clusters(message_embeddings,
                      n_neighbors,
                      n_components,
                      min_dist,
                      min_cluster_size,
                      random_state = None):
    """
    Generate HDBSCAN cluster object after reducing embedding dimensionality with UMAP
    """
    
    umap_embeddings = (umap.UMAP(n_neighbors=n_neighbors, 
                                n_components=n_components,
                                min_dist=min_dist,
                                metric='manhattan', 
                                random_state=42)
                            .fit_transform(message_embeddings))

    clusters = hdbscan.HDBSCAN(min_cluster_size = min_cluster_size,
                               metric='euclidean', 
                               cluster_selection_method='eom').fit(umap_embeddings)

    return clusters

def score_clusters(clusters, prob_threshold = 0.05):
    """
    Returns the label count and cost of a given cluster supplied from running hdbscan
    """
    
    cluster_labels = clusters.labels_
    label_count = len(np.unique(cluster_labels))
    total_num = len(clusters.labels_)
    cost = (np.count_nonzero(clusters.probabilities_ < prob_threshold)/total_num)
    
    return label_count, cost

def random_search(embeddings, space, num_evals):
    """
    Randomly search hyperparameter space and limited number of times 
    and return a summary of the results
    """
    
    results = []
    
    for i in trange(num_evals):
        n_neighbors = random.choice(space['n_neighbors'])
        n_components = random.choice(space['n_components'])
        min_dist = random.choice(space['min_dist'])
        min_cluster_size = random.choice(space['min_cluster_size'])
        
        clusters = generate_clusters(embeddings, 
                                     n_neighbors = n_neighbors, 
                                     n_components = n_components, 
                                     min_dist = min_dist,
                                     min_cluster_size = min_cluster_size, 
                                     random_state = 42)
    
        label_count, cost = score_clusters(clusters, prob_threshold = 0.05)
                
        results.append([i, n_neighbors, n_components, min_dist, min_cluster_size, 
                        label_count, cost])
    
    result_df = pd.DataFrame(results, columns=['run_id', 'n_neighbors', 'n_components','min_dist', 
                                               'min_cluster_size', 'label_count', 'cost'])
    
    return result_df.sort_values(by='cost')

space = {
        'n_neighbors': [15, 25, 50, 75, 100, 300, 500, 1000, 2000],
        'n_components': [2],
        'min_dist': [0, 0.1],
        'min_cluster_size': range(5, 110, 5),
    }
results = random_search(pc_encoded, space, 100)

#   500 n_neighbor, 0.1 min_dist, manhattan metrix, and min_cluster_size=80 appears to give the minimize the cost (0.03144) at n_component=2


# UMAP embedding
embedding = umap.UMAP(
    n_neighbors=500,
    min_dist=0.1,
    n_components=2,
    metric='manhattan',
    random_state=42,
).fit_transform(pc_encoded)

# HDBSCAN
labels = hdbscan.HDBSCAN(
    min_cluster_size=80
).fit_predict(embedding)

clustered = (labels >= 0)
unique_clusters = np.unique(labels)

# Visualize
plt.figure(figsize=(9, 6))
plt.scatter(embedding[~clustered, 0],
            embedding[~clustered, 1],
            color=(0.5, 0.5, 0.5),
            s=0.1,
            alpha=0.5)
scatter = plt.scatter(embedding[clustered, 0],
            embedding[clustered, 1],
            c=labels[clustered],
            s=0.1,
            cmap='Spectral');
legend_labels = [f'{cluster}' for cluster in unique_clusters][1:]
num_clusters = len(unique_clusters)-1
legend_handles = [plt.Line2D([], [], marker='o', markersize=8, color='white', markerfacecolor=plt.cm.get_cmap('Spectral')(cluster/(num_clusters-1))) for cluster in unique_clusters[1:]]
plt.legend(legend_handles, legend_labels, loc='center left', bbox_to_anchor=(1.0, 0.5), ncol=2, fontsize=8, title='Cluster')
plt.title('500 neighbors, 0.1 min dist, manhattan, HDBSCAN')
plt.tight_layout()
plt.show()

# Add cluster assignments
vl_clustered = vl_df_full.copy()
vl_clustered['Cluster'] = labels.tolist()
vl_clustered.to_csv('avg_foldenrichment_full_clustered.csv')

# Get F-score overlay on UMAP for each ligand
ligand = 'Ndes'
ligands = ['4Hy', 'Ndes', 'End', 'Tam', 'Phlo', 'Nar', 'Quin', 'Nal']
embed_vl_map = {embed_index:vl_index for embed_index, vl_index in enumerate(vl_df_full.index.tolist())}
vl_embed_map = {vl_index:embed_index for embed_index, vl_index in enumerate(vl_df_full.index.tolist())}

for ligand in ligands:
    sorted_df = vl_df_full.sort_values(f'FLC_{ligand}')
    sorted_indices = np.array([vl_embed_map[index] for index in sorted_df.index])
    sorted_embedding = embedding[sorted_indices]
    plt.scatter(
        sorted_embedding[:, 0], 
        sorted_embedding[:, 1], 
        c=sorted_df[f'FLC_{ligand}'],
        cmap = 'Purples',
        vmin=0, 
        vmax=4,
        s=0.25)
    plt.colorbar()
    plt.title(f'500 neighbors, 0.1 min dist, manhattan, {ligand}')
    plt.tight_layout()
    plt.show()

# Average functional score of cluster
# Chosen clusters had the highest average functional score with minimum 15 hits (score => 1.5)
vl_clustered = pd.read_csv('avg_foldenrichment_full_clustered.csv', index_col=0)

#To get cluster counts
cluster_counts = dict(Counter(vl_clustered['Cluster']))
cluster_counts = {key: cluster_counts[key] for key in sorted(cluster_counts)}
plt.bar(cluster_counts.keys(), cluster_counts.values())
for i, v in enumerate(cluster_counts.values()):
    plt.text(i-1, v + 100, str(v), ha='center', va='bottom', rotation=90, fontsize=9)
plt.xticks([x for x in range(-1,23)], rotation=90)
plt.ylim((0, 4500))
plt.ylabel('Number of Variants')
plt.xlabel('Cluster')
plt.tight_layout()
plt.savefig('clustersizes.pdf')
plt.show()

# To get average score of clusters
cluster_average = vl_clustered.iloc[:, :9]
cluster_average['Cluster'] = vl_clustered['Cluster']
cluster_average = cluster_average.groupby('Cluster').mean()
ligands = ['4Hy', 'End', 'Ndes', 'Tam', 'Phlo', 'Nar', 'Quin', 'Nal', 'EllA']
fig, axes = plt.subplots(nrows=9, ncols=1, figsize=(10, 30), sharex=True)
for i, lig in enumerate(ligands):
    axes[i].bar(cluster_average.index, cluster_average[f'FLC_{lig}'])
    axes[i].set_title(lig, fontsize=25, loc='right')
    axes[i].set_xlabel('Cluster')
    axes[i].set_ylabel('Average F-Score')
    axes[i].set_xticks([x for x in range(-1, 23)],labels=[x for x in range(-1, 23)], rotation=90)
    axes[i].tick_params(axis='both', which='major', labelsize=22)
plt.tight_layout()
plt.savefig('cluster_averagefscores.pdf')
plt.show()

# To get average score of clusters after 1.5 minimum fscore filter
filtered_df = vl_clustered.iloc[:, :9][vl_clustered.iloc[:,:9] >=1.5]
filtered_df['Cluster'] = vl_clustered['Cluster']
filtered_av = filtered_df.groupby('Cluster').mean()
ligands = ['4Hy', 'End', 'Ndes', 'Tam', 'Phlo', 'Nar', 'Quin', 'Nal', 'EllA']
fig, axes = plt.subplots(nrows=9, ncols=1, figsize=(10, 30), sharex=True)
for i, lig in enumerate(ligands):
    axes[i].bar(filtered_av.index, filtered_av[f'FLC_{lig}'])
    axes[i].set_title(lig, fontsize=25, loc='right')
    axes[i].set_xlabel('Cluster')
    axes[i].set_ylabel('Average F-Score')
    axes[i].set_xticks([x for x in range(-1, 23)],labels=[x for x in range(-1, 23)], rotation=90)
    axes[i].tick_params(axis='both', which='major', labelsize=22)
plt.tight_layout()
plt.savefig('cluster_averagefscores_min1.5.pdf')
plt.show()

counts_df = filtered_df.groupby('Cluster').count()
q1 = np.quantile(counts_df.values.flatten(), 0.25) # 11.75 is the first quantile. 15 is a good threshold

# To get number of variants in each cluster that pass the 1.5 minimum fscore threshold
ligands = ['4Hy', 'End', 'Ndes', 'Tam', 'Phlo', 'Nar', 'Quin', 'Nal', 'EllA']
fig, axes = plt.subplots(nrows=9, ncols=1, figsize=(10, 30), sharex=True)
for i, lig in enumerate(ligands):
    axes[i].bar(counts_df.index, counts_df[f'FLC_{lig}'])
    axes[i].set_title(lig, fontsize=25, loc='right')
    axes[i].set_xlabel('Cluster')
    axes[i].set_ylabel('Number of Variants (\u226515 F-score)')
    axes[i].set_xticks([x for x in range(-1, 23)],labels=[x for x in range(-1, 23)], rotation=90)
    axes[i].tick_params(axis='both', which='major', labelsize=22)
plt.tight_layout()
plt.savefig('clustersizes_min1.5.pdf')
plt.show()

average_df = filtered_df.groupby('Cluster').apply(lambda x: x.mean())

## Ranking based on highest average functional score with a minimum of 15 hits
for ligand in ligands:
    avg_drop_df = average_df.drop(average_df.index[counts_df.loc[:, f'FLC_{ligand}'] < 15])
    avg_drop_df = avg_drop_df.sort_values(f'FLC_{ligand}', ascending=False)
    try:
        avg_drop_df = avg_drop_df.drop(-1)
    except:
        pass
    top_five = avg_drop_df.index[:3].tolist()
    for index, cluster in enumerate(top_five):
        cluster_df = vl_clustered.loc[(vl_clustered['Cluster'] == cluster) & (vl_clustered[f'FLC_{ligand}'] >= 1.5)]
        cluster_df.to_csv(f'{ligand}_cluster{cluster+1}_rank{index+1}_1.5_hits_15min_averagesort.csv')


# If heatmap is to be made by calculating probabilities of each mutant position weighted by functional score 
def composition(iterable):
    ''' iterable contains another iterable of which to count'''
    totallen = len(iterable)
    composition = list(zip(*iterable))
    composition_counts = {i+1: Counter(ls) for i, ls in enumerate(composition)}
    for key, inner_dict in composition_counts.items():
        for inner_key, value in inner_dict.items():
            composition_counts[key][inner_key] = value/totallen
    return composition_counts
    
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
        
def get_prob(ref_dict, exp_dict):
    ''' returns a pandas Series containing the raw probibilities of each amino acid at each mutant position in the sequences of exp_dict '''
    indexls = []
    valuels = []
    #out_dict = copy.deepcopy(ref_dict)
    for outerkey, outervalue in ref_dict.items():
        for innerkey, innervalue in outervalue.items():
            try:
                valuels.append(exp_dict[outerkey][innerkey])
                indexls.append(f'{outerkey}_{innerkey}')
                #out_dict[outerkey][innerkey] = exp_dict[outerkey][innerkey]/innervalue
            except:
                valuels.append(0)
                #out_dict[outerkey][innerkey] = 0
    return pd.Series(valuels, index=indexls)

def get_enrichprob(ref_dict, exp_dict):
    ''' returns a pandas Series containing the probibilities of each amino acid at each mutant position in the sequences of exp_dict normalized to the frequency of the amino acid at that position in the ref_dict '''
    indexls = []
    valuels = []
    for outerkey, outervalue in ref_dict.items():
        for innerkey, innervalue in outervalue.items():
            try:
                valuels.append(exp_dict[outerkey][innerkey]/innervalue)
                indexls.append(f'{outerkey}_{innerkey}')
            except:
                valuels.append(0)
    return pd.Series(valuels, index=indexls)

def log2_transform_with_zeros(df):
    min_non_zero = df[df>0].min().min()
    df.replace(0, min_non_zero, inplace=True)
    df = np.log2(df)
    return df


directory = "/Volumes/Jackie/TtgR_RNASeq/umap"
all_files = {}
for file in os.listdir(directory):
    filename = os.fsdecode(file)
    if filename.endswith("_averagesort.fasta"):
        splitname = filename.split('_')
        ligand = splitname[0]
        cluster = splitname[1]
        rank = int(splitname[2][4:])
        all_files[(ligand, cluster, rank)] = get_freqs(filename)

#obtain lig_ref_weights from initial_weights.py
ref_counts = lig_ref_weights

ligands = ['4Hy', 'End', 'Ndes', 'Tam', 'Phlo', 'Nar', 'Quin', 'Nal']
lig_dfs = []
for ligand in ligands:
    cluster_FE = {}
    subset = {key: value for key, value in all_files.items() if key[0] == ligand}
    sorteditems = sorted(subset.items(), key=lambda x: x[0][2])
    sorted_subset = dict(sorteditems)
    for ligand_cluster, counter in sorted_subset.items():
        #cluster_FE[ligand_cluster[1]] = get_prob(ref_counts, counter)
        cluster_FE[ligand_cluster[1]] = get_enrichprob(ref_counts[ligand], counter)
    lig_dfs.append(pd.DataFrame(cluster_FE))
lig_dfs = list(map(lambda x: log2_transform_with_zeros(x), lig_dfs))
numcols = len(lig_dfs)


fig, axes = plt.subplots(ncols=numcols+1, figsize=(6*numcols, 25))
for i, df in enumerate(lig_dfs):
    df = df.iloc[:, :3]
    #df.to_csv(f'{ligands[i]}_enrichmentscores.csv')
    ax = axes[i]
    heatmap = sns.heatmap(df, cmap='coolwarm', linewidths=3, linecolor='white', ax=ax, cbar=False)
    ax.tick_params(axis='x', which='both', labelsize=30, rotation=90)
    ax.tick_params(axis='y', which='both', labelsize=30, rotation=0)
    ax.set_title(ligands[i], fontsize=40)
    ax.set_aspect('equal')
cbar_ax = fig.add_axes([0.90, 0.15, 0.02, 0.7])
cbar = plt.colorbar(heatmap.collections[0], cax=cbar_ax)
cbar.ax.tick_params(labelsize=40)
axes[numcols].remove()
cbar.set_label('log2(Enrichment)', fontsize=60, rotation=-90, labelpad=70)
plt.tight_layout(pad=3)
plt.show()
