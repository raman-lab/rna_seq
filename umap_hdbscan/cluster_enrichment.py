import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection

#Script to identify enriched clusters for each ligand

cluster_df = pd.read_csv('avg_foldenrichment_full_clustered.csv')

clusters = list(map(lambda x: x+1, sorted(cluster_df['Cluster'].unique())[1:]))
ligands = ['FLC_4Hy', 'FLC_End', 'FLC_Ndes', 'FLC_Tam', 'FLC_Quin', 'FLC_Nal', 'FLC_Phlo', 'FLC_Nar']

xvalues, yvalues = np.meshgrid(np.arange(len(ligands)), np.arange(len(clusters)))

avgvalues = np.array([]) # average values of each cluster-ligand combination (>=1.5 F-Score)
hitsperc = np.array([]) # %hit (F-score >= 1.5)
numhits = np.array([]) # Number of hits (F-score >= 1.5)

for cluster in clusters:
    subset = cluster_df.loc[cluster_df['Cluster'] == cluster-1]
    totalvariants = subset.shape[0]
    for ligand in ligands:
        subset2 = subset.loc[subset[ligand] >= 1.5]
        numhits = np.append(numhits, subset2.shape[0])
        hitsperc = np.append(hitsperc, subset2.shape[0]/totalvariants)
        avgvalues = np.append(avgvalues, subset2.loc[:, ligand].mean())

avgscores_dict = {ligand: None for ligand in ligands}
hitperc_dict = {ligand: None for ligand in ligands}
numhits_dict = {ligand: None for ligand in ligands}
hitsperc2 = [0] * len(hitsperc)
for i, ligand in enumerate(ligands):
    chunk = hitsperc[i::8]
    chunk2 = avgvalues[i::8]
    chunk3 = numhits[i::8]
    hitperc_dict[ligand] = chunk
    avgscores_dict[ligand] = chunk2
    numhits_dict[ligand] = chunk3
    #chunk_sum = sum(chunk)
    chunkavg = [ratio/max(chunk)/2 for ratio in chunk]
    largest = max(chunkavg)
    for x, item in enumerate(chunkavg):
        hitsperc2[i + x * 8] = item

hiteperc_df = pd.DataFrame(hitperc_dict)
hiteperc_df.index = clusters
hiteperc_df.to_csv('cluster_hitpercentage_fscore1.5.csv')
avgscores_df = pd.DataFrame(avgscores_dict)
avgscores_df.index = clusters
avgscores_df.to_csv('cluster_avgscore_fscore1.5.csv')
numhits_df = pd.DataFrame(numhits_dict)
numhits_df.index = clusters
numhits_df.to_csv('cluster_numhits_fscore1.5.csv')
#avgvalues_capped = [4 if x >= 4 else x for x in avgvalues.flatten()]

fig, ax = plt.subplots()
ax.set_aspect('equal')

circles = [plt.Circle((j,i), radius=r, lw=0) for r, j, i in zip(hitsperc2, xvalues.flat, yvalues.flat)]
col = PatchCollection(circles, array=avgvalues, cmap="plasma", match_original=True)
col.set_clim(vmin=0, vmax=max(avgvalues))
ax.add_collection(col)

ax.set(xticks=np.arange(len(ligands)), yticks=np.arange(len(clusters)),
       xticklabels=ligands, yticklabels=clusters)
ax.set_xticklabels(ligands, rotation=90)
ax.set_xticks(np.arange(len(ligands)+1)-0.5, minor=True)
ax.set_yticks(np.arange(len(clusters)+1)-0.5, minor=True)
ax.grid(which='major')

colorbar = fig.colorbar(col, ax=ax)
# Adjust the position of the color bar and label


# Add legend for circles
ax.add_patch(plt.Circle((9, 5), 0.05, color='black', clip_on=False, lw=0))
ax.add_patch(plt.Circle((9, 3.5), 0.125, color='black', clip_on=False, lw=0))
ax.add_patch(plt.Circle((9, 2), 0.25, color='black', clip_on=False, lw=0))
ax.add_patch(plt.Circle((9, 0.5), 0.5, color='black', clip_on=False, lw=0))

ax.text(10.75, 5, 0.1, fontsize=10, ha='center', va='center')
ax.text(10.75, 3.5, 0.25, fontsize=10, ha='center', va='center')
ax.text(10.75, 2, 0.5, fontsize=10, ha='center', va='center')
ax.text(10.75, 0.5, 1, fontsize=10, ha='center', va='center')

plt.tight_layout()
colorbar.set_label('Avg. F-Score', rotation=-90, labelpad=15, fontsize=10)
colorbar.ax.set_position([0.8, 0.45, 0.05, 0.5]) 
plt.savefig('umapcluster_circleheatmap2.pdf')
plt.show()
