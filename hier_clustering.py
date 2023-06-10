#python script for comparison of ligands via hierarchical clustering

import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.cluster as sc
import scipy.spatial as ss
import sklearn.metrics as sm
import sklearn.impute as si
import itertools
import sys
sys.setrecursionlimit(100000)


def read_input(file,merged):
	x = pd.read_csv(file,header=0,index_col=0)
	if merged == False:
		x['average']=x.mean(axis=1)
		return x['average']
	else:
		return x

def merge_inputs(in_list,names_list):
	concat = pd.concat(in_list,axis=1)
	concat.columns = names_list
	concat = concat.dropna(axis=0,how='any')
	return concat

def impute_data(merged_df,missing_val):
	if missing_val == '0':
		missing_val = np.int64(0)
	elif missing_val == 'nan':
		missing_val = np.nan
	else:
		pass
#	merged_df = merged_df.replace(to_replace=[0],value=np.nan)
	knn_impute = si.KNNImputer(n_neighbors=5,weights='distance',
		missing_values=missing_val)
	imputed_df = knn_impute.fit_transform(merged_df)
	print(imputed_df)
	imputed_df = pd.DataFrame(imputed_df,columns=merged_df.columns,index=merged_df.index)
	return imputed_df

def agg_cluster(merged_df,testing,outname,hmeth,hdist,max_cluster,hcmap,extension,no_log,
	hide_dendro):
	scores = {}
	metr = ['euclidean','cityblock','chebyshev','correlation','minkowski','sqeuclidean','cosine','jensenshannon','braycurtis']
	meth = ['single','average','complete','median','ward','weighted']
#	meth=['single','average']
	types = list(itertools.product(meth,metr))
	limit = 25
	if no_log == True:
		pass
	else:
		merged_df = np.log2(merged_df,out=np.zeros_like(merged_df),where=(merged_df!=0))
	if testing == True:
		temp_dict  = {}
		for x in types:
			type = '{0}_{1}'.format(x[0],x[1])
			print(type)
			try:
				h_linkage = sc.hierarchy.linkage(merged_df,method=x[0],metric=x[1])
			except:
				temp_dict[type]=0
			else:
				pair_dist = ss.distance.pdist(merged_df,metric=x[1])
				coph_corr,coph_mat = sc.hierarchy.cophenet(h_linkage,pair_dist)
				temp_dict[type]=coph_corr
		print(temp_dict)
		testing_df = pd.DataFrame.from_dict(temp_dict,orient='index')
		testing_df.columns=['Cophenetic Coefficient']
		testing_df.to_csv('{0}_cophenetic_coeff.csv'.format(outname))
		bar = sns.barplot(data=testing_df,x=testing_df.index,y='Cophenetic Coefficient',
			color='#C3B1E1')
		locs, labels = plt.xticks()
		bar.set_xticklabels(labels, rotation=45,ha='right',size=6)
		plt.tight_layout()
		plt.savefig('{0}_cophenetic_coeff.pdf'.format(outname))
		plt.clf()
		sorted_test = testing_df.sort_values(['Cophenetic Coefficient'],ascending=False)
		sorted_test = sorted_test[sorted_test.iloc[:,0]>0.70]
		if sorted_test.shape[0]>0:
			pass
		else:
			print("No acceptable method/metric via cophenetic score")
			quit()
		score_dict = {}
		score_dict['ncluster']=list(range(2,limit+1))
		for x in sorted_test.index.values:
			xtypes = x.split('_')
			print(xtypes)
			temp_score = []
			h_linkage = sc.hierarchy.linkage(merged_df,method=xtypes[0],metric=xtypes[1])
			for y in range(2,limit+1):
				clusters = sc.hierarchy.fcluster(h_linkage,criterion='maxclust',t=y)
				score = sm.silhouette_score(merged_df,clusters,metric=xtypes[1])
				temp_score.append(score)
			score_dict['{0}_{1}'.format(xtypes[0],xtypes[1])] = temp_score
		score_df = pd.DataFrame.from_dict(score_dict)
		score_df = score_df.set_index(['ncluster'])
		line = sns.lineplot(data=score_df)
		plt.tight_layout()
		plt.savefig('{0}_nclusters.pdf'.format(outname))
		plt.clf()
		score_df.to_csv('{0}_silhouette.csv'.format(outname))
	else:
		print('{0} {1}'.format(hmeth,hdist))
		print(merged_df)
		h_linkage = sc.hierarchy.linkage(merged_df,method=hmeth,metric=hdist)
		v_linkage = sc.hierarchy.linkage(merged_df.transpose(),
			method=hmeth,metric=hdist)
		xlabels = [x.split('_')[-1] for x in merged_df.columns]
		merged_df.columns = xlabels
		clustering = sc.hierarchy.fcluster(h_linkage,criterion='maxclust',t=max_cluster)
		#setting colors for clusters in the clustermap
		lut = dict(zip(set(np.sort(clustering)),sns.color_palette('tab10',max_cluster)))
		row_colors = pd.DataFrame(clustering)[0].map(lut)
		colormap = sns.color_palette(hcmap, as_cmap=True)
		if no_log == True:
			test = sns.clustermap(merged_df,row_linkage=h_linkage,col_linkage=v_linkage,
				cmap=colormap,row_colors = [row_colors],dendrogram_ratio=(0.2,0.2),
				vmin=0,figsize=(12,8),robust=False,
				cbar_kws=dict(orientation='horizontal')
				)

		else:
			test = sns.clustermap(merged_df,row_linkage=h_linkage,col_linkage=v_linkage,
				cmap=colormap,row_colors = [row_colors],dendrogram_ratio=(0.2,0.2),
				center=0,figsize=(12,8),robust=False,
				cbar_kws=dict(orientation='horizontal')
				)
		axcbar = test.ax_cbar
		x0, _y0, _w, _h = test.cbar_pos
		axcbar.locator_params(nbins=5)
		ax = test.ax_heatmap
		ax.set(yticklabels=[])
#		ax.set(xticklabels=xlabels)
		ax.tick_params(right=False)
		ax.tick_params(axis='both', which='major', labelsize=18)
		axcbar.tick_params(labelsize=18)
		ax.set_ylabel('Variants',rotation=-90,fontsize=24,labelpad=14.0)
		ax.set_xlabel('Ligands',fontsize=24)
		if hide_dendro == True:
			test.ax_row_dendrogram.set_visible(False)
			test.ax_col_dendrogram.set_visible(False)
			axcbar.set_position([0.5, 0.9, 
				test.ax_row_dendrogram.get_position().width, 0.05])
		else:
			axcbar.set_position([x0, 0.9, 
				test.ax_row_dendrogram.get_position().width, 0.05])
		axcbar.set_title('Normalized Fold Enrichment',fontsize=18)
		for spine in ax.spines.values():
			spine.set_visible(True)
			spine.set_color('black')
		plt.savefig('{0}_{1}_{2}_{3}_cluster.{4}'.format(outname,
			hmeth,hdist,max_cluster,extension),bbox_inches='tight',transparent=True)
		merged_df['cluster_id']=clustering
		merged_df = merged_df.sort_values(['cluster_id'])
		print(merged_df)
		merged_df.to_csv('{0}_{1}_{2}_{3}_clusterlab.csv'.format(outname,
			hmeth,hdist,max_cluster))
		return merged_df
	
def create_boxplots(merged_df,max_cluster,outname):
	for x in range(1,max_cluster+1):
		temp_df = merged_df[merged_df['cluster_id']==x]
		fig,ax = plt.subplots(1,1)
		fig.set_size_inches(8,6)
		sns.violinplot(data=temp_df.iloc[:,:-1],cut=0,scale='count',inner='box',
			ax=ax,color='#a6bddb',fliersize=1,saturation=1)
		ax.set_ylabel('Normalized Fold Enrichment')
		ax.set_xlabel('Ligands')
#		ax.set_xticklabels(temp_df.columns.values[:-1],fontsize=8,rotation=45,ha='right')
		plt.tight_layout()
		plt.savefig('{0}_{1}_vplot.pdf'.format(outname,x))
		plt.clf()		

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='PCA analysis of variants across conditions')
	parser.add_argument('--input','-i',required=True,nargs='*',
		help='Input files - 1 for each condition')
	parser.add_argument('--suffix','-s',required=False,default='',
		help='Suffix to remove for distinguishing inputs')
	parser.add_argument('--merged','-m',required=False,action='store_true',default=False,
		help='Call this flag if the input is already merged. Ignores suffix')
	parser.add_argument('--impute','-p',required=False,default=None,
		choices=[None,'nan','0'],
		help='Call this flag if there is missing data that must be imputed. Uses a KNN imputation from SciKit Learn. Missing data is represented by either 0 or np.nan')
	parser.add_argument('--outname','-o',required=False,default=None,
		help='Outfile prefix')
	parser.add_argument('--testing','-t',required=False,default=False,action='store_true',
		help='Call this flag to determine optimal distance metric and method parameters')
	parser.add_argument('--method','-d',required=False,default='ward',
		choices=['single','average','complete','median','ward','centroid','weighted'],
		help='Choose clustering method. Default is ward')
	parser.add_argument('--metric','-c',required=False,default='euclidean',
		choices=['euclidean','cityblock','chebyshev','correlation','minkowski','sqeuclidean','cosine','jensenshannon','braycurtis'],
		help='Choose distance metric. Default is euclidean')
	parser.add_argument('--max_cluster','-mc',required=False,default=8,
		help='Specify cluster #')
	parser.add_argument('--heatmap_cmap','-hc',required=False,default='crest',
		help='Heatmap colormap')
	parser.add_argument('--remove_bad','-r',required=False,default=False,
		action='store_true',
		help='Call this flag to remove all variants that do not perform at a threshold above wildtype. Threshold dictated by --remove_threshold')
	parser.add_argument('--remove_threshold','-rt',required=False,default=1.0,
		help='Define threshold for --remove_bad. Default is 1.0. Value is for variant/wt')
	parser.add_argument('--extension','-e',required=False,default='pdf',
		help='Specify graph format. Default is pdf. May want to choose png or tiff')
	parser.add_argument('--no_log','-nl',required=False,default=False,action='store_true',
		help='Call this flag to prevent log2 transformation')
	parser.add_argument('--hide_dendro','-hd',required=False,default=False,
		action='store_true',
		help='Call this flag to hide dendrograms')
	parsed = parser.parse_args()
	if parsed.outname is None:
		outname = parsed.input[0].split('/')[-1].split('.')[0]
	else:
		outname = parsed.outname
	if parsed.merged == False:
		name_list = [x.removesuffix('_{0}'.format(parsed.suffix)) for x in parsed.input]
		in_list = [read_input(x,False) for x in parsed.input]
		merged_df = merge_inputs(in_list,name_list)
	else:
		merged_df = read_input(parsed.input[0],True)
	if parsed.impute is not None:
		merged_df = impute_data(merged_df,parsed.impute)
		merged_df.to_csv('{0}_impute.csv'.format(outname))
	else:
		pass
	if parsed.remove_bad == True:
		merged_df = merged_df[(merged_df>np.float64(parsed.remove_threshold)).any(axis=1)]
		merged_df.to_csv('{0}_threshold.csv'.format(outname))
	else:
		pass
	merged_df = agg_cluster(merged_df,parsed.testing,outname,parsed.method,parsed.metric,
		int(parsed.max_cluster),parsed.heatmap_cmap,parsed.extension,parsed.no_log,
		parsed.hide_dendro)
	if parsed.testing == False:
		create_boxplots(merged_df,int(parsed.max_cluster),outname)
	else:
		pass
	

