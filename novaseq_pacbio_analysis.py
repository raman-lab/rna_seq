'''
script to select pacbio barcodes from a novaseq dataset.
PacBio barcodes should be in the format:
mutant,barcode1,barcode2,barcode3,...
Novaseq data should be in the format:
barcode,fold_induction
'''
import argparse
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from collections import OrderedDict

def read_barcodes(barcode_file):
	x = pd.read_csv(barcode_file,header=0,index_col=None,dtype='str')
	x.rename(columns={x.columns[0]: 'mutants'},inplace=True)
	return x

def read_fis(fi_file):
	x = pd.read_csv(fi_file)
	return x

def read_control_names(cn_file):
	x = pd.read_csv(cn_file,header=None)
	x.columns=['labels','barcode']
	x.set_index(x.columns[1],inplace=True)
	return x	

def intersection(barcode_file,fi_file,use_col,verbose,use_median):
	out_dict = OrderedDict()
	metrics_df = pd.DataFrame()
	for row in range(0,barcode_file.shape[0]):
		bcs = pd.DataFrame(barcode_file.iloc[row,1:])
		bcs.dropna(how='any',axis=0,inplace=True)
		bcs.columns=['barcode']
		intersection = pd.merge(bcs,fi_file,how='inner',on=['barcode'])
		print intersection
		if verbose == True:
			print "Variant: ",barcode_file.iat[row,0]
			print "Variant Barcodes: ",bcs.shape[0]
			print "Pacbio_Novaseq_shared_barcodes: ",intersection.shape[0]
		else:
			pass
		out_dict[barcode_file.iat[row,0]]=intersection.iloc[:,use_col].tolist()
		if use_median == False:
			temp_df = pd.DataFrame(
				[[barcode_file.iat[row,0],np.nanmean(intersection.iloc[:,use_col]),np.nanstd(intersection.iloc[:,use_col])]])
		else:
			temp_df = pd.DataFrame(
				[[barcode_file.iat[row,0],np.nanmedian(intersection.iloc[:,use_col]),np.nanstd(intersection.iloc[:,use_col])]])
		metrics_df = metrics_df.append(temp_df,ignore_index=True)
	if use_median == False:
		metrics_df.columns = ['variant','average','std_dev']
	else:
		metrics_df.columns = ['variant','median','std_dev']
	out_df = pd.DataFrame.from_dict(out_dict,orient='index')
	return out_df, metrics_df

def control_intersection(indict):
	dfs = [indict[x] for x in indict.keys()]
	mdf = pd.concat(dfs,axis=1,join='inner')
	mdf.reset_index(inplace=True)
	return mdf
	
def sum_input(in_dict):
	sum_dict = OrderedDict()
	for k in in_dict:
		sum_dict[k] = in_dict[k].sum(axis=1,skipna=True,numeric_only=True)
	sum_df = pd.DataFrame.from_dict(sum_dict,orient='columns')
	return sum_df

def sum_ratios(sum_df,inner,outer,total_df):
	inner = [x.split(',') for x in inner]
	outer = [x.split(',') for x in outer]
	in_lst = [(int(x[0]),int(x[1])) for x in inner]
	out_lst = [(int(x[0]),int(x[1])) for x in outer]
	file_list = sum_df.columns.values.tolist()
	first_ratio_df = pd.DataFrame()
	final_ratio_df = pd.DataFrame()
	if inner is not None:
		for x in in_lst:
			first_ratio = sum_df[file_list[x[0]]]/sum_df[file_list[x[1]]]
			total_ratio = total_df[file_list[x[1]]]/total_df[file_list[x[0]]]
			first_ratio_df[file_list[x[0]]]=first_ratio*total_ratio[0]
	if outer is not None:
		for x in out_lst:
			second_ratio=first_ratio_df[file_list[x[0]]]/first_ratio_df[file_list[x[1]]]
			final_ratio_df[file_list[x[0]]]=second_ratio
	return sum_df,first_ratio_df,final_ratio_df
		
def normalize(sample_df, control_ratio_df):
	final_avg = np.mean(control_ratio_df.iat[0,0])
	normalized_sample_df = pd.DataFrame()
	for x in sample_df.columns.values:
		normalized_sample_df[x]=sample_df[x]/final_avg
	normalized_sample_df.index.name='variant'
	return normalized_sample_df

def plot_averages(out_df,name):
	fig,ax = plt.subplots(1,1)
	fig.set_size_inches(8,6)
	xvals = np.arange(0,out_df.shape[0])
	ax.bar(xvals+0.5,out_df.iloc[:,1],yerr=out_df.iloc[:,2],
		width=1,edgecolor='black',color='#77DD77',ecolor='black',
		linewidth=1,capsize=5,align='center',
		error_kw={'elinewidth': 1,'ecapthick': 2}
		)
	ax.set_xlim(xmin=0)
	ax.set_xticks(xvals+0.5)
	ax.set_xticklabels(out_df.iloc[:,0],rotation=45,ha='right')
	ax.set_xlabel('Variants')
	ax.set_ylabel('NovaSeq_FI')
	ax.spines['right'].set_visible=False
	ax.spines['top'].set_visible=False
	plt.tight_layout()
	plt.savefig('{0}_averages.pdf'.format(name))

def plot_box(out_df,name):
	fig,ax = plt.subplots(1,1)
	fig.set_size_inches(8,6)
	to_plot = [np.array(out_df.iloc[x,:].dropna(how='any',axis=0)) for x in range(0,out_df.shape[0])]
	ax.boxplot(to_plot,sym='.')
	ax.set_xticklabels(out_df.index.values,rotation=45,ha='right')
	plt.tight_layout()
	plt.savefig('{0}_boxplot.pdf'.format(name))

if __name__ == '__main__':
	parser = argparse.ArgumentParser(prefix_chars='-+',description='''
		script to select pacbio barcodes from a novaseq dataset.
		PacBio barcodes should be in the format:
		mutant,barcode1,barcode2,barcode3,...
		Novaseq data should be in the format:
		barcode,fold_induction
		''')
	parser.add_argument('--barcode_file','-b',required=True,
		help='Barcode input file from PacBio analysis')
	parser.add_argument('--fi_file','-f',required=True,
		help='Fold induction file from NovaSeq analysis')
	parser.add_argument('--out','-o',required=False,default=None,
		help='Specify outfile name')
	parser.add_argument('--use_col','-u',required=False,default=1,
		help='0-indexed column in the fi_file to use during pacbio-novaseq matching')
	parser.add_argument('--sum','-s',required=False,default=False,action='store_true',
		help='''
		Call this flag to sum up the counts of all barcodes corresponding to a variant
		This requires the use of the ++totals flag and will fail if not provided.
		''')
	parser.add_argument('++sum_totals','+t',required=False,default=None,
		help='Specify the location of a *_totals.csv containing total reads')
	parser.add_argument('++sum_inner','+i',required=False,default=None,nargs='*',
		help='''
		Comma separated numbers of indices in --input that go together. These will be ratio'd first (first number is numerator, second is denominator) before the outer pairs. If there are multiple pairs, separate them by a space. Ex: '1,2 3,4'. Use this flag if a ratio of the totals is required.
		'''
		)
	parser.add_argument('++sum_outer','+o',required=False,default=None,nargs='*',
		help='''
		Comma separated numbers of indices in --input that go together. These will be ratio'd (first number is numerator, second is denominator) after the inner pairs. Please use only numbers found in the NUMERATOR of the inner set. If there are multiple pairs, separate them by a space. Use this flag if a ratio of the inners is required
		'''
		)
	parser.add_argument('++sum_controls','+c',required=False,default=None,
		help='''
		Call this flag to normalize the summed data to a set of known controls. CSV file with name,barcode (no header).
		'''
		)
	parser.add_argument('--median','-m',required=False,default=False,action='store_true',
		help='Call this flag to force median analysis instead of mean')
	parsed = parser.parse_args()
	barcode_df = read_barcodes(parsed.barcode_file)
#	print barcode_df
	fi_df = read_fis(parsed.fi_file)
#	print fi_df
	if parsed.out == None:
		basename = parsed.barcode_file.split('/')[-1].split('.')[0]
	else:
		basename = parsed.out
	if parsed.sum == False:
		merge_df,metrics_df = intersection(barcode_df,fi_df,int(parsed.use_col),
			True,parsed.median)
#		print merge_df
		plot_box(merge_df,basename)
		merge_df.to_csv('{0}_fi_data.csv'.format(basename),index=False)
		metrics_df.to_csv('{0}_metrics.csv'.format(basename),index=False)
	else:
		total_df = read_fis(parsed.sum_totals)
		print total_df
		fi_dict = OrderedDict()		
#		print fi_df
		for x in range(1,fi_df.shape[1]):
			merge_df, metrics_df = intersection(barcode_df,fi_df,x,False,parsed.median)
			fi_dict[fi_df.columns[x]]=merge_df
		print fi_dict.keys()
		sum_input_df = sum_input(fi_dict)
		sum_df,first_ratio_df,final_ratio_df = sum_ratios(
			sum_input_df,parsed.sum_inner,parsed.sum_outer,total_df)		
		
		if parsed.sum_controls is not None:
			control_names = read_control_names(parsed.sum_controls)
			fi_df.set_index(fi_df.columns[0],inplace=True)
			temp_dict = {
				'sums': fi_df,
				'control_names': control_names
				}
			control_data = control_intersection(temp_dict)
			print control_data
			summed_controls = pd.DataFrame(control_data.sum(
				axis=0,skipna=True,numeric_only=True))
#			print summed_controls
			control_sum,control_first_ratio,control_final_ratio = sum_ratios(
				 summed_controls.transpose(),parsed.sum_inner,parsed.sum_outer,total_df)
			print control_final_ratio
			normalized = normalize(final_ratio_df,control_final_ratio)
			print normalized
			plot_box(normalized,basename)
			normalized.reset_index(inplace=True)
			normalized.to_csv('{0}_normalized.csv'.format(basename),index=False)
			
		else:
			plot_box(sum_df,basename)
		sum_df.to_csv('{0}_sum_data.csv'.format(basename))
		if first_ratio_df.empty == False:
			first_ratio_df.to_csv('{0}_first_ratio.csv'.format(basename),index=False)
		if final_ratio_df.empty == False:
			final_ratio_df.to_csv('{0}_final_ratio.csv'.format(basename),index=False)
		

	