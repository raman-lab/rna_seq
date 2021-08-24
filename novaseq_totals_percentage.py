#script to analyze differences in barcode counts as a function of threshold
import dask.dataframe as dd
import dask
import pandas as pd
import numpy as np
import argparse
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

def read_combined_input(input_list):
	x = dd.read_csv(input_list,header=0)
	x = x.set_index(x.columns[0])
	return x

def delayed_shell(input_list,cutoff):
	results = []
	for input in input_list:
		base = input.split('/')[-1].split('.')[0]
		input_df = delayed_input(input,cutoff,base)
		results.append(input_df)
	return results
	
@dask.delayed
def delayed_input(input,cutoff,base):
	x = dd.read_csv(input,header=None)
	x = x.dropna(how='any')
	#,dtype={0:np.int64,1:'str'})
	x = x.set_index(x.columns[1])
	x = x[x[0]>=cutoff]
	x.columns = [base]
	return x

def intersection(indict):
	dfs = [indict[x] for x in indict.keys()]
	for x in range(0,len(dfs)):
		dfs[x].columns = [list(indict.keys())[x]]
		print(dfs[x].npartitions)
	mdf = dd.concat(dfs,axis=1,join='inner',interleave_partitions=True)
	return mdf

def combined_curation(in_df,cutoffs):
	results = []
	for cutoff in cutoffs:
		print cutoff
		results.append(delayed_cutoff(in_df,cutoff))
	return results

@dask.delayed
def delayed_cutoff(in_df,cutoff):
	col_cur = []
	for col in in_df.columns:
		col_cur.append('{0}>={1}'.format(col,cutoff))
	curation = ' & '.join(col_cur)
	out_df = in_df.query(curation)
	return out_df
	
def undelayed_cutoff(in_df,cutoff):
	col_cur = []
	for col in in_df.columns:
		col_cur.append('{0}>={1}'.format(col,cutoff))
	curation = ' & '.join(col_cur)
	out_df = in_df.query(curation)
	return out_df

def calc_perc_change(df_list):
	totals_list = []
	for x in df_list:
		totals_list.append(x.sum(axis=0))
	perc_change = []
	for x in range(1,len(totals_list)):
#		print x
		perc_change.append(1-totals_list[x]/totals_list[x-1])
	return perc_change

def plot_change(perc_change,cutoffs,outname):
	fig,ax = plt.subplots(1,1)
	for y in perc_change.columns:
		ax.scatter(cutoffs[1:],perc_change[y])
	ax.set_xlabel('Cutoff')
	ax.set_ylabel('%_change')
	ax.set_xscale('log')
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	ax.legend()
	plt.tight_layout()
	plt.savefig(outname+'_total_change.pdf')
	plt.clf()
	

def all_funcs(file_lst,cutoffs,combined,outname):
	cutoffs = [int(x) for x in cutoffs.split(',')]
	if combined == False:
		cutoff_dfs = []
		for cut in cutoffs:
			in_df_lst = dask.compute(delayed_shell(file_lst,cut))
			concat_count_df = intersection_dfs(in_df_lst[0])
			cutoff_dfs.append(concat_count_df)
	else:
#		print file_lst
		concat_count_df = read_combined_input(file_lst[0])
#		print concat_count_df.compute()
		cutoff_dfs = dask.compute(combined_curation(concat_count_df,cutoffs))[0]
#	print cutoff_dfs
	perc_change = calc_perc_change(cutoff_dfs)
	perc_change = list(dask.compute(*perc_change))
	perc_change_df = pd.DataFrame(columns=perc_change[0].index)
	perc_change_df = perc_change_df.append(perc_change,ignore_index=True)
#	for x in perc_change:
#		perc_change_df = perc_change_df.append(x,ignore_index=True)
	print perc_change_df.shape
	print len(cutoffs[1:])
	perc_change_df.index = cutoffs[1:]
	perc_change_df.to_csv('{0}_totals_perc_change.csv'.format(outname))
	plot_change(perc_change_df,cutoffs,outname)
	
		

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='script to plot changes in totals')
	parser.add_argument('--input','-i',required=True,nargs='*',
		help= 'Specify input files')
	parser.add_argument('--cutoffs','-c',required=False,
		default='1,5,10,50,100,250,500,750,1000,2500,5000,7500,10000',
		help = '''
		Dictate cutoff list separated by commas. 
		Default list is: 1,5,10,50,100,250,500,750,1000,2500,5000,7500,10000
		'''
		)
	parser.add_argument('--combined','-co',required=False,default=False,
		action='store_true',
		help='Call this flag if the input is already a combined input file. Best case is if that input file has no prior curation on read counts')
	parser.add_argument('--out','-o',required=False,default=None,
		help='Specify outfile prefix')
	
	parsed = parser.parse_args()
	if parsed.out is None:
		outname = parsed.input[0].split('/')[-1].split('.')[0]
	else:
		outname = parsed.out
	all_funcs(parsed.input,parsed.cutoffs,parsed.combined,outname)