#Script to plot NGS read count averages across multiple replicates

import pandas as pd
import numpy as np
import argparse
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

def read_input(input):
	if input.lower().endswith('.xlsx'):
		test= pd.read_excel(input,header=0,index_col=None,dtype='str')
	elif input.lower().endswith('.csv'):
		test= pd.read_csv(input,header=0,index_col=None,dtype='str')
	elif input.lower().endswith('.tsv'):
		test= pd.read_csv(input,header=0,index_col=None,sep='\t',dtype='str')
	else:
		print "Use an acceptable input format"
		quit()
	for col in test.columns:
		if col == 'count':
			test[col]=pd.to_numeric(test[col])
		else:
			pass
	return test

def get_ratios(in_df):
	total=np.sum(in_df.iloc[:,2])
	in_df['ratio']=in_df.iloc[:,2]/total
	return in_df

def get_averages(in_list):
	test = pd.DataFrame()
	test['Mutants']=in_list[0].iloc[:,0]
	for x in range(0,len(in_list)):
		test[x]=in_list[x]['ratio'].copy()
	test['average']=np.mean(test,axis=1)
	test['stdev']=np.std(test,axis=1)
	return test

def plot_averages(avg_df,name):
	fig,ax=plt.subplots(1,1)
	xvals=np.arange(avg_df.shape[0])
	ax.bar(xvals+0.5,avg_df['average'],yerr=avg_df['stdev'],linewidth=1,capsize=5,ec='black')
	ax.set_xticks(xvals+0.5)
	ax.set_xlabel('Mutant')
	ax.set_ylabel('Percent Count')
	ax.set_xticklabels(avg_df['Mutants'],ha='right',rotation=45)
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	plt.suptitle('NGS_count_ratio')
	plt.tight_layout()
	plt.savefig(name+'_count_ratio.pdf')
	plt.close(fig)

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Script to plot output from 7b_miseq_match.py')
	parser.add_argument('--input','-i',required=True,nargs='*',
		help='A number of input files from separate 7b_miseq_match.py analyses')
	parsed = parser.parse_args()
	name=parsed.input[0].split('/')[-1].split('.')[0]
	input_list = [read_input(x) for x in parsed.input]
	input_list = [get_ratios(x) for x in input_list]
	in_df = get_averages(input_list)
	plot_averages(in_df,name)
	
