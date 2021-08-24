#Script to analyze MiSeq data based on DNA oligos (not translated reads)
#Requires a Run on a CSV file generated from the following BASH script:
'''
awk '{print substr($0,9,length($0)-9)}' good_reads.csv > good_reads_trimmed.csv
sort good_reads_trimmed.csv | uniq -c | sort -nr > good_reads_trimmed_counts.csv
'''

import pandas as pd
import numpy as np
import argparse
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

def read_csv(input):
	return pd.read_csv(input)

def curate_csv(in_df,five_con,three_con,count):
	count_df = in_df[in_df['count']>=count]
	curated_df = count_df[count_df['oligo'].str.contains(
					str(five_con),regex=False)].copy()
#	print curated_df
	no_five = count_df[~count_df['oligo'].str.contains(
					str(five_con),regex=False)].copy()
	temp_curated_df = curated_df[curated_df['oligo'].str.contains(
					str(three_con),regex=False)].copy()
	temp_curated_df['oligo']=temp_curated_df['oligo'].str.extract(
		r'{0}(.*){1}'.format(five_con,three_con),expand=False)
#	print temp_curated_df
	no_three = curated_df[~curated_df['oligo'].str.contains(
					str(three_con),regex=False)].copy()
	final_curated_df = temp_curated_df.groupby(['oligo']).sum().reset_index()
	no_three_df = no_three.groupby(['oligo']).sum().reset_index()
	no_five_df = no_five.groupby(['oligo']).sum().reset_index()
	print final_curated_df.columns.values
	final_curated_df.sort_values(['count'],ascending=False,inplace=True)
	final_curated_df = final_curated_df[['count','oligo']]
	no_five_df = no_five_df[['count','oligo']]
	no_three_df = no_three_df[['count','oligo']]
	no_five_df.sort_values(['count'],ascending=False,inplace=True)
	no_three_df.sort_values(['count'],ascending=False,inplace=True)
	total = pd.DataFrame(
		[[np.sum(final_curated_df.iloc[:,0]),'total']],columns=['count','oligo'])
	return curated_df,no_five_df,no_three_df,final_curated_df,total

def plot_lengths(in_df,name):
	fig,ax = plt.subplots(1,1)
	str_len = in_df['oligo'].astype(str).map(len)
	ax.hist(str_len)
	plt.savefig('{0}_seqlen_hist.pdf'.format(name))
	plt.close()

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description=
	'''
	Script to analyze MiSeq data based on DNA oligos (not translated reads)
	Requires a Run on a CSV file generated from the following BASH script:
	awk '{print substr($0,9,length($0)-9)}' good_reads.csv > good_reads_trimmed.csv
	sort good_reads_trimmed.csv | uniq -c | sort -nr > good_reads_trimmed_counts.csv
	'''
	)
	parser.add_argument('--input','-i',required=True,help='CSV input file')
	parser.add_argument('--five','-5',required=True,
		help='5` constant region to search for')
	parser.add_argument('--three','-3',required=True,
		help="3` constant region")
	parser.add_argument('--count','-c',required=False,default=1,
		help='Minimum read count to be considered')
	parsed = parser.parse_args()
	
	in_df = read_csv(parsed.input)
	name = parsed.input.split('/')[-1].split('.')[0]
	curated_df,no_five_df,no_three_df,final_curated_df,total = curate_csv(
		in_df,parsed.five,parsed.three,int(parsed.count))
	no_five_df.to_csv('{0}_nofive.csv'.format(name),index=False)
	no_three_df.to_csv('{0}_nothree.csv'.format(name),index=False)
	final_curated_df.to_csv('{0}_curated.csv'.format(name),index=False)
	curated_df.to_csv('{0}_intermediate.csv'.format(name),index=False)
	total.to_csv('{0}_curated_total.csv'.format(name),index=False)
	plot_lengths(in_df,'{0}_orig'.format(name))
	plot_lengths(curated_df,'{0}_curated'.format(name))