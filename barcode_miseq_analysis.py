#script to plot barcode analysis based on shell scripts
#assumes files are in the order of: count,barcode

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import argparse
import scipy.stats

def read_input(input):
	return pd.read_csv(input,header=0,index_col=None)

def curate_input(in_df,name,bar_len):
	barcode_df = in_df[in_df['barcode'].str.len()==bar_len]
	other_df = in_df[in_df['barcode'].str.len()!=bar_len]
	other_df['barcode_length']=other_df['barcode'].str.len()
	other_df['barcode_length']=other_df['barcode_length'].apply(pd.to_numeric)
	other_df['barcode_length'] = other_df['barcode_length'].replace(np.nan,0)
	return barcode_df,other_df
	
def plot_barcodes(barcode_df,name,percentage,length):
	fig,ax = plt.subplots(1,1)
	if percentage==False:
		ax.hist(barcode_df.iloc[:,0].tolist(),bins=100)
	else:
		total=np.float64(np.sum(barcode_df.iloc[:,0]))
		ax.hist(barcode_df.iloc[:,0].tolist(),bins=100,density=True)
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	ax.set_ylabel('Abundance')
	ax.set_xlabel('barcode_count')
	if percentage == False:
		plt.savefig('{0}_{1}N_barcode_hist.pdf'.format(name,length),transparent=True)
	else:
		plt.savefig('{0}_{1}N_barcode_density.pdf'.format(name,length),transparent=True)
	plt.close()

def plot_other(other_df,name,percentage,length):
	fig,ax = plt.subplots(1,1)
	if percentage == False:
		ax.hist(other_df.iloc[:,2].tolist(),bins=100)
	else:
		ax.hist(other_df.iloc[:,2].tolist(),bins=100,density=True)
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	ax.set_ylabel('Abundance')
	ax.set_xlabel('barcode_length')
	if percentage == False:
		plt.savefig('{0}_not{1}N_length_hist.pdf'.format(name,length),transparent=True)
	else:
		plt.savefig('{0}_not{1}N_length_density.pdf'.format(name,lenth),transparent=True)
	plt.close()

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Script to analyze barcodes based off of get_barcode_counts.sh')
	parser.add_argument('--input','-i',required=True,
		help='Input file that has barcode counts')
	parser.add_argument('--percentage','-p',required=False,default=False,
		action='store_true',help='Call this flag to plot densities instead of absolute counts')
	parser.add_argument('--other','-o',required=False,default=True,action='store_false',
		help='Call this flag to NOT plot histograms of lengths other than 16')
	parser.add_argument('--length','-l',required=False,default=16,
		help='Define barcode length. Default is 16')
	parsed = parser.parse_args()
	name=parsed.input.split('/')[-1].split('.')[0]
	bar_len = int(parsed.length)
	in_df = read_input(parsed.input)
	barcode_df,other_df = curate_input(in_df,name,bar_len)
	barcode_df.to_csv('{0}_{1}N_barcodes.csv'.format(name,bar_len))
	other_df.to_csv('{0}_not{1}N_barcodes.csv'.format(name,bar_len))
	plot_barcodes(barcode_df,name,parsed.percentage,bar_len)
	if parsed.other == True:
		plot_other(other_df,name,parsed.percentage,bar_len)
	else:
		pass