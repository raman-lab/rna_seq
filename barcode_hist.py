#Python script to plot barcode counts as a histogram.

import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

def read_input(input,cutoff,need_head):
	if need_head == False:
		temp = pd.read_csv(input,header=0,index_col=None)
	else:
		temp = pd.read_csv(input,header=None,index_col=None)
		temp.columns=['count','barcode']
	temp = temp[temp['count']>=cutoff]
	return temp

def plot_input(in_df,name):
	fig,ax = plt.subplots(1,1)
	ax.hist(in_df.iloc[:,0],bins=100)
	ax.set_xlabel('Barcode_count')
	ax.set_ylabel('Abundance')
	plt.tight_layout()
	plt.savefig('{0}_hist.pdf'.format(name))
	plt.clf()

def plot_by_barcode (in_df,name):
	in_df.sort_values(by=['count'],axis=0,ascending=True,inplace=True)
	fig,ax = plt.subplots(1,1)
	xvals = np.arange(0,in_df.shape[0])
	ax.bar(xvals+0.5,in_df.iloc[:,0],width=0.8)
	ax.set_xlabel('Barcode')
	ax.set_ylabel('Read Count')
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	plt.tight_layout()
	plt.savefig('{0}_bycount.pdf'.format(name))

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Script to plot barcode counts as a histogram')
	parser.add_argument('--input','-i',required=True,
		help='Input count,barcode csv file')
	parser.add_argument('--outname','-o',required=False,default=None,
		help='Specify outfile name')
	parser.add_argument('--cutoff','-c',required=False,default=0,
		help='Minimum value to plot')
	parser.add_argument('--need_head','-n',required=False,default=False,
		action='store_true',help='Call this flag if the file needs a header')
	parsed = parser.parse_args()
	if parsed.outname == None:
		basename = parsed.input.split('/')[-1].split('.')[0]
	else:
		basename = parsed.outname
	in_df = read_input(parsed.input,int(parsed.cutoff),parsed.need_head)
	plot_input(in_df,basename)
#	plot_by_barcode(in_df,basename)