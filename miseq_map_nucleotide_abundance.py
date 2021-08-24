'''
Script to be used after miseq_spacer_map.py (*_duplicated_barcodes.csv) or on any file with a list of strings and a desire to plot the percentages of each letter at every position. Requires a target string length to be defined. Specifically designed for DNA sequences. 
'''
import argparse
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from collections import OrderedDict


def read_input(input_file):
	x = pd.read_csv(input_file,header=0,index_col=None)
	return x
	
def parse_reads(input_df,tar_len):
	out_dict = OrderedDict()
	z = []
	x = input_df[input_df[input_df.columns[0]].str.len()==tar_len]
	for l in range(0,x.shape[0]):
		z = z+[x.iat[l,0] for num in range(0,int(x.iat[l,2]))]
	y = zip(*z)
	for pos in range(0,len(y)):
		count_a = np.float64(y[pos].count('A'))/np.float64(len(y[pos]))
		count_c = np.float64(y[pos].count('C'))/np.float64(len(y[pos]))
		count_g = np.float64(y[pos].count('G'))/np.float64(len(y[pos]))
		count_t = np.float64(y[pos].count('T'))/np.float64(len(y[pos]))
		out_dict[pos]=[count_a,count_c,count_g,count_t]
	return out_dict

def plot_percentages(in_dict,name):
	xvals = [x+1 for x in in_dict.keys()]
	yvals = zip(*in_dict.values())
	labels = ['A','C','G','T']
	fig,ax = plt.subplots(1,1)
	for y in range(0,len(yvals)):
		ax.plot(xvals,yvals[y],'o-',label=labels[y])
	ax.set_xlabel('Position')
	ax.set_ylabel('Percentage')
	plt.legend()
	plt.tight_layout()
	plt.savefig('{0}_nt_perc.pdf'.format(name))
	plt.clf()

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='''
	Script to be used after miseq_spacer_map.py (*_duplicated_barcodes.csv) or on any file with a list of strings and a desire to plot the percentages of each letter at every position. Requires a target string length to be defined. Specifically designed for DNA sequences. 
	''')
	parser.add_argument('--input','-i',required=True,
		help='Input csv file with the strings in the first column')
	parser.add_argument('--length','-l',required=False,default=16,
		help='Desired length of the strings. Default is 16')
	parser.add_argument('--output','-o',required=False,default=None,
		help='Output file name. Defaults to input+Suffixes')
	parsed = parser.parse_args()
	if parsed.output == None:
		outname = parsed.input.split('/')[-1].split('.')[0]
	else:
		outname = parsed.output
	in_df = read_input(parsed.input)
	graph_dict = parse_reads(in_df,int(parsed.length))
	plot_percentages(graph_dict,outname)