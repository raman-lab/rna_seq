#python script to plot totals from miseq_spacer_map.py
import argparse
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from collections import OrderedDict

def read_totals(file):
	return pd.read_csv(file,header=0)

def plot_bars(indf,usecol,suffix,basename):
	xvals = np.arange(0,indf.shape[0]-1)
	fig,ax = plt.subplots(1,1)
	fig.set_size_inches(8,6)
	ax.bar(xvals+0.5,indf.iloc[:-1,usecol])
	ax.set_xticks(xvals+0.5)
	ax.set_xlabel('Variants')
	ax.set_ylabel(suffix)
	ax.set_xticklabels(indf.iloc[:-1,0],rotation=45,ha='right')
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	plt.tight_layout()
	plt.savefig('{0}_{1}.pdf'.format(basename,suffix))
	plt.close()

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='script to plot total file from miseq_spacer_map.py')
	parser.add_argument('--input','-i',required=True,
		help='Input file')
	parser.add_argument('--outname','-o',required=False,default=None,
		help='Specify outfile name')
	parsed = parser.parse_args()
	if parsed.outname == None:
		basename = parsed.input.split('/')[-1].split('.')[0]
	else:
		basename = parsed.outname
	indf = read_totals(parsed.input)
	plot_bars(indf,1,'numbarcodes',basename)
	plot_bars(indf,2,'totalreads',basename)