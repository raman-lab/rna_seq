#script to create box plots from multiple input files.
import argparse
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from collections import OrderedDict

def read_input(input):
	test = pd.read_csv(input,header=0)
	test = test.dropna(axis=0,how='any')
	test.set_index(test.columns[0],inplace=True)
	return test

def parse_input(input_list,suffix):
	out_dict = OrderedDict()
	for file in input_list:
		file_list = file.split('_{0}'.format(suffix))[0].split('_')
		mut = file_list[-1]
		out_dict[mut]=read_input(file)
	return out_dict

def box_plot(out_dict,basename,notext,use_range):
	hfont = {'fontname': 'Helvetica'}
	to_plot = [np.array(out_dict[x].iloc[:,0]) for x in out_dict]
	fig,ax = plt.subplots(1,1)
	fig.set_size_inches(8,6)
#	print plot_vals
	if use_range == False:
		ax.boxplot(to_plot,sym='.')
	else:
		ax.boxplot(to_plot,whis='range')
	if notext == False:
		ax.set_xticklabels(out_dict.keys(),fontsize=8,rotation=45,ha='right',**hfont)
	else:
		ax.tick_params(labelbottom=False)
	plt.tight_layout()
	if notext==False:
		plt.savefig('{0}_boxplot_labeled.pdf'.format(basename),transparent=True)
	else:
		plt.savefig('{0}_boxplot.pdf'.format(basename),transparent=True)
	plt.close()


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Script to make multiple input files into a boxplot')
	parser.add_argument('--input','-i',nargs='*',required=True,
		help='Input *rsq_raw files. Can take multiple')
	parser.add_argument('--title','-t',default='fold_enrichment',
		help='Title for the output boxplot')
	parser.add_argument('--notext','-n',required=False,default=False,action='store_true',
		help='Call this flag to remove text')
	parser.add_argument('--range','-r',required=False,default=False,action='store_true',
		help='Call this flag to force whiskers to the range of the data.')
	parser.add_argument('--suffix','-s',required=False,default=None,
		help='Specify suffix from files that should be removed prior to variant identification')
	parser.add_argument('--outname','-o',required=False,default=None,
		help='Outfile name')
	
	parsed = parser.parse_args()
	if parsed.outname is None:
		basename = parsed.input[0].split('/')[-1].split('.')[0]
	else:
		basename = parsed.outname
	in_dict = parse_input(parsed.input,parsed.suffix)
	box_plot(in_dict,basename,parsed.notext,parsed.range)
	