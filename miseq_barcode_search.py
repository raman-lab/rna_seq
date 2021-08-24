#script to search for a 16N barcode given a provided constant flanking region.

import argparse
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import scipy.stats
import pandas as pd
from collections import OrderedDict


def search_barcode(sequence,five_con,three_con):
	five_con = five_con.upper()
	three_con = three_con.upper()
	if five_con in sequence and three_con in sequence:
		five_con_start = sequence.find(five_con)
		print five_con_start
		three_con_start = sequence.find(three_con)
#		print three_con_start
		barcode_start = five_con_start+len(five_con)
		barcode = sequence[barcode_start:three_con_start]
		if len(barcode) == 16:
			return barcode
		else:
#			print len(barcode)
			return None
	else:
		return None
		
def read_file(input,five_con,three_con):
	out_dict = OrderedDict()
	out_dict['No_barcode']=0
	counter = 0
	with open(input,'r') as f:
		for line in f:
			counter +=1
			print counter
			barcode = search_barcode(line,five_con,three_con)
			if barcode is not None:
				if barcode in out_dict.keys():
					out_dict[barcode]+=1
				else:
					out_dict[barcode]=1
			else:
				out_dict['No_barcode']+=1
	return out_dict

def plot_barcode_dist(out_df,name):
	fig,(ax1,ax2) = plt.subplots(2,1,sharex=True)
	graph1 = ax1.hist(out_df.iloc[:,0].tolist(),bins=100)
	graph2 = ax2.hist(out_df.iloc[:,0].tolist(),bins=100)
	ymax = scipy.stats.mode(out_df.iloc[:,0].tolist())[1][0]
	ax1.set_ylim(ymax*0.9,ymax)
	ax2.set_ylim(0,20)
	ax1.spines['bottom'].set_visible(False)
	ax1.spines['right'].set_visible(False)
	ax1.spines['top'].set_visible(False)
	ax2.spines['right'].set_visible(False)
	ax2.spines['top'].set_visible(False)
	plt.savefig('test.pdf',transparent=True)
	plt.close()
	

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='script to search for a 16N barcode given a provided constant flanking region.')
	parser.add_argument('--input','-i',required=True,
		help='Input read file')
	parser.add_argument('--five_con','-5',required=True,
		help='Five prime constant region around barcode')
	parser.add_argument('--three_con','-3',required=True,
		help='Three prime constant region around barcode')
	parsed = parser.parse_args()
	name=parsed.input.split('/')[-1].split('.')[0]
	out_dict = read_file(parsed.input,parsed.five_con,parsed.three_con)
	out_df = pd.DataFrame.from_dict(out_dict,orient='index')
	print out_df.iloc[:,0]
	out_df.to_csv('{0}_barcodes.csv'.format(name))
	plot_barcode_dist(out_df,name)
