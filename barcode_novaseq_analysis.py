#Script to analyze novaseq data.

import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from collections import OrderedDict

def read_input(input,cutoff):
	x = pd.read_csv(input)
	if x.shape[0]==3:
		x.drop(x.columns[0],axis=1,inplace=True)
	else:
		pass
	x.set_index(x.columns[1],inplace=True)
	x.columns = [input.split('/')[-1].split('.')[0]]
	x.astype('float64').dtypes
	print [np.float64(np.sum(x.iloc[:,0]))]
	print x.shape
	y = x[x[x.columns[0]]>=cutoff]
	print y.shape
	total = [np.float64(np.sum(y.iloc[:,0]))]
	print total
	return y,total

def read_control_names(cn_file):
	x = pd.read_csv(cn_file,header=None)
	x.columns=['labels','barcode']
	x.set_index(x.columns[1],inplace=True)
	return x
	
def intersection(indict):
	dfs = [indict[x] for x in indict.keys()]
	mdf = pd.concat(dfs,axis=1,join='inner')
	mdf.reset_index(inplace=True)
	return mdf


def ratios(mdf,inner,outer,file_list,total_dict):
	file_list = [x.split('/')[-1].split('.')[0] for x in file_list]
	if inner is not None:
		inner = [x.split(',') for x in inner]
		in_lst = [(int(x[0]),int(x[1])) for x in inner]
	if outer is not None:
		outer = [x.split(',') for x in outer]
		out_lst = [(int(x[0]),int(x[1])) for x in outer]
	first_ratio_df = pd.DataFrame()
	final_ratio_df = pd.DataFrame()
	first_ratio_df['barcode'] = mdf['barcode']
	final_ratio_df['barcode'] = mdf['barcode']
	if inner is not None:
		for x in in_lst:
			first_ratio = sum_df[file_list[x[0]]]/sum_df[file_list[x[1]]]
			total_ratio = total_df[file_list[x[1]]]/total_df[file_list[x[0]]]
			first_ratio_df[file_list[x[0]]]=first_ratio*total_ratio[0]
	if outer is not None:
		for x in out_lst:
			second_ratio=first_ratio_df[file_list[x[0]]]/first_ratio_df[file_list[x[1]]]
			final_ratio_df[file_list[x[0]]]=second_ratio
	return first_ratio_df,final_ratio_df

def normalize(sample_df, control_ratio_df):
	final_avg = np.mean(control_ratio_df.iloc[:,2])
	final_std = np.std(control_ratio_df.iloc[:,2])
	normalized_sample_df = pd.DataFrame()
	normalized_sample_df['barcode']=sample_df['barcode']
	for x in sample_df.columns.values:
		if x != 'barcode':
			normalized_sample_df[x]=sample_df[x]/final_avg
		else:
			pass
	return normalized_sample_df

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Script to analyze novaseq data')
	required = parser.add_argument_group('Required')
	required.add_argument('--input','-i',required=True,nargs='*',
		help='1 or more input files for analysis')
	parser.add_argument('--inner','-n',required=False,nargs='*',default=None,
		help='''
		Comma separated numbers of indices in --input that go together. These will be ratio'd first (first number is numerator, second is denominator) before the outer pairs. If there are multiple pairs, separate them by a space. Ex: '1,2 3,4'
		''')
	parser.add_argument('--outer','-o',required=False,nargs='*',default=None,
		help='''
		Comma separated numbers of indices in --input that go together. These will be ratio'd (first number is numerator, second is denominator) after the inner pairs. Please use only numbers found in the NUMERATOR of the inner set. If there are multiple pairs, separate them by a space.
		''')
	parser.add_argument('--control_names','-c',required=False,default=None,
		help='CSV file with name,barcode (no header)')
	parser.add_argument('--outname', '-on',required=False,default=None,
		help='Outfile name')
	parser.add_argument('--cutoff','-t',required=False,default=1,
		help='Minimum number of times a sample needs to be observed to be counted.')
	
	
	parsed = parser.parse_args()
	
	in_dict = OrderedDict()
	total_dict = OrderedDict()
	if parsed.outname is None:
		outname = parsed.input[0].split('/')[-1].split('.')[0]
	else:
		outname = parsed.outname
	for file in parsed.input:
		base = file.split('/')[-1].split('.')[0]
		in_dict[base],total_dict[base]=read_input(file,int(parsed.cutoff))
	x = intersection(in_dict)
	first_ratio_df,final_ratio_df = ratios(x,parsed.inner,
			parsed.outer,parsed.input,total_dict)
	x.to_csv('{0}_raw_counts.csv'.format(outname),index=False)
	print total_dict.keys()
	total_df = pd.DataFrame.from_dict(total_dict,orient='columns')
	total_df.to_csv('{0}_totals.csv'.format(outname),index=False)
	if parsed.inner is not None:
		first_ratio_df.to_csv('{0}_first_ratio.csv'.format(outname),index=False)
		if parsed.outer is not None:
			final_ratio_df.to_csv('{0}_final_ratio.csv'.format(outname),index=False)
	if parsed.control_names is not None:
		control_names = read_control_names(parsed.control_names)
		final_ratio_df.set_index(final_ratio_df.columns[0],inplace=True)
		temp_dict = {
			'ratios': final_ratio_df,
			'control_names': control_names
			}
		control_data = intersection(temp_dict)
		if control_data.shape[0]!=control_names.shape[0]:
			print "Some controls missing"
			quit()
		else:
			pass
		control_data.to_csv('{0}_control_data.csv'.format(outname),index=False)
		if parsed.outer is not None:
			final_ratio_df.reset_index(inplace=True)
			normalized_ratio_df = normalize(final_ratio_df,control_data)
			normalized_ratio_df.to_csv(
				'{0}_normalized_ratio.csv'.format(outname),index=False)
	else:
		pass
		