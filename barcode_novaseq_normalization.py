#Script to normalize a sample barcode dataset to a control barcode dataset
#THis is not normally used in the pacbio/novaseq workflow
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from collections import OrderedDict

def read_controls(control_file):
	x = pd.read_csv(control_file,header=None)
	x.columns=['labels','barcode']
	x.set_index(x.columns[1],inplace=True)
	return x

def read_samples(sample_file):
	return pd.read_csv(sample_file,header=0)

def read_data(input):
	x = pd.read_csv(input)
	x.drop(x.columns[0],axis=1,inplace=True)
	x.set_index(x.columns[1],inplace=True)
	x.columns = [input.split('/')[-1].split('.')[0]]
	x.astype('float64').dtypes
	return x

def get_controls(control_df,control_names):
	dfs = [control_df,control_names]
	mdf = pd.concat(dfs,axis=1,join='inner')
	mdf.drop(['labels'],axis=1,inplace=True)
	return mdf

def intersection(indict):
	dfs = [indict[x] for x in indict.keys()]
	mdf = pd.concat(dfs,axis=1,join='inner')
	mdf.reset_index(inplace=True)
	return mdf

def ratios(mdf,inner,outer,file_list,total_dict):
	file_list = [x.split('/')[-1].split('.')[0] for x in file_list]
	inner = [x.split(',') for x in inner]
	outer = [x.split(',') for x in outer]
	in_lst = [(int(x[0]),int(x[1])) for x in inner]
	out_lst = [(int(x[0]),int(x[1])) for x in outer]
	first_ratio_df = pd.DataFrame()
	final_ratio_df = pd.DataFrame()
	first_ratio_df['barcode'] = mdf['barcode']
	final_ratio_df['barcode'] = mdf['barcode']
	for x in in_lst:
		first_ratio = mdf[file_list[x[0]]]/mdf[file_list[x[1]]]
		total_ratio = total_dict[file_list[x[1]]]/total_dict[file_list[x[0]]]
		first_ratio_df[file_list[x[0]]]=first_ratio*total_ratio
	for x in out_lst:
		second_ratio=first_ratio_df[file_list[x[0]]]/first_ratio_df[file_list[x[1]]]
		final_ratio_df[file_list[x[0]]]=second_ratio
	return first_ratio_df,final_ratio_df

def normalize(sample_df, control_ratio_df):
	final_avg = np.mean(control_ratio_df.iloc[:,1])
	final_std = np.mean(control_ratio_df.iloc[:,1])
	normalized_sample_df = pd.DataFrame()
	normalized_sample_df['barcode']=sample_df['barcode']
	for x in sample_df.columns.values:
		if x != 'barcode':
			normalized_sample_df[x]=sample_df[x]/final_avg
		else:
			pass
	return normalized_sample_df
	

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Script to normalize a sample barcode dataset to a control barcode dataset')
	required = parser.add_argument_group('Required')
	required.add_argument('--input','-i',required=True,
		help='Ratio data to be normalized by the control data')
	required.add_argument('--inner','-n',required=True,nargs='*',
		help='''
		Comma separated numbers of indices in --control_data that go together. These will be ratio'd first (first number is numerator, second is denominator) before the outer pairs. If there are multiple pairs, separate them by a space. Ex: '1,2 3,4'
		''')
	required.add_argument('--control_data','-cd',required=True,nargs='*',
		help='A set of control data that will be analyzed according to --inner and --outer')
	required.add_argument('--control_names','-cn',required=True,
		help='CSV file with name,barcode (no header)')
	parser.add_argument('--outer','-o',required=False,nargs='*',default='0,2',
		help='''
		Comma separated numbers of indices in --control_data that go together. These will be ratio'd (first number is numerator, second is denominator) after the inner pairs. Please use only numbers found in the NUMERATOR of the inner set. If there are multiple pairs, separate them by a space.
		''')
	parsed = parser.parse_args()
	outname = parsed.input.split('/')[-1].split('.')[0]
	in_dict = OrderedDict()
	total_dict = OrderedDict()
	control_names = read_controls(parsed.control_names)
#	print control_names
	sample_df = read_samples(parsed.input)
	for file in parsed.control_data:
		base = file.split('/')[-1].split('.')[0]
		temp = read_data(file)
		in_dict[base]=get_controls(temp,control_names)
		#print in_dict[base]
		total_dict[base] = np.float64(np.sum(in_dict[base].iloc[:,0]))
#	print in_dict
	x = intersection(in_dict)
#	print x
	first_ratio_df,final_ratio_df = ratios(x,parsed.inner,
		parsed.outer,parsed.control_data,total_dict)
	normalized_df = normalize(sample_df,final_ratio_df)
	normalized_df.to_csv('{0}_normalized.csv'.format(outname),index=False)
	final_ratio_df.to_csv('control_data_output.csv',index=False)
	

		