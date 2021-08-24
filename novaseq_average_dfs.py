#Script to average multiple replicate sums and return average representation of each variant as a % of the total.

import argparse
import pandas as pd
import numpy as np
import decimal
from collections import OrderedDict

def read_input(input_file):
	temp = pd.read_csv(input_file,dtype='str',index_col=None)
	out_df = pd.read_csv(input_file,index_col=None)
	out_df.iloc[:,0]=temp.iloc[:,0]
#	out_df.columns = range(0,out_df.shape[1])
	out_df = out_df.rename(columns={out_df.columns[0]: 'variant'})
	return out_df

def read_tots(tot):
	x = pd.read_csv(tot)
	return x

def make_ratios(in_dict,tot_dict):
	out_dict = OrderedDict()
	for k in in_dict:
		in_dict[k].set_index(['variant'],inplace=True)
#		print in_dict[k]
		out_dict[k] = in_dict[k].div(tot_dict[k].iloc[0])
		out_dict[k].columns = range(0,out_dict[k].shape[1])
	return out_dict
	
def average_dfs(in_dict):
	panel = pd.Panel(in_dict)
	mean_df = panel.mean(axis=0)
	std_df = panel.std(axis=0)
	return mean_df, std_df

def ratios(mdf,inner,outer):
	inner = [x.split(',') for x in inner]
	in_lst = [(int(x[0]),int(x[1])) for x in inner]
	first_ratio_df = pd.DataFrame(index=mdf.index)
	final_ratio_df = pd.DataFrame(index=mdf.index)
	for x in in_lst:
		first_ratio_df[x[0]] = mdf.iloc[:,x[0]]/mdf.iloc[:,x[1]]
	if outer is not None:
		outer = [x.split(',') for x in outer]
		out_lst = [(int(x[0]),int(x[1])) for x in outer]
		for x in out_lst:
			second_ratio=first_ratio_df.iloc[:,x[0]]/first_ratio_df.iloc[:,x[1]]
			final_ratio_df[x[0]]=second_ratio
	return first_ratio_df,final_ratio_df

def propagate_error(mdf,first_avgs_df,stdevs,inner,outer):
	inner = [x.split(',') for x in inner]
	in_lst = [(int(x[0]),int(x[1])) for x in inner]

	first_ratio_df = pd.DataFrame(index=mdf.index)
	final_ratio_df = pd.DataFrame(index=mdf.index)
	for x in in_lst:
		first_ratio_df[x[0]]=np.float64(
			np.sqrt(((mdf.iloc[:,x[0]]/mdf.iloc[:,x[1]])**2)*(((stdevs.iloc[:,x[0]]/mdf.iloc[:,x[0]])**2)+((stdevs.iloc[:,x[1]]/mdf.iloc[:,x[1]])**2))))
	if outer is not None:
		outer = [x.split(',') for x in outer]
		out_lst = [(int(x[0]),int(x[1])) for x in outer]
		for x in out_lst:
			final_ratio_df[x[0]]=np.float64(
			np.sqrt(((first_avgs_df.iloc[:,x[0]]/first_avgs_df.iloc[:,x[1]])**2)*(((first_ratio_df.iloc[:,x[0]]/first_avgs_df.iloc[:,x[0]])**2)+((first_ratio_df.iloc[:,x[1]]/first_avgs_df.iloc[:,x[1]])**2))))
	else:
		pass
	return first_ratio_df,final_ratio_df

def stitch_dfs(df1,df2):
	df3 = pd.DataFrame(index=df1.index)
	for x in range(0,df1.shape[1]):
		df3['{0}_avg'.format(x)]=df1.iloc[:,x]
		df3['{0}_std'.format(x)]=df2.iloc[:,x]
	return df3
if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Script to average multiple replicate sums and return average representation of each variant as a % of the total.')
	parser.add_argument('--input','-i',required=True,nargs='*',
		help='1 or more input files that must have the same shape')
	parser.add_argument('--totals','-t',required=False,nargs='*',default=None,
		help='Total files to normalize all rows to. Omit if you just want the averages of the inputs')
	parser.add_argument('--inner','-n',required=True,nargs='*',
		help='''
		Comma separated numbers of indices in --input that go together. These will be ratio'd first (first number is numerator, second is denominator) before the outer pairs. If there are multiple pairs, separate them by a space. Ex: '1,2 3,4'
		''')
	parser.add_argument('--outer','-o',required=False,nargs='*',default=None,
		help='''
		Comma separated numbers of indices in --input that go together. These will be ratio'd (first number is numerator, second is denominator) after the inner pairs. These values must be from the result of the inner ratio dataframe. For example, '0,1' will divide the first inner pair by the second inner pair. If there are multiple pairs, separate them by a space.
		''')
	parser.add_argument('--outname','-on',required=False,
		help='Define outfile naming')
	parsed = parser.parse_args()
	if parsed.outname is None:
		outname = parsed.input[0].split('/')[-1].split('.')[0]
	else:
		outname = parsed.outname
	in_dict = OrderedDict()
	for x in range(0,len(parsed.input)):
		in_dict[x] = read_input(parsed.input[x])
	tot_dict = OrderedDict()
#	print in_dict[0]
	if parsed.totals is not None:
#		print 'totals'
		for x in range(0,len(parsed.totals)):
			tot_dict[x] = read_tots(parsed.totals[x])
		norm_dict = make_ratios(in_dict,tot_dict)
	else:
		norm_dict = in_dict
	print norm_dict[0]
	mean_df, std_df = average_dfs(norm_dict)
	mean_df.to_csv('test_mean.csv')
	std_df.to_csv('test_std.csv')
	first_mean_df,final_mean_df = ratios(mean_df,parsed.inner,parsed.outer)
	first_std_df,final_std_df = propagate_error(
		mean_df,first_mean_df,std_df,parsed.inner,parsed.outer)
	combined_first = stitch_dfs(first_mean_df,first_std_df)
	combined_first.to_csv('{0}_first_ratio.csv'.format(outname))
	combined_ratio = stitch_dfs(mean_df,std_df)
	combined_ratio.to_csv('{0}_ratio.csv'.format(outname))
	if parsed.outer is not None:
		combined_final = stitch_dfs(final_mean_df,final_std_df)
		combined_final.to_csv('{0}_final_ratio.csv'.format(outname))
	else:
		pass
	

	
	