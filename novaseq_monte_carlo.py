#Monte carlo sampling the NGS reads to recalculate fold enrichment
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from collections import OrderedDict
import argparse

def read_csv(input):
	x = pd.read_csv(input,header=None,index_col=False)
	x.columns = ['sequence']
	return x

def read_mapped_barcodes(barcode_file):
	x = pd.read_csv(barcode_file,header=0,index_col=None,dtype='str')
	x.rename(columns={x.columns[0]: 'variants'},inplace=True)
	return x

def read_control_names(cn_file):
	x = pd.read_csv(cn_file,header=None)
	x.columns=['labels','barcode']
	x.set_index(x.columns[1],inplace=True)
	return x	
	
def read_diagram(diagram_file):
	with open(diagram_file,'r') as f:
		diag_list = [line.rstrip() for line in f]
	diag_list = [x.upper() for x in diag_list]
	return diag_list

def diag_indices(diag_list,variables):
	vars = [x.upper() for x in variables]
	var_list = [diag_list.index(x) for x in vars]
	con_list = [x for x in range(0,len(diag_list)) if x not in var_list]
	return (var_list,con_list)

def get_barcodes(count_df,diag_list,indices,ignore_constant):
	out_dict=OrderedDict()
	fasta_list = []
	og_count_df = count_df.copy()
	count_df = count_df.dropna(how='any',axis=0)
	#initial curation for constant regions:
	if ignore_constant == True:
		no_con_df = og_count_df
	else:
		for con in indices[1]:
			count_df = pd.DataFrame(count_df[count_df.iloc[:,0].str.contains(
				diag_list[con],regex=False)])
			print 'constant_region_{0}: {1}'.format(con,count_df.shape)
		no_con_df = og_count_df.drop(count_df.index)
	
#This for loop cycles through the fasta list and searches for the sequence of the fasta in the count_df dataframe. It also sorts the values of the df based on count and appends them to out_dict
	print og_count_df.shape
	print 'Count_DF_initial: {0}'.format(count_df.shape)
	barcode_ndx = diag_list.index('BARCODE')
	constant_indices = [barcode_ndx-1,barcode_ndx+1]
	out_dict = OrderedDict()
	barcode_dict = OrderedDict()
	out_df = pd.DataFrame()
	temp = count_df.iloc[:,0].str.extract(
		r'{0}(.*){1}'.format(
			diag_list[constant_indices[0]],diag_list[constant_indices[1]]),expand=False)
	if 'UMI' in diag_list:
		umi_ndx = diag_list.index('UMI')
		umi_indices = [umi_ndx-1,umi_ndx+1]
		temp_umi = count_df.iloc[:,0].str.extract(
			r'{0}(.*){1}'.format(
			diag_list[umi_indices[0]],diag_list[umi_indices[1]]),expand=False)
	if 'UMI' in diag_list:
		out_df['barcode']=temp_umi+'XXXX'+temp
	else:
		out_df['barcode']='XXXX'+temp
	return out_df

def get_barcode_counts(in_df,cutoff):
	out_dict = OrderedDict()
	totals_dict = OrderedDict()
	curated_unique_barcodes = OrderedDict()
	counts = in_df.groupby(['barcode']).size().reset_index(name='count')
	umi_barcodes = pd.DataFrame(counts['barcode'].str.split('XXXX',expand=True))
#	print umi_barcodes
	umi_barcodes.rename(columns={0:'umi',1:'barcode'},inplace=True)
	umi_barcode_counts = umi_barcodes.groupby(
		['barcode']).size().reset_index(name='count')
	y = umi_barcode_counts[umi_barcode_counts['count']>cutoff]
	y.set_index(['barcode'],inplace=True)
	x = np.sum(y,axis=0)
	return y,x

def intersection(indict):
	dfs = [indict[x] for x in indict.keys()]
	for x in range(0,len(dfs)):
		dfs[x].rename(columns={'count':indict.keys()[x]},inplace=True)
	mdf = pd.concat(dfs,axis=1,join='inner')
	mdf.reset_index(inplace=True)
	return mdf

def control_intersection(bc_df,count_df):
	concat_count_df_cp = count_df.copy()
	concat_count_df_cp.set_index(concat_count_df_cp.columns[0],inplace=True)
	control_dict = {
		'count': concat_count_df_cp,
		'control': bc_df
		}
	control_data = intersection(control_dict)
	return control_data

def get_mapped_barcodes(barcode_file,fi_file):
	out_dict = OrderedDict()
	for row in range(0,barcode_file.shape[0]):
		bcs = pd.DataFrame(barcode_file.iloc[row,1:])
		bcs.dropna(how='any',axis=0,inplace=True)
		bcs.columns=['barcode']
		intersection = pd.merge(bcs,fi_file,how='inner',on=['barcode'])
		out_dict[barcode_file.iat[row,0]]=intersection
	return out_dict
	
def ratios(mdf,inner,outer,file_list,total_df):
	file_list = [x.split('/')[-1].split('.')[0] for x in file_list]
	print mdf
	print file_list
	print total_df
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
			numerator = mdf[file_list[x[0]]]*total_df.loc['count',file_list[x[1]]]
			denominator = mdf[file_list[x[1]]]*total_df.loc['count',file_list[x[0]]]
			first_ratio_df[file_list[x[0]]]=numerator/denominator
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

def sum_ratios(sum_df,inner,outer,total_df):
	inner = [x.split(',') for x in inner]
	outer = [x.split(',') for x in outer]
	in_lst = [(int(x[0]),int(x[1])) for x in inner]
	out_lst = [(int(x[0]),int(x[1])) for x in outer]
	file_list = sum_df.columns.values.tolist()
	first_ratio_df = pd.DataFrame()
	final_ratio_df = pd.DataFrame()
	if inner is not None:
		for x in in_lst:
			numerator = sum_df[file_list[x[0]]]*total_df.loc['count',file_list[x[1]]]
			denominator = sum_df[file_list[x[1]]]*total_df.loc['count',file_list[x[0]]]
			first_ratio_df[file_list[x[0]]]=numerator/denominator
	if outer is not None:
		for x in out_lst:
			second_ratio=first_ratio_df[file_list[x[0]]]/first_ratio_df[file_list[x[1]]]
			final_ratio_df[file_list[x[0]]]=second_ratio
	return sum_df,first_ratio_df,final_ratio_df

def monte_carlo_read_sampling(read_df,nsample):
	nsample = np.float64(nsample)
#	print read_df.shape[0]
	if nsample >=1.0:
		nchoice = np.int64(nsample)
	else:
		nchoice = np.int64(nsample*read_df.shape[0])
#	print nchoice
	mc_df = pd.DataFrame(np.random.choice(read_df.iloc[:,0],size=nchoice,replace=True))
	return mc_df

def all_funcs(input_lst,map_file,barcode_file,take_sum,bc_controls,outname,
		take_median,diagram_vars,cutoff,inner,outer):
	in_df_dict = OrderedDict()
	mapping_df = read_mapped_barcodes(barcode_file)
	total_dict = OrderedDict()
	variables = diagram_vars.split(',')
	diag_list = read_diagram(map_file)
	diag_ndxs = diag_indices(diag_list,variables)
	for input in input_lst:
		base = input.split('/')[-1].split('.')[0]
		input_df = read_csv(input)
		barcodes_df = get_barcodes(input_df,diag_list,diag_ndxs,False)
		in_df_dict[base],total_dict[base] = get_barcode_counts(barcodes_df,cutoff)
	concat_count_df = intersection(in_df_dict)
	if bc_controls is not None:
		control_names_df = read_control_names(bc_controls)
		control_ind_df = control_intersection(control_names_df,concat_count_df)
		control_ind_df.to_csv('{0}_control_fold_enrichment.csv'.format(outname))
	else:
		pass
	mapped_dict = get_mapped_barcodes(mapping_df,concat_count_df)
	total_df = pd.DataFrame.from_dict(total_dict,orient='columns')
	out_dict = OrderedDict()
	for x in mapped_dict:
		if take_sum == False:
			first_ratio_df,final_ratio_df = ratios(mapped_dict[x],inner,outer,input_lst,
				total_df)
			if bc_controls is not None:
				norm_df = normalize(final_ratio_df,control_ind_df)
			else:
				norm_df = final_ratio_df
			print 'norm_df: \n',norm_df
			if take_median == False:
				out_dict[x] = [norm_df.iloc[:,1].mean(axis=0),norm_df.iloc[:,1].std(axis=0)]
				colnames = ['mean','stdev']
			else:
				out_dict[x] = [norm_df.iloc[:,1].median(axis=0),norm_df.iloc[:,1].std(axis=0)]
				colnames = ['median','stdev']
		else:
			sum_df = pd.DataFrame(mapped_dict[x].sum(
				axis=0,skipna=True,numeric_only=True)).transpose()
			sum_df,first_ratio_df,final_ratio_df = sum_ratios(sum_df,inner,outer,total_df)
			out_dict[x] = [final_ratio_df.iat[0,0],0]
			colnames = ['sum_ratio','stdev']
	out_df = pd.DataFrame.from_dict(out_dict,orient='index')
	out_df.columns = colnames
	if outname is None:
		outname = input_lst[0].split('/')[-1].split('.')[0]
	out_df.to_csv('{0}_fold_enrichment.csv'.format(outname))
	total_df.to_csv('{0}_total_reads.csv'.format(outname))
	concat_count_df.to_csv('{0}_barcode_counts.csv'.format(outname))
		
def monte_carlo_funcs(input_lst,map_file,barcode_file,take_sum,bc_controls,outname,
		take_median,diagram_vars,cutoff,inner,outer,mc_sample,mc_cycles):
	variables = diagram_vars.split(',')
	diag_list = read_diagram(map_file)
	diag_ndxs = diag_indices(diag_list,variables)
	in_df_lst = [read_csv(input) for input in input_lst]
	mapping_df = read_mapped_barcodes(barcode_file)
	print mapping_df
	mc_out_df = pd.DataFrame()
	mc_std_df = pd.DataFrame()
	mc_out_df['variants'] = mapping_df['variants']
	mc_std_df['variants'] = mapping_df['variants']
	for x in range(0,mc_cycles):
		in_df_dict = OrderedDict()
		total_dict = OrderedDict()
		for input in range(0,len(input_lst)):
			base = input_lst[input].split('/')[-1].split('.')[0]
			mc_in_df = monte_carlo_read_sampling(in_df_lst[input],mc_sample)
			barcodes_df = get_barcodes(mc_in_df,diag_list,diag_ndxs,False)
			in_df_dict[base],total_dict[base] = get_barcode_counts(barcodes_df,cutoff)
		concat_count_df = intersection(in_df_dict)
		if bc_controls is not None:
			control_names_df = read_control_names(bc_controls)
			control_ind_df = control_intersection(control_names_df,concat_count_df)
		mapped_dict = get_mapped_barcodes(mapping_df,concat_count_df)
		total_df = pd.DataFrame.from_dict(total_dict,orient='columns')
		out_dict = OrderedDict()
		for y in mapped_dict:
			if take_sum == False:
				first_ratio_df,final_ratio_df = ratios(
					mapped_dict[y],inner,outer,input_lst,total_df)
				if bc_controls is not None:
					norm_df = normalize(final_ratio_df,control_ind_df)
				else:
					norm_df = final_ratio_df
				if take_median == False:
					out_dict[y] = [norm_df.iloc[:,1].mean(axis=0),norm_df.iloc[:,1].std(axis=0)]
					colnames = ['mean','stdev']
				else:
					out_dict[y] = [norm_df.iloc[:,1].median(axis=0),norm_df.iloc[:,1].std(axis=0)]
					colnames = ['median','stdev']
			else:
				sum_df = pd.DataFrame(mapped_dict[y].sum(
					axis=0,skipna=True,numeric_only=True)).transpose()
				sum_df,first_ratio_df,final_ratio_df = sum_ratios(
					sum_df,inner,outer,total_df)
				out_dict[y] = [final_ratio_df.iat[0,0],0]
				colnames = ['sum_ratio','stdev']
		
		out_df = pd.DataFrame.from_dict(out_dict,orient='index')
#		print out_df
		out_df.reset_index(inplace=True)
		mc_out_df[x] = out_df.iloc[:,1]
#		print mc_out_df
		mc_std_df[x] = out_df.iloc[:,2]
	mc_out_df.to_csv('{0}_mc_averages.csv'.format(outname))
	mc_std_df.to_csv('{0}_mc_stds.csv'.format(outname))
	final_out_df = pd.DataFrame()
	final_out_df['variants']=mapping_df['variants']
	final_out_df['mean']=mc_out_df.iloc[:,1:].mean(axis=1)
	final_out_df['stdev']=mc_std_df.iloc[:,1:].mean(axis=1)
	final_out_df.to_csv('{0}_mc_final.csv'.format(outname))

		
if __name__ == '__main__': 
	parser = argparse.ArgumentParser(prefix_chars='-+',description='''
	Script to analyze NovaSeq data for V5 Single plasmid systems. 
	Compatible with UMI-tagged amplicons. 
	Requires the following:
		good_reads.csv output from quality_filter_and_translate.cpp
		Amplicon map with constant regions, UMI (optional), and barcode
		Barcode DF with mapped barcodes per each variant (Dict with barcodes in
			a row and variants in the first column)
	''')
	required = parser.add_argument_group('Required')
	required.add_argument('--input','-i',required=True,nargs='*',
		help='Input *_good_reads.csv file')
	required.add_argument('--map','-m',required=True,
		help='Amplicon map. Each line should have a different region of the amplicon')
	required.add_argument('--barcodes','-b',required=True,
		help='Mapped barcode file from Miseq or PacBio mapping')
	parser.add_argument('--monte_carlo','-mc',required=False,default=None,
		help='Specify Monte Carlo sample size. Default is None -- no MC analysis performed. If the value is a decimal, will take a percentage of the respective read count.')
	parser.add_argument('--sum','-s',required=False,default=False,action='store_true',
		help='Call this flag to force summed analysis of variant performance')
	parser.add_argument('--controls','-c',required=False,default=None,action='store_true',
		help='CSV file with name,barcode (no header)')
	parser.add_argument('--outname','-on',required=False,default=None,
		help='Specify outfile prefix')
	parser.add_argument('--median','-me',required=False,default=False,action='store_true',
		help='Call this flag to force median analysis instead of mean')
	parser.add_argument('--variables','-v',required=False,default='umi,barcode',
		help='Comma-separated list of variables in the map file')
	parser.add_argument('--cutoff','-cu',required=False,default=3,
		help='Read count threshold to curate good_reads.csv')
	parser.add_argument('--inner','-n',required=False,nargs='*',default=None,
		help='''
		Comma separated numbers of indices in --input that go together. These will be ratio'd first (first number is numerator, second is denominator) before the outer pairs. If there are multiple pairs, separate them by a space. Ex: '1,2 3,4'
		''')
	parser.add_argument('--outer','-o',required=False,nargs='*',default=None,
		help='''
		Comma separated numbers of indices in --input that go together. These will be ratio'd (first number is numerator, second is denominator) after the inner pairs. Please use only numbers found in the NUMERATOR of the inner set. If there are multiple pairs, separate them by a space.
		''')
	parser.add_argument('++mc_cycles','+c',required=False,default=100,
		help='Specify number of cycles for Monte Carlo sampling')
	parser.add_argument('++bc_counts','+b',required=False,default=False,
		action='store_true',
		help='Indicate a barcode_counts.csv file that can be used instead of multiple input files.')
	parser.add_argument('++shell','+s',required=False,default=False,action='store_true',
		help='Call this flag if input files are from barcode counts in the format of barcode,count from bash analysis (no header).')
		
	
	parsed = parser.parse_args()
	if parsed.outname is None:
		outname = input_lst[0].split('/')[-1].split('.')[0]
	else:
		outname = parsed.outname
	if parsed.bc_counts == False:
		if parsed.monte_carlo is None:
			all_funcs(parsed.input,parsed.map,parsed.barcodes,parsed.sum,
				parsed.controls,outname,parsed.median,parsed.variables,
				int(parsed.cutoff),parsed.inner,parsed.outer)
		else:
			monte_carlo_funcs(parsed.input,parsed.map,parsed.barcodes,parsed.sum,
				parsed.controls,outname,parsed.median,parsed.variables,
				int(parsed.cutoff),parsed.inner,parsed.outer,parsed.monte_carlo,
				int(parsed.mc_cycles))
	else:
		pass
	