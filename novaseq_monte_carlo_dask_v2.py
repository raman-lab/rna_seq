#Monte carlo sampling the NGS reads to recalculate fold enrichment.
#Uses python3 for the latest version of dask

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from collections import OrderedDict
import dask
import dask.dataframe as dd
import argparse

def read_csv(input,base):
	x = dd.read_csv(input,header=None)
	x.columns = [base]
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

def read_concat_barcodes(file_lst):
	x = dd.read_csv(file_lst[0],header=0)
	x = x.set_index(x.columns[0])
	y = pd.read_csv(file_lst[1],header=0,index_col=0)
	return x,y


def delayed_shell(input_list,cutoff):
	results = []
	for input in input_list:
		base = input.split('/')[-1].split('.')[0]
		input_df = delayed_input(input,cutoff,base)
		results.append(input_df)
	return results
	
@dask.delayed
def delayed_input(input,cutoff,base):
	x = dd.read_csv(input,header=None)
	x = x.dropna(how='any')
	#,dtype={0:np.int64,1:'str'})
	x = x.set_index(x.columns[1])
	x = x[x[0]>=cutoff]
	x.columns = [base]
	return x
	
def get_barcodes(count_df,diag_list,indices,ignore_constant):
	out_dict=OrderedDict()
	fasta_list = []
	og_count_df = count_df.copy()
	count_df = count_df.dropna(how='any')
	#initial curation for constant regions:
	if ignore_constant == True:
		no_con_df = og_count_df
	else:
		for con in indices[1]:
			count_df = count_df[count_df.iloc[:,0].str.contains(
				diag_list[con],regex=False)]
			#this print step might be slow
			print('constant_region_{0}: {1}'.format(con,count_df.compute().shape))
		#Drop does not work for indices, only axis=1
#		no_con_df = og_count_df.drop(count_df.index)
	
#This for loop cycles through the fasta list and searches for the sequence of the fasta in the count_df dataframe. It also sorts the values of the df based on count and appends them to out_dict
#	print('Count_DF_initial: {0}'.format(count_df.shape))
	barcode_ndx = diag_list.index('BARCODE')
	constant_indices = [barcode_ndx-1,barcode_ndx+1]
	out_dict = OrderedDict()
	barcode_dict = OrderedDict()
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
		out_df=temp_umi+'XXXX'+temp
	else:
		out_df='XXXX'+temp
	out_df = out_df.to_frame()
	out_df.columns = ['barcode']
	return out_df

def get_barcode_counts(in_df,cutoff):
	out_dict = OrderedDict()
	totals_dict = OrderedDict()
	curated_unique_barcodes = OrderedDict()
#	print(in_df.compute())
	counts = in_df.groupby(['barcode']).size().reset_index()
	counts.columns = ['barcode','count']
	umi_barcodes = counts['barcode'].str.partition('XXXX')
	umi_barcodes = umi_barcodes.drop([1],axis=1)
	umi_barcodes.columns = ['umi','barcode']
#	print(umi_barcodes)
	umi_barcode_counts = umi_barcodes.groupby(
		['barcode']).size().reset_index()
	y = umi_barcode_counts[umi_barcode_counts[0]>=cutoff]
	y = y.set_index('barcode')
	x = np.sum(y,axis=0).compute()
	return y,x

def intersection(indict):
	dfs = [indict[x] for x in indict.keys()]
#	print(dfs[0])
	for x in range(0,len(dfs)):
		dfs[x].columns = [list(indict.keys())[x]]
		print(dfs[x].npartitions)
	mdf = dd.concat(dfs,axis=1,join='inner',interleave_partitions=True)
	#this line is required for some reason
#	mdf = mdf.repartition(npartitions=1)
#	mdf = dd.from_pandas(mdf.compute(),npartitions=1)
	return mdf

def intersection_dfs(dfs):
	mdf = dd.concat(dfs,axis=1,join='inner')
#	mdf = dd.from_pandas(mdf.compute(),npartitions=1)
	return mdf

def control_intersection(bc_df,count_df):
#	concat_count_df_cp = count_df.set_index(count_df.columns[0])
	print(bc_df)
#	print(count_df.compute())
	control_data = dd.concat([bc_df,count_df],axis=1,join='inner')
	control_data = control_data.compute()
	return control_data

def delayed_intersection(dfs):
	mdf = dd.concat(dfs,axis=1,join='inner')
	return mdf

def get_mapped_barcodes(barcode_file,fi_file):
	out_dict = OrderedDict()
	for row in range(0,barcode_file.shape[0]):
		bcs = pd.DataFrame(barcode_file.iloc[row,1:])
		bcs.dropna(how='any',axis=0,inplace=True)
#		print(bcs)
		bcs = bcs.drop_duplicates()
		bcs.set_index(bcs.columns[0],inplace=True)
#		print(barcode_file.iat[row,0])
		print(bcs.shape)
		intersection = dd.concat([bcs,fi_file],axis=1,join='inner')
		out_dict[barcode_file.iat[row,0]]=intersection.compute()
		'''
		intersection = dd.merge(bcs,fi_file,left_index=True,right_index=True,how='inner')
		#for some reason merging here duplicates the rows
		intersection = intersection.compute()
		intersection.reset_index(inplace=True)
		intersection = intersection.drop_duplicates(keep='first')
		intersection = intersection.set_index(intersection.columns[0])
		out_dict[barcode_file.iat[row,0]]=intersection
		'''
	return out_dict
	
def ratios(mdf,inner,outer,file_list,total_df,start_val):
	file_list = [x for x in mdf.columns[start_val:]]
	if inner is not None:
		inner = [x.split(',') for x in inner]
		in_lst = [(int(x[0]),int(x[1])) for x in inner]
	if outer is not None:
		outer = [x.split(',') for x in outer]
		out_lst = [(int(x[0]),int(x[1])) for x in outer]
	first_ratio_df = pd.DataFrame(index=mdf.index)
	final_ratio_df = pd.DataFrame(index=mdf.index)
	if inner is not None:
		for x in in_lst:
			numerator = mdf[file_list[x[0]]]*total_df.loc[0,file_list[x[1]]]
			print(numerator)
			denominator = mdf[file_list[x[1]]]*total_df.loc[0,file_list[x[0]]]
			print(denominator)
			first_ratio_df[file_list[x[0]]]=np.float64(numerator/denominator)
	if outer is not None:
		for x in out_lst:
			second_ratio=np.float64(
				first_ratio_df[file_list[x[0]]]/first_ratio_df[file_list[x[1]]])
			final_ratio_df[file_list[x[0]]]=second_ratio
	return first_ratio_df,final_ratio_df

def normalize(sample_df, control_ratio_df):
	final_avg = np.mean(control_ratio_df.iloc[:,0])
	final_std = np.std(control_ratio_df.iloc[:,0])
	normalized_sample_df = pd.DataFrame(index=sample_df.index)
	for x in sample_df.columns.values:
		if x != 'barcode':
			normalized_sample_df[x]=sample_df[x]/final_avg
		else:
			pass
	return normalized_sample_df

def sum_ratios(sum_df,inner,outer,total_df,start_val):
	inner = [x.split(',') for x in inner]
	outer = [x.split(',') for x in outer]
	in_lst = [(int(x[0]),int(x[1])) for x in inner]
	out_lst = [(int(x[0]),int(x[1])) for x in outer]
	file_list = [x for x in sum_df.columns[start_val:]]
	first_ratio_df = pd.DataFrame(index=sum_df.index)
	final_ratio_df = pd.DataFrame(index=sum_df.index)
	if inner is not None:
		for x in in_lst:
			numerator = sum_df[file_list[x[0]]]*total_df.loc[0,file_list[x[1]]]
			denominator = sum_df[file_list[x[1]]]*total_df.loc[0,file_list[x[0]]]
			first_ratio_df[file_list[x[0]]]=numerator/denominator
	if outer is not None:
		for x in out_lst:
			second_ratio=first_ratio_df[file_list[x[0]]]/first_ratio_df[file_list[x[1]]]
			final_ratio_df[file_list[x[0]]]=second_ratio
	return sum_df,first_ratio_df,final_ratio_df


#----monte carlo read sampling functions
def monte_carlo_read_sampling(read_df,nsample):
	nsample = np.float64(nsample)
#	print(read_df.shape[0])
	if nsample >=1.0:
		nchoice = np.int64(nsample)
	else:
		nchoice = np.int64(nsample*read_df.shape[0])
#	print(nchoice)
	mc_df = dd.from_array(np.random.choice(read_df.iloc[:,0],size=nchoice,replace=True))
	return mc_df

#@dask.delayed
def add_col(df):
	out = df[df.columns[0]].apply(lambda x: np.random.randint(0,x))

#----monte carlo read sampling functions end


#----monte carlo barcode sampling functions
def bars_from_mapped(mapping_df):
	out_lst = []
	names = []
	for x in range(0,mapping_df.shape[0]):
		temp = dd.from_pandas(pd.DataFrame(mapping_df.iloc[x,1:]),npartitions=1)
#		print temp
		out_lst.append(temp)
		names.append(mapping_df.iat[x,0])
	return out_lst,names

def bar_rsample(dd_cols,nsample):
	samples = []
	for x in dd_cols:
		sample = bar_shuffler(x,nsample)
		samples.append(sample)
	return samples

@dask.delayed
def bar_shuffler(indf,nsample):
	shuffles = []
	for x in range(0,nsample):
		shuffles.append(indf.sample(frac=1))
	return shuffles

def bar_split_dfs(shuf_lst,nsplit):
	out_lst = []
	for dfs in shuf_lst:
		bruh = bar_partition_dfs(dfs,nsplit)
		out_lst.append(bruh)
	return out_lst

@dask.delayed
def bar_partition_dfs(dfs,nsplit):
#	nsplit = nsplit-1
	bruh = []
	for df in dfs:
		df = df.dropna(how='any')
		dfsize = int(len(df)/nsplit)
#		print(df)
#		print(dfsize)
		try:
			partitioned = dd.concat([df.iloc[i*dfsize:(i+1)*dfsize].copy() for i in range(nsplit)],axis=1)
		except ValueError:
			print('bar_partition_dfs ValueError')
			print(df)
			print(dfsize)
			quit()
		else:
			partitioned.columns = range(nsplit)
			bruh.append(partitioned)
	return bruh

def bar_get_counts(df_list,val_df):
	out = []
	for x in df_list:
		y = grab_bar_counts(dask.compute(*x),val_df)
		out.append(y)
	return out

@dask.delayed
def grab_bar_counts(df_list,val_df):
	out_dfs = []
	for x in df_list:
		temp_dfs = []
		for y in x.columns:
#			z = dd.from_pandas(x[y].to_frame(),npartitions=1)
			z = pd.DataFrame(x[y])
			t = z.dropna(how='any',axis=0)
			t = t.set_index(t.columns[0])
			try: 
				temp = dd.concat([t,val_df],join='inner',axis=1)
			except ValueError:
				print('Got valueerror')
				print(x)
				print(z)
				print(t)
#				print(val_df)
				quit()
			else:
				temp = temp.drop_duplicates(keep='first')
				temp_dfs.append(temp.sum(axis=0))
		sum_dfs = dd.concat(temp_dfs,join='outer',axis=1)
		out_dfs.append(sum_dfs)
	return out_dfs
	
#----monte carlo barcode sampling functions end


#----Bootstrap error functions------
def bootse(df):
	means = df.apply(np.mean,axis=0)
	se = np.std(means)
	return means,se

def jack_sample(arr):
	out_df = pd.DataFrame()
	for x in range(0,arr.shape[0]):
		out_df[x] = np.concatenate((arr[:x],arr[x+1:]),axis=None)
	return out_df
	
def jackse(df):
	means = np.mean(df)
	jse =  np.sqrt(np.float64((means.shape[0]-1.0)/means.shape[0])*np.sum((means-np.mean(means))**2))
	return means,jse

def bca_calc(arr,nboot,alpha,converge,threshold):
#	print arr
	if converge == False:
		theta_hat = np.nanmean(arr)
		nx = arr.shape[0]
#		bsamp = get_sample(arr,nboot)
		bsamp = pd.DataFrame(arr)
#		print bsamp.shape
		bmeans,bse = bootse(bsamp)
#		print bmeans
		jsamp = jack_sample(arr)
		jmeans,jse = jackse(jsamp)
	else:
		theta_hat = np.nanmean(arr)
		nx = arr.shape[0]
		bsamp = get_sample(arr,nboot)
		jsamp = jack_sample(arr)
		jmeans,jse = jackse(jsamp)
		while True:
			old_bmeans,old_bse = bootse(bsamp)
			bsamp = pd.concat([bsamp,get_sample(arr,nboot)],axis=1)
			bmeans,bse = bootse(bsamp)
			if np.absolute(np.mean(old_bmeans)-np.mean(bmeans))/np.mean(bmeans) < threshold:
	#			print np.absolute(old_bmeans-bmeans)/bmeans
				nboot = bsamp.shape[1]
				break
			else:
				pass
#	print bsamp.shape
#	print bse, jse
	if bse == 0 or jse == 0:
		alpha1=alpha
		alpha2=1-alpha
		z0 = None
		ahat = None
	else:
		z0 = decimal.Decimal(
			scipy.stats.norm.ppf(np.sum(bmeans<theta_hat)/np.float64(nboot))
			)
		atop = decimal.Decimal(np.sum((np.mean(jmeans)-jmeans)**3.0))
		abot = decimal.Decimal(6.0*(((jse**2.0)*(nx/(nx-1.0)))**(3/2.0)))
		ahat = atop/abot
#		print ahat,abot
		alpha1 = scipy.stats.norm.cdf(np.float64(
			z0+(z0+decimal.Decimal(scipy.stats.norm.ppf(alpha)))/(
			1-ahat*(z0+decimal.Decimal(scipy.stats.norm.ppf(alpha))))))
		alpha2 = scipy.stats.norm.cdf(np.float64(
			z0+(z0+decimal.Decimal(scipy.stats.norm.ppf(1.0-alpha)))/(
			1-ahat*(z0+decimal.Decimal(scipy.stats.norm.ppf(1-alpha))
			))))
#	print alpha1,alpha2
	confpoint = np.quantile(bmeans,[alpha1,alpha2])
	percentiles = [alpha1,alpha2]
	u=jmeans-theta_hat
	return confpoint,z0,ahat,u,bmeans,bse,bsamp,percentiles

def ci_calc(arr,alpha):
	return np.quantile(arr,[alpha,1.0-alpha])

def format_out(perc,confpoint,z0,ahat,bse,alpha):
	out_df = pd.DataFrame()
	out_df['actual_percentile']=[alpha,1.0-alpha]
	out_df['adjusted_percentiles']=perc
	out_df['confidence_interval']=confpoint
	out_df['bias_correction_z0']=[z0,None]
	out_df['acceleration_factor_ahat']=[ahat,None]
	out_df['bootstrap_standard_error']=[bse,None]
	return out_df
#----bootstrap error funcs end-----


#----visualization-----
def box_plot(out_df,basename,notext,use_range):
#	hfont = {'fontname': 'Helvetica'}
	to_plot = [out_df[x].to_list() for x in out_df.columns]
	fig,ax = plt.subplots(1,1)
	fig.set_size_inches(8,6)
#	print plot_vals
	if use_range == False:
		ax.boxplot(to_plot,sym='.')
	else:
		ax.boxplot(to_plot,whis='range')
	if notext == False:
		ax.set_xticklabels(out_df.columns,fontsize=8,rotation=45,ha='right')
	else:
		ax.tick_params(labelbottom=False)
	plt.tight_layout()
	if notext==False:
		plt.savefig('{0}_boxplot_labeled.pdf'.format(basename),transparent=True)
	else:
		plt.savefig('{0}_boxplot.pdf'.format(basename),transparent=True)
	plt.close()




#----visualization end-----
def all_funcs(input_lst,map_file,barcode_file,take_sum,bc_controls,outname,
		take_median,diagram_vars,cutoff,inner,outer,bc_input,shell):
	if outname is None:
		outname = input_lst[0].split('/')[-1].split('.')[0]
	else:
		pass
	mapping_df = read_mapped_barcodes(barcode_file)
	variables = diagram_vars.split(',')
	diag_list = read_diagram(map_file)
	diag_ndxs = diag_indices(diag_list,variables)
	if bc_input == False and shell == False:
		in_df_dict = OrderedDict()
		total_dict = OrderedDict()
		for input in input_lst:
			base = input.split('/')[-1].split('.')[0]
			input_df = read_csv(input,base)
			barcodes_df = get_barcodes(input_df,diag_list,diag_ndxs,False)
			in_df_dict[base],total_dict[base] = get_barcode_counts(barcodes_df,cutoff)
		concat_count_df = intersection(in_df_dict)
		total_df = pd.DataFrame.from_dict(total_dict,orient='columns')
	elif bc_input == True and shell == False:
		concat_count_df,total_df = read_concat_barcodes(input_lst)
	elif bc_input == False and shell == True:
		in_df_dict = OrderedDict()
		total_dict = OrderedDict()
		in_df_lst = dask.compute(delayed_shell(input_lst,cutoff))
		concat_count_df = intersection_dfs(in_df_lst[0])
		name_list = [input.split('/')[-1].split('.')[0] for input in input_lst]
		total_df = pd.DataFrame([[x.sum(axis=0).compute()[0] for x in in_df_lst[0]]],
			columns=name_list)
		'''
		for input in input_lst:
			print('importing {0}'.format(input))
			base = input.split('/')[-1].split('.')[0]
			input_df = read_shell_barcodes(input,cutoff)
			in_df_dict[base] = input_df
			total_dict[base] = np.sum(input_df,axis=0).compute()
		concat_count_df = intersection(in_df_dict)
		print(concat_count_df.npartitions)
		total_df = pd.DataFrame.from_dict(total_dict,orient='columns')
		'''
		dask.compute(concat_count_df)[0].to_csv('test.csv')
		print('wrote_concat')
#	print(concat_count_df.compute())
	if bc_controls is not None:
		control_names_df = read_control_names(bc_controls)
		control_counts_df = control_intersection(control_names_df,concat_count_df)
		control_counts_df.to_csv('{0}_control_counts.csv'.format(outname))
		control_first_ratio_df,control_ind_df = ratios(control_counts_df,inner,
			outer,input_lst,total_df,1)
		print('control\n',control_ind_df)
		control_ind_df.to_csv('{0}_control_fold_enrichment.csv'.format(outname))
	else:
		pass
	mapped_dict = get_mapped_barcodes(mapping_df,concat_count_df)
#	print(mapped_dict[mapped_dict.keys()[0]])
	out_dict = OrderedDict()
	for x in mapped_dict:
		print('calculating for {0}: '.format(x))
		if take_sum == False:
			first_ratio_df,final_ratio_df = ratios(mapped_dict[x],inner,outer,input_lst,
				total_df,0)
			if bc_controls is not None:
				norm_df = normalize(final_ratio_df,control_ind_df)
			else:
				norm_df = final_ratio_df
#			print('norm_df: \n',norm_df)
			if take_median == False:
				out_dict[x] = [norm_df.iloc[:,0].mean(axis=0),norm_df.iloc[:,0].std(axis=0)]
				colnames = ['mean','stdev']
			else:
				out_dict[x] = [norm_df.iloc[:,0].median(axis=0),norm_df.iloc[:,0].std(axis=0)]
				colnames = ['median','stdev']
		else:
			sum_df = pd.DataFrame(mapped_dict[x].sum(
				axis=0,skipna=True,numeric_only=True)).transpose()
			sum_df,first_ratio_df,final_ratio_df = sum_ratios(sum_df,inner,outer,
				total_df,0)
			out_dict[x] = [final_ratio_df.iat[0,0],0]
			colnames = ['sum_ratio','stdev']
		final_ratio_df.to_csv('{0}_{1}_ratios.csv'.format(outname,x))
		mapped_dict[x].to_csv('{0}_{1}_barcode_counts.csv'.format(outname,x))
	out_df = pd.DataFrame.from_dict(out_dict,orient='index')
	out_df.columns = colnames
	out_df.to_csv('{0}_fold_enrichment.csv'.format(outname))
	total_df.to_csv('{0}_total_reads.csv'.format(outname))
	concat_count_df = concat_count_df.compute()
	concat_count_df.to_csv('{0}_barcode_counts.csv'.format(outname))
	
def monte_carlo_funcs(input_lst,map_file,barcode_file,take_sum,bc_controls,outname,
		take_median,diagram_vars,cutoff,inner,outer,mc_cycles,mc_sample,shell):
	variables = diagram_vars.split(',')
	diag_list = read_diagram(map_file)
	diag_ndxs = diag_indices(diag_list,variables)
	in_df_lst = [read_csv(input) for input in input_lst]
	mapping_df = read_mapped_barcodes(barcode_file)
	print(mapping_df)
	mc_out_df = pd.DataFrame(index=mapping_df.index)
	mc_std_df = pd.DataFrame(index=mapping_df.index)
	if shell == False:
		for x in range(0,mc_cycles):
			in_df_dict = OrderedDict()
			total_dict = OrderedDict()
			if shell == False:
				for input in range(0,len(input_lst)):
					base = input_lst[input].split('/')[-1].split('.')[0]
					mc_in_df = monte_carlo_read_sampling(in_df_lst[input],mc_sample)
					barcodes_df = get_barcodes(mc_in_df,diag_list,diag_ndxs,False)
					in_df_dict[base],total_dict[base] = get_barcode_counts(
						barcodes_df,cutoff)
			else:
				pass
				#need to fill this in
			concat_count_df = intersection(in_df_dict)
			if bc_controls is not None:
				control_names_df = read_control_names(bc_controls)
				control_counts_df = control_intersection(control_names_df,concat_count_df)
				control_ind_df = ratios(control_counts_df,inner,outer,input_lst,total_df,1)
			mapped_dict = get_mapped_barcodes(mapping_df,concat_count_df)
			total_df = pd.DataFrame.from_dict(total_dict,orient='columns')
			out_dict = OrderedDict()
			for y in mapped_dict:
				if take_sum == False:
					first_ratio_df,final_ratio_df = ratios(
						mapped_dict[y],inner,outer,input_lst,total_df,0)
					if bc_controls is not None:
						norm_df = normalize(final_ratio_df,control_ind_df)
					else:
						norm_df = final_ratio_df
					if take_median == False:
						out_dict[y] = [norm_df.iloc[:,0].mean(axis=0),norm_df.iloc[:,0].std(axis=0)]
						colnames = ['mean','stdev']
					else:
						out_dict[y] = [norm_df.iloc[:,0].median(axis=0),norm_df.iloc[:,0].std(axis=0)]
						colnames = ['median','stdev']
				else:
					sum_df = pd.DataFrame(mapped_dict[y].sum(
						axis=0,skipna=True,numeric_only=True)).transpose()
					sum_df,first_ratio_df,final_ratio_df = sum_ratios(
						sum_df,inner,outer,total_df,0)
					out_dict[y] = [final_ratio_df.iat[0,0],0]
					colnames = ['sum_ratio','stdev']
		
		out_df = pd.DataFrame.from_dict(out_dict,orient='index')
#		print(out_df)
		out_df.reset_index(inplace=True)
		mc_out_df[x] = out_df.iloc[:,0]
#		print(mc_out_df)
		mc_std_df[x] = out_df.iloc[:,1]
	mc_out_df.to_csv('{0}_mc_averages.csv'.format(outname))
	mc_std_df.to_csv('{0}_mc_stds.csv'.format(outname))
	final_out_df = pd.DataFrame(index=mapping_df.index)
	final_out_df['variants']=mapping_df['variants']
	final_out_df['mean']=mc_out_df.iloc[:,1:].mean(axis=1)
	final_out_df['stdev']=mc_std_df.iloc[:,1:].mean(axis=1)
	final_out_df.to_csv('{0}_mc_final.csv'.format(outname))

def monte_carlo_barcode_sample(input_lst,map_file,barcode_file,take_sum,
	bc_controls,outname,take_median,diagram_vars,cutoff,inner,outer,bc_input,mc_cycles,
	mc_sample,shell):
	if outname is None:
		outname = input_lst[0].split('/')[-1].split('.')[0]
	else:
		pass
	mapping_df = read_mapped_barcodes(barcode_file)
	variables = diagram_vars.split(',')
	diag_list = read_diagram(map_file)
	diag_ndxs = diag_indices(diag_list,variables)
	if bc_input == False and shell == False:
		in_df_dict = OrderedDict()
		total_dict = OrderedDict()
		for input in input_lst:
			base = input.split('/')[-1].split('.')[0]
			input_df = read_csv(input,base)
			barcodes_df = get_barcodes(input_df,diag_list,diag_ndxs,False)
			in_df_dict[base],total_dict[base] = get_barcode_counts(barcodes_df,cutoff)
		concat_count_df = intersection(in_df_dict)
		total_df = pd.DataFrame.from_dict(total_dict,orient='columns')
	elif bc_input == True and shell == False:
		concat_count_df,total_df = read_concat_barcodes(input_lst)
#		print(concat_count_df.compute())
	elif bc_input == False and shell == True:
		in_df_dict = OrderedDict()
		total_dict = OrderedDict()
		in_df_lst = dask.compute(delayed_shell(input_lst,cutoff))
		concat_count_df = intersection_dfs(in_df_lst[0])
		print(concat_count_df.head)
		name_list = [input.split('/')[-1].split('.')[0] for input in input_lst]
		total_df = pd.DataFrame([[x.sum(axis=0).compute()[0] for x in in_df_lst[0]]],
			columns=name_list)
	if bc_controls is not None:
		print(bc_controls)
		control_names_df = read_control_names(bc_controls)
		control_counts_df = control_intersection(control_names_df,concat_count_df)
		control_counts_df.to_csv('{0}_control_counts.csv'.format(outname))
		control_first_ratio_df,control_ind_df = ratios(control_counts_df,inner,
			outer,input_lst,total_df,1)
		print('control\n',control_ind_df)
		control_ind_df.to_csv('{0}_control_fold_enrichment.csv'.format(outname))
	else:
		pass
	bar_cols,variant_names = bars_from_mapped(mapping_df)
#	print(mapping_df)
	bar_shuffled = dask.compute(*bar_rsample(bar_cols,mc_cycles))
	print(len(bar_shuffled))
	print(len(bar_shuffled[0]))
	print(bar_shuffled[0][0])
	split_shuffled = dask.compute(*bar_split_dfs(bar_shuffled,mc_sample))
	print(split_shuffled[0][0].compute())
	split_counted = dask.compute(*bar_get_counts(split_shuffled,concat_count_df))
#	zipped_muts = zip(*split_counted)
	bootstrap_dict = OrderedDict()
	boot_iterations_df = pd.DataFrame()
	for z in range(0,len(split_counted)):
		out_dict = OrderedDict()
		x = dask.compute(*split_counted[z])
		for y in range(0,len(x)):
			first_ratio_df,final_ratio_df = ratios(x[y].transpose(),inner,outer,
				input_lst,total_df,0)
			if bc_controls is not None:
				norm_df = normalize(final_ratio_df,control_ind_df)
			else:
				norm_df = final_ratio_df
#			print('norm_df: \n',norm_df)
			if take_median == False:
				out_dict[y] = [norm_df.iloc[:,0].mean(axis=0)]
				colnames = ['mean']
			else:
				out_dict[y] = [norm_df.iloc[:,0].median(axis=0)]
				colnames = ['median']
		out_df = pd.DataFrame.from_dict(out_dict,orient='index')
		confpoint,z0,ahat,u,bmeans,bse,bsamp,perc = bca_calc(
			out_df.iloc[:,0],1,
			0.05,False,0.001)
		bootstrap_dict[variant_names[z]] = [np.mean(out_df.iloc[:,0]),confpoint[0],confpoint[1]]
		boot_iterations_df[variant_names[z]] = out_df.iloc[:,0]
	bootstrap_df = pd.DataFrame.from_dict(bootstrap_dict,orient='index')
	bootstrap_df.to_csv('{0}_bootstrap_ci.csv'.format(outname))
	boot_iterations_df.to_csv('{0}_bootstrap_data.csv'.format(outname))
	box_plot(boot_iterations_df,outname,False,False)
		
			
if __name__ == '__main__': 
	parser = argparse.ArgumentParser(prefix_chars='-+',description='''
	Written for Python3
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
		help='Specify Monte Carlo sample size. Default is None -- no MC analysis performed. If the value is a decimal, will take a percentage of the respective read count. If using shell-analysis input, this must be a decimal. If barcode sampling is called, then this is the number of groups to downsample barcodes into.')
	parser.add_argument('--sum','-s',required=False,default=False,action='store_true',
		help='Call this flag to force summed analysis of variant performance')
	parser.add_argument('--controls','-c',required=False,default=None,
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
		help='Indicate a single barcode_counts.csv file that can be used instead of multiple input files. The first input must be the barcode_counts.csv file and the second must be a totals.csv file')
	parser.add_argument('++shell','+s',required=False,default=False,action='store_true',
		help='Call this flag if input files are from barcode counts in the format of barcode,count from bash analysis (no header).')
	parser.add_argument('--mc_bar','-mb',required=False,action='store_true',default=False,
		help='Call this flag to force monte carlo sampling of barcodes. If this flag is called, the --monte_carlo flag will dictate the number of barcode groups and ++mc_cycles flag will be the number of times to iterate.')
		
	
	parsed = parser.parse_args()
	if parsed.outname is None:
		outname = parsed.input[0].split('/')[-1].split('.')[0]
	else:
		outname = parsed.outname

	if parsed.monte_carlo is None and parsed.mc_bar == False:
		all_funcs(parsed.input,parsed.map,parsed.barcodes,parsed.sum,
			parsed.controls,outname,parsed.median,parsed.variables,
			int(parsed.cutoff),parsed.inner,parsed.outer,parsed.bc_counts,parsed.shell)
	elif parsed.monte_carlo is not None and parsed.mc_bar == True:
		monte_carlo_barcode_sample(parsed.input,parsed.map,parsed.barcodes,parsed.sum,
			parsed.controls,outname,parsed.median,parsed.variables,
			int(parsed.cutoff),parsed.inner,parsed.outer,
			parsed.bc_counts,int(parsed.mc_cycles),int(parsed.monte_carlo),parsed.shell)
	elif parsed.monte_carlo is not None and parsed.mc_bar == False:
		monte_carlo_funcs(parsed.input,parsed.map,parsed.barcodes,parsed.sum,
			parsed.controls,outname,parsed.median,parsed.variables,
			int(parsed.cutoff),parsed.inner,parsed.outer,parsed.monte_carlo,
			int(parsed.mc_cycles))
	else:
		pass
