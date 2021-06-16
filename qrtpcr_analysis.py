#Script to analyze the qRT-PCR output given a set of wells as input.

import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from collections import OrderedDict
import matplotlib.patches as mpatch
from matplotlib.lines import Line2D
import re

def read_map(plate_map):
	x= pd.read_csv(plate_map,header=0,sep='\t',index_col=0)
	x = x.reset_index(drop=True)
	return x

def read_input(input):
	ct_df = pd.read_excel(input,sheet_name=0,header=7,index_col=None)
	raw_df = pd.read_excel(input,sheet_name=1,header=7,index_col=None)
	return ct_df,raw_df

def read_bif_input(input):
	ct_df = pd.read_excel(input,header=19,index_col=None)
	temp_col = []
	for x in ct_df['Well'].tolist():
		y= re.findall(r"[^\W\d_]+|\d+",x)
		z = y[0]+y[1].lstrip('0')
		temp_col.append(z)
	ct_df['Well'] = temp_col
	ct_df.dropna(how='any',axis=1,inplace=True)
	return ct_df


def calc_error(std_arr):
	out_arr = np.array([])
	for col in range(0,std_arr.shape[1]):
		std_sum = np.sum(std_arr[:,col]**2)
		error = np.sqrt(std_sum)
		out_arr = np.append(out_arr,error)
	return out_arr

def dd_err(std_arr):
	std_sum = np.sum(std_arr**2)
	error = np.sqrt(std_sum)
	return error

def brep_error(std_arr):
	std_sum = np.sum(std_arr**2)
	std_sqrt = np.sqrt(std_sum)
	error=(1.0/std_arr.shape[0])*std_sqrt
	return error

def exp_error(dd_df,row):
	return np.log(2)*np.exp2(-dd_df.iat[row,0])*dd_df.iat[row,1]
	
def average_treps(plate_df,ct_df,biorad):
	samples = plate_df.iloc[:,1]
	trep_well_groups=OrderedDict()
	out_dict = OrderedDict()
	for sam in samples:
		sample_name = '||'.join(sam.split('||')[:-1])
		if sample_name not in out_dict.keys():
			out_dict[sample_name] = []
			trep_well_groups[sample_name]=[]
#	print out_dict
#	print ct_df.iloc[:,0]
	for k in out_dict:
		test=pd.DataFrame(
			plate_df[plate_df['Sample_ID'].str.contains(k,regex=False)].copy())
		
		trep_well_groups[k]=test.iloc[:,0].tolist()
		indices = []
		for well in trep_well_groups[k]:
			print well
			well_ndx, = np.where(ct_df.iloc[:,0]==well)
			indices.append(well_ndx[0])
		if biorad == False:
			trep_avg = np.mean(ct_df.iloc[indices,6])
			trep_std = np.std(ct_df.iloc[indices,6])
		else:
			trep_avg = np.mean(ct_df.iloc[indices,3])
			trep_std = np.std(ct_df.iloc[indices,3])
		out_dict[k]=[trep_avg,trep_std]
	out_df = pd.DataFrame.from_dict(out_dict,orient='index')
	out_df.columns = ['avg','stdev']
#	print trep_well_groups
	trep_well_df = pd.DataFrame.from_dict(trep_well_groups,orient='index')
	print out_df
	return out_df,trep_well_df

def average_breps(trep_df):
	trep_df['samples']=trep_df.index.values
	out_dict = OrderedDict()
	for sam in trep_df['samples']:
		sample_name = '||'.join(sam.split('||')[:-1])
		if sample_name not in out_dict.keys():
			out_dict[sample_name] = []
	for k in out_dict:
		test=trep_df[trep_df['samples'].str.contains(k,regex=False)]
		brep_avg = np.mean(test.iloc[:,0])
		brep_std = brep_error(test.iloc[:,1])
		out_dict[k] = [brep_avg,brep_std]
	out_df = pd.DataFrame.from_dict(out_dict,orient='index')
	out_df.columns=['avg','stdev']
	print out_df
	return out_df
		
def calc_ddct(brep_df,sam_list,lig_list,control):
	brep_df['samples']=brep_df.index.values
	mean_dict = OrderedDict()
	std_dict = OrderedDict()
	sam_list = sam_list.split(',')
	tot_sam = sam_list + [control]
	print tot_sam
	lig_list = lig_list.split(',')
	for sam in brep_df['samples']:
		sample_lst=sam.split('||')
		for s in tot_sam:
			base='_'.join(sample_lst[0].split('_')[:-1])
			name='{0}_{1}'.format(base,sample_lst[-1])
			if name not in mean_dict.keys():
				mean_dict[name]=[x for x in lig_list]
				std_dict[name]=[x for x in lig_list]
			else:
				pass
			if '{0}'.format(s) == sample_lst[-1]:
				ndx,=np.where(brep_df['samples']==sam)
				for l in lig_list:
					if l in sample_lst[0]:
						n = lig_list.index(l)
						mean_dict[name][n]=brep_df.iat[ndx[0],0]
						std_dict[name][n]=brep_df.iat[ndx[0],1]
					else:
						pass
			else:
				pass
#			print mean_dict
	mean_df = pd.DataFrame.from_dict(mean_dict,orient='index')
	mean_df.columns=lig_list
	std_df = pd.DataFrame.from_dict(std_dict,orient='index')
	std_df.columns=lig_list
	print mean_df
	print std_df
	mean_diff_dict = OrderedDict()
	std_diff_dict = OrderedDict()
	for ndx in mean_df.index.values:
		if control in ndx:
			control_ndx, = np.where(mean_df.index.values==ndx)
		else:
			pass
	print mean_df.index.values.tolist()
	for x in range(0,mean_df.shape[0]):
#		print x
#		print mean_df.index.values.tolist()[x]
		mean_diff_dict[mean_df.index.values.tolist()[x]]=mean_df.iloc[x,:]-mean_df.iloc[control_ndx[0],:]
		std_arr = np.vstack((std_df.iloc[x,:],std_df.iloc[control_ndx[0],:]))
#		print std_arr
		std_diff_dict[mean_df.index.values.tolist()[x]]=calc_error(std_arr)	
	mean_d_df = pd.DataFrame.from_dict(mean_diff_dict,orient='index')
	std_d_df = pd.DataFrame.from_dict(std_diff_dict,orient='index')
	mean_d_df = mean_d_df.reindex(std_d_df.index,copy=True)
	dd_dict = OrderedDict()
	print mean_d_df
	print std_d_df
	for row in range(0,mean_d_df.shape[0]):
		mean_diff = mean_d_df.iat[row,1]-mean_d_df.iat[row,0]
		std_err = dd_err(std_d_df.iloc[row,:])
		dd_dict[mean_df.index.values.tolist()[row]]=[mean_diff,std_err]
	dd_df = pd.DataFrame.from_dict(dd_dict,orient='index')
	print dd_df
	exp_df = pd.DataFrame(index=['exp2','std'])
	for row in range(0,dd_df.shape[0]):
		exp_df[dd_df.index.values[row]] = np.array([np.exp2(-dd_df.iat[row,0]),
			exp_error(dd_df,row)])
	exp_df.columns=dd_df.index.values
	print exp_df
	return dd_df,exp_df
		
def plot_treps(trep_well_df,raw_df,lig_list,samples_list,name,plot_all):
	lig_list = lig_list.split(',')
	samples_list = samples_list.split(',')
	fig,ax = plt.subplots(1,1)
	avgs_df = pd.DataFrame()
	stdevs_df = pd.DataFrame()
#	print trep_well_df.index.values
	for k in range(0,trep_well_df.shape[0]):
		out_df = pd.DataFrame()
		for y in range(0,trep_well_df.shape[1]):
			test=pd.DataFrame(raw_df[raw_df['Well']==trep_well_df.iat[k,y]].copy())
			test.reset_index(inplace=True,drop=True)
			out_df[y]=test.iloc[:,4]
		avgs_df[trep_well_df.index[k]]=np.mean(out_df,axis=1)
		stdevs_df[trep_well_df.index[k]]=np.std(out_df,axis=1)
	fig,ax = plt.subplots(1,1)
	fig.set_size_inches(8,6)
	shapes=['o','v','s','P','D','*','h']
	colors=['red','blue','green','purple','orange','yellow','pink','black']
	brep_dict = OrderedDict()
#	print avgs_df.columns.values
	for lig in lig_list:
		cols = []
		for val in avgs_df.columns.values:
#			print val
			if lig in val:
				cols.append(val)
		for col in cols:
			col_lst = col.split('||')
			if plot_all == False:
				if samples_list[0] in col_lst[1]:
					if col_lst[-1] not in brep_dict.keys():
	#					print col_lst
						brep_dict[col_lst[-1]]=[col]
					else:
						brep_dict[col_lst[-1]].append(col)
				else:
					pass
			else:
				if col_lst[-1] not in brep_dict.keys():
	#				print col_lst
					brep_dict[col_lst[-1]]=[col]
				else:
					brep_dict[col_lst[-1]].append(col)
#		print brep_dict['brep_1'][0]
	legends = []
	patches = []
	repeats = []
	mrepeats = []
	mkrs = []
	for k in range(0,len(brep_dict.keys())):
		y = brep_dict.keys()[k]
		print brep_dict[y]
#			print y
		for col in range(0,len(brep_dict[y])):
			test = brep_dict[y][col].split('||')
			lab = '{0}_{1}'.format(test[0].split('_')[-1],test[1])
			if lab not in repeats:
				repeats.append(lab)
				patches.append(mpatch.Patch(color=colors[col],label=lab))
			else:
				pass
			if test[-1] not in mrepeats:
				mrepeats.append(test[-1])
				mkrs.append(Line2D([0],[0],marker=shapes[k],label=test[-1],mfc='black',
					mec='black',ls=''))
			else:
				pass
			ax.errorbar(range(1,avgs_df.shape[0]+1),
				avgs_df[brep_dict[y][col]],fmt=shapes[k],
				yerr=stdevs_df[brep_dict[y][col]],ecolor=colors[col],
				mec=colors[col],mfc=colors[col],ms=4)
	
	markerlegend = plt.legend(handles=mkrs)
	firstlabel = plt.gca().add_artist(markerlegend)
	plt.legend(handles=patches,loc='lower right')
	ax.set_ylabel(r'$\Delta$Rn')
	ax.set_xlabel('Cycles')
	plt.suptitle(name)
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	plt.tight_layout()
	plt.savefig('{0}_dRn.pdf'.format(name),transparent=True)
	plt.close()
		
				
				
		
			
	
		
		
		



	
	


if __name__ == '__main__':
	parser=argparse.ArgumentParser(description='TBF')
	required = parser.add_argument_group('Required')
	required.add_argument('--plate_map','-m',required=True,
		help='Plate map file')
	required.add_argument('--input','-i',required=True,
		help='Output from qRTPCR run')
	parser.add_argument('--samples','-s',required=False,default='GFP',
		help='Comma separated list of samples (not control) labels as found in plate_map')
	parser.add_argument('--control','-c',required=False,default='rRNA',
		help='Endogenous control sample name')
	parser.add_argument('--ligs','-l',required=False,default='DMSO,Nar',
		help='Comma separated list of <Vehicle>,<Ligand> labels as found in plate_map')
	parser.add_argument('--plot_all','-p',required=False,default=False,
		action='store_true',help='Call this flag to plot all values')
	parser.add_argument('--biorad','-b',required=False,default=False,action='store_true',
		help='Call this flag if input file is derived from BioRad instruments. This assumes that the dRn values are not provided and no graphs will be made.')
	
	parsed = parser.parse_args()
	base=parsed.input.split('/')[-1].split('.')[0]
	map_df = read_map(parsed.plate_map)
	if parsed.biorad == False:
		ct_df,raw_df = read_input(parsed.input)
	else:
		ct_df = read_bif_input(parsed.input)
	trep_avgs_df,trep_well_df = average_treps(map_df,ct_df,parsed.biorad)
	brep_df = average_breps(trep_avgs_df)
	dd_df,exp_df = calc_ddct(brep_df,parsed.samples,parsed.ligs,
		parsed.control)
#	print raw_df
#	print trep_well_df
	
	trep_well_df.to_csv('{0}_trep_wells.csv'.format(base))
	trep_avgs_df.to_csv('{0}_trep_avgs.csv'.format(base))
	brep_df.to_csv('{0}_brep_avgs.csv'.format(base))
	dd_df.to_csv('{0}_ddct.csv'.format(base))
	exp_df.to_csv('{0}_exp2.csv'.format(base))
	if parsed.biorad == False:
		plot_treps(trep_well_df,raw_df,parsed.ligs,parsed.samples,base,parsed.plot_all)

	