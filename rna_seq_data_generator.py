import pandas as pd
import numpy as np
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.SeqIO import FastaIO
from Bio import pairwise2
from collections import OrderedDict
import argparse

def read_diagram(diagram_file):
	with open(diagram_file,'r') as f:
		diag_list = [line.rstrip() for line in f]
	diag_list = [x.upper() for x in diag_list]
	return diag_list

def read_input(input_file):
	return pd.read_csv(input_file,
		names=['variant','average','stdev'],header=None,index_col=None,
		dtype={'variant':str,'average':np.float64,'stdev':np.float64})

def diag_indices(diag_list,variables):
	vars = [x.upper() for x in variables]
	var_list = [diag_list.index(x) for x in vars]
	con_list = [x for x in range(0,len(diag_list)) if x not in var_list]
	return (var_list,con_list)

def make_barcodes(in_df,num_bar,bar_length):
	bases = ['A','G','T','C']
	out_dict = OrderedDict()
	for x in range(0,in_df.shape[0]):
		out_dict[in_df.iat[x,0]] = [''.join(y for y in np.random.choice(
			bases,bar_length)) for z in range(0,num_bar)]
	out_df = pd.DataFrame.from_dict(out_dict,orient='index')
	out_df.reset_index(inplace=True)
	out_df.rename(columns={out_df.columns[0]:'variant'},inplace=True)
	return out_df

def predict_counts(in_df,barcode_df,dna_counts,rna_basal_counts,read_stds):
	seq_shape = barcode_df.iloc[:,1:].shape
	basal_dna_counts = np.round(np.random.normal(loc=dna_counts,
		scale=dna_counts*read_stds,size=seq_shape))
#	print basal_dna_counts
	basal_dna_total = np.sum(basal_dna_counts)
#	print basal_dna_total
	print 'basal_dna_ratio: \n',basal_dna_counts/basal_dna_total
	ind_dna_counts = np.round(np.random.normal(loc=dna_counts,
		scale=dna_counts*read_stds,size=seq_shape))
	ind_dna_total = np.sum(ind_dna_counts)
	basal_rna_count_arr = np.round(np.random.normal(loc=rna_basal_counts,
		scale=rna_basal_counts*read_stds,size=seq_shape))
	basal_rna_total = np.sum(basal_rna_count_arr)
	print 'basal_rna_ratio: \n',basal_rna_count_arr/basal_rna_total
	basal_ratio = (basal_rna_count_arr/basal_rna_total)/(basal_dna_counts/basal_dna_total)
	print 'basal_ratio: \n',basal_ratio
	fe_ratio = np.random.normal(loc=in_df.iloc[:,1],scale=in_df.iloc[:,2])
	fe_ratio = fe_ratio.reshape(in_df.shape[0],1)
	print 'fe_ratio: \n',fe_ratio
	ind_ratio = fe_ratio*basal_ratio
	print 'ind_ratio: \n',ind_ratio
#	print ind_dna_counts
#	print ind_dna_total
	print 'ind_dna_ratio: \n',ind_dna_counts/ind_dna_total
	ind_rna_count_perc = ind_ratio*(ind_dna_counts/ind_dna_total)
	print "ind_rna_count_perc: \n",ind_rna_count_perc
	print ind_rna_count_perc/np.linalg.norm(ind_rna_count_perc,keepdims=True)
	#I have to normalize the ind_rna_count_perc array to its total because the sum of the ind_rna_count_perc dataframe is >1 and does not make sense given it should show the fraction of reads per barcode.
	test = ind_rna_count_perc/np.linalg.norm(ind_rna_count_perc,keepdims=True)
	ind_rna_total = np.round(np.random.normal(loc=np.sum(basal_rna_count_arr),
		scale = np.sum(basal_rna_count_arr)*0.2))
#	print ind_rna_total
	ind_rna_count_arr = np.round(test*ind_rna_total)
	print "ind_rna_count_arr: \n",ind_rna_count_arr
	return basal_dna_counts,ind_dna_counts,basal_rna_count_arr,ind_rna_count_arr,fe_ratio
	
def make_variants(umi_len,ind_count_arr,barcodes_df,diag_list):
	bases = ['A','G','T','C']
	ind_count_df = pd.DataFrame(ind_count_arr)
	barcode_loc = diag_list.index('BARCODE')
	print diag_list
	out_df = pd.DataFrame()
	for x in range(0,barcodes_df.shape[0]):
		bc_only = barcodes_df.iloc[x,1:].copy()
		print bc_only
		#Need to double check the following line works on the ind_count_arr
		bc_counts = ind_count_df.iloc[x,:]
		print bc_counts
		bc_read_df = pd.DataFrame(np.repeat(bc_only,np.int64(bc_counts)))
		print bc_read_df
		bc_read_df.reset_index(drop=True,inplace=True)
		bc_read_df.rename(columns={bc_read_df.columns[0]:0},inplace=True)
		print bc_read_df.shape
		print np.sum(bc_counts,axis=0)
		umi_lst = [''.join(np.random.choice(bases,umi_len)) for y in range(
			0,bc_read_df.shape[0])]
		umi_df = pd.DataFrame(umi_lst)
		print umi_df.shape
		temp_df = pd.DataFrame(np.repeat('',bc_read_df.shape[0]))
		for y in range(0,len(diag_list)):
			if diag_list[y]=='BARCODE':
				temp_df = temp_df+bc_read_df
			elif diag_list[y]=='UMI':
				temp_df = temp_df+umi_df
			else:
				temp_df = temp_df+pd.DataFrame(
					np.repeat(diag_list[y],bc_read_df.shape[0]))
				#print temp_df
		out_df = out_df.append(temp_df)
		print temp_df.iloc[0,:]
		print temp_df.shape
	print out_df.shape
	print np.sum(ind_count_arr)
	return out_df
	
		
		

if __name__ == '__main__':
	parser = argparse.ArgumentParser(prefix_chars='-+',description='''
		Script that creates NGS reads for analysis
		''')
	required = parser.add_argument_group('Required')
	required.add_argument('--input','-i',required=True,
		help='CSV file with variant names, fold induction values, standard deviations')
	required.add_argument('--diagram','-d',required=True,
		help='Diagram of the desired reads')
	parser.add_argument('--barcodes','-b',required=False,default=1000,
		help='Number of barcodes per variant. Can vary by up to 10% based on Gaussian.')
	parser.add_argument('--barcode_len','-bl',required=False,default=16,
		help='Barcode length')
	parser.add_argument('++dna_counts','+d',required=False,default=5000,
		help='Define DNA count mean')
	parser.add_argument('--read_stds','-r',required=False,default=0.1,
		help='Read standard deviations. Give as a decimal.')
	parser.add_argument('++rna_counts','+r',required=False,default=50000,
		help='Define RNA basal count mean')
	parser.add_argument('--umi_len','-u',required=False,default=10,
		help='UMI length')
	parser.add_argument('--outname','-o',required=False,default=None,
		help='Specify outfile prefix')
	parsed = parser.parse_args()
	if parsed.outname is None:
		outname = parsed.input.split('/')[-1].split('.')[0]
	else:
		outname = parsed.outname
	diagram_list = read_diagram(parsed.diagram)
	in_df = read_input(parsed.input)
	bar_df = make_barcodes(in_df,int(parsed.barcodes),int(parsed.barcode_len)) 
	print bar_df
	
	basal_dna_counts,ind_dna_counts,basal_rna_count_arr,ind_rna_count_arr,fe_ratio = predict_counts(
		in_df,bar_df,int(parsed.dna_counts),
		int(parsed.rna_counts),np.float64(parsed.read_stds))
	fe_ratio_df = pd.DataFrame(fe_ratio,index=bar_df.iloc[:,0])
	fe_ratio_df.to_csv('{0}_fold_expression.csv'.format(outname))
	basal_dna_count_df = pd.DataFrame(basal_dna_counts)
	basal_rna_count_df = pd.DataFrame(basal_rna_count_arr)
	ind_dna_count_df = pd.DataFrame(ind_dna_counts)
	ind_rna_count_df = pd.DataFrame(ind_rna_count_arr)
	basal_dna_count_df.to_csv('{0}_basal_dna_count.csv'.format(outname))
	basal_rna_count_df.to_csv('{0}_basal_rna_count.csv'.format(outname))
	ind_dna_count_df.to_csv('{0}_ind_dna_count.csv'.format(outname))
	ind_rna_count_df.to_csv('{0}_ind_rna_count.csv'.format(outname))
	
	basal_dna_reads = make_variants(
		int(parsed.umi_len),basal_dna_counts,bar_df,diagram_list)
	basal_rna_reads = make_variants(
		int(parsed.umi_len),basal_rna_count_arr,bar_df,diagram_list)
	ind_dna_reads = make_variants(
		int(parsed.umi_len),ind_dna_counts,bar_df,diagram_list)
	ind_rna_reads = make_variants(
		int(parsed.umi_len),ind_rna_count_arr,bar_df,diagram_list)
	out_df = pd.DataFrame()
	out_df = out_df.append([basal_dna_reads,basal_rna_reads,ind_dna_reads,ind_rna_reads])
	out_df.to_csv('{0}_good_reads_total.csv'.format(outname),
		index_label=False,index=False,header=False)
	bar_df.to_csv('{0}_barcodes.csv'.format(outname),index=False)
	basal_dna_reads.to_csv('{0}_basal_dna_reads.csv'.format(outname),
		index_label=False,index=False,header=False)
	basal_rna_reads.to_csv('{0}_basal_rna_reads.csv'.format(outname),
		index_label=False,index=False,header=False)
	ind_dna_reads.to_csv('{0}_ind_dna_reads.csv'.format(outname),
		index_label=False,index=False,header=False)
	ind_rna_reads.to_csv('{0}_ind_rna_reads.csv'.format(outname),
		index_label=False,index=False,header=False)
	
	