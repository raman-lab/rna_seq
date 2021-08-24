'''
Script to match Miseq_spacer reads. Modified from pacbio_match.py script. Requires two given constant regions surrounding the barcode and a fasta file of desired gene variants (Can be partial sequences).
CSV inputs are just a list of all possible sequences

'''
import argparse
import matplotlib
matplotlib.use('agg')
import pandas as pd
import numpy as np
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.SeqIO import FastaIO
from Bio import pairwise2
from collections import OrderedDict

def read_fasta(fasta_file):
	the_fastas = SeqIO.parse(fasta_file,'fasta')
	return the_fastas

def read_csv(input):
	x = pd.read_csv(input,header=None,index_col=False)
	x.columns = ['sequence']
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

def dna_match_fastas(count_df,the_fastas,diag_list,indices,ignore_constant):
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
	final_out_dict = OrderedDict()
	for rec in the_fastas:
		fasta_list.append(rec)
		print rec.name
		out_df=pd.DataFrame(
			count_df[count_df.iloc[:,0].str.contains(
				str(rec.seq.upper()),regex=False)].copy())
		count_df = pd.DataFrame(
			count_df[~count_df.iloc[:,0].str.contains(
				str(rec.seq.upper()),regex=False)].copy())
		print 'Count_DF_final: {0}'.format(count_df.shape)
		print 'Out_DF: {0}'.format(out_df.shape)
#		out_df.drop_duplicates(subset=['sequence'],keep='first')
#		print 'Out_DF_drop_duplicates: {0}'.format(out_df.shape)
		if out_df.shape[0]>0:
			out_dict[rec.name]=out_df
		else:
			pass
		final_out_dict[rec.name]=[rec.seq,out_df.shape[0]]
#This for loop cycles through out dict, matches them to fasta_list, and then searches for the full 5'_leader+target_sequence+3'_ender
	total_reads = og_count_df.shape[0]
	total_muts = np.sum([final_out_dict[k][1] for k in final_out_dict.keys()])
	final_out_dict['total_oligos'] = [None,total_muts]
	final_out_dict['total_reads']=[None,total_reads]
	final_out_dict['no_constants']=[None,no_con_df.shape[0]]
#	final_out_dict['Unmatched']=[None,count_df.shape[0]-total_muts]
	return out_dict,final_out_dict,count_df
	

def get_barcodes(in_dict,diag_list):
	barcode_ndx = diag_list.index('BARCODE')
	constant_indices = [barcode_ndx-1,barcode_ndx+1]
	out_dict = OrderedDict()
	barcode_dict = OrderedDict()
	for x in in_dict:
		print x
		temp = in_dict[x].iloc[:,0].str.extract(
			r'{0}(.*){1}'.format(diag_list[constant_indices[0]],diag_list[constant_indices[1]]),expand=False)
		barcode_dict[x] = temp
		out_dict[x] = pd.unique(temp).tolist()
	#This for loop searches for the two constant regions around the barcode and pulls anything between the constant regions into out_dict
	#Out_dict contains only unique barcodes while barcode_dict contains all barcodes for a particular variant
	return out_dict,barcode_dict


def get_barcode_counts(barcode_dict,cutoff,unique_barcodes,strict):
	out_dict = OrderedDict()
	totals_dict = OrderedDict()
	curated_unique_barcodes = OrderedDict()
	strict_out_dict = OrderedDict()
	for k in barcode_dict:
		if strict == False:
			barcode_dict[k] = pd.DataFrame(barcode_dict[k])
			barcode_dict[k].columns = ['barcode']
		else:
			pass
		y = barcode_dict[k].groupby(
			[barcode_dict[k].columns[0]]).size().reset_index(name='count')
		y = y[['count','barcode']]
		y = y[y[y.columns[0]]>=cutoff]
		out_dict[k] = y
		print "{0} unique barcodes: {1}".format(k,y.shape[0])
		temp = pd.DataFrame(unique_barcodes[k],columns=['barcode'])
#		temp.columns=['barcode']
		xx = pd.merge(
			temp,y,how='inner',on=['barcode'])
		curated_unique_barcodes[k] = xx['barcode'].tolist()
		print len(curated_unique_barcodes[k])
		totals_dict[k]=[y.shape[0],np.sum(y[y.columns[0]])]
	if strict == True:
		#Out_dict contains all count,barcode greater than the cutoff
		#totals_dict is the total reads of barcodes that pass read cutoff
		#Curated_unique_barcodes is a dictionary whose keys are the variants and whose values are a list of barcodes
		return out_dict,totals_dict,curated_unique_barcodes
	else:
		for k in out_dict:
			out_dict[k].drop(labels=['count'],axis=1,inplace=True)
	#		print out_dict[k]
		for k in barcode_dict:
			strict_out_dict[k] = pd.merge(
				barcode_dict[k],out_dict[k],how='inner',on=['barcode'])
		return strict_out_dict,curated_unique_barcodes

#Isolates barcodes found mapped to only a single variant.
def get_unique_barcodes(in_dict,full_barcode_dict):
	x = pd.DataFrame([])
	out_dict = OrderedDict()
	unique_bars = OrderedDict()
	for l in in_dict.keys():
		if len(in_dict[l])>0:
			x = x.append(in_dict[l])
		else:
			pass
	x.columns = ['barcode']
	print 'total_barcodes: ',x.shape[0]
	#Total barcodes are all barcodes found in all variants
	y = x.groupby([x.columns[0]]).size().reset_index(name='count')
	print 'unique_barcodes: ',y.shape[0]
	#Unique barcodes across all variants
	singles = y[y['count']==1]
	print 'singles: ',singles.shape[0]
	#Singles are unique barcodes found in only one variant
	for k in in_dict.keys():
		print k
		if len(in_dict[k])>0:
			temp = pd.DataFrame(in_dict[k])
			temp.columns=['barcode']
			#zz is the intersection of the singles list and the unique barcodes in in_dict that map to a particular variant
			zz = pd.merge(singles,temp,how='inner',on=['barcode'])
#			print zz
			print singles.shape[0]
			out_dict[k]=zz['barcode'].tolist()
			aa = pd.DataFrame(full_barcode_dict[k])
	#		print full_barcode_dict[k]
			aa.columns = ['barcode']
#			print 'aa'
#			print aa
			#bb is the intersection of the singles list and all barcodes in the full_barcode_dict
			bb = pd.merge(singles,aa,how='inner',on=['barcode'])
			print bb.shape
			unique_bars[k] = bb
		else:
			pass
	duplicated_bars = y[y['count']>1]
	duplicated_bars.columns=['barcode','shared_variants']
	print('-----')
	dup_count_df = pd.DataFrame(index=duplicated_bars['barcode'].tolist())
	for k in full_barcode_dict.keys():
#		print full_barcode_dict[k]
		temp = pd.DataFrame(full_barcode_dict[k])
		temp.columns = ['barcode']
#		print temp
		cc = pd.merge(temp,duplicated_bars,how='inner',on=['barcode'])
		dd = cc.groupby([cc.columns[0]]).size().reset_index(name=k)
		dd.set_index(['barcode'],inplace=True)
		dup_count_df = dup_count_df.join(dd,how='outer')
#	print 'test'
	dup_sum_df = pd.DataFrame(dup_count_df.sum(axis=1,skipna=True),columns=['count'])
	duplicated_bars.set_index(['barcode'],inplace=True)
	duplicated_bars = duplicated_bars.join(dup_sum_df,how='outer')
	'''
	#out_dict contains the unique barcodes associated with a variant
	#unique_bars contains all barcodes that are singles that map to a variant (multiple counts of a barcode should be represented multiple times)
	'''
	return out_dict,duplicated_bars,unique_bars,dup_count_df

def determine_ratio_curation(barcode_count_dict,
		unique_barcodes,dup_count_df,ratio_float,duplicated_bars):
	max_df = dup_count_df.copy().max(axis=1,numeric_only=True)
	div_df = dup_count_df.copy().div(max_df,axis=0)
	true_df = div_df>ratio_float
	true_sum = pd.DataFrame(true_df.sum(axis=1))
	passing_bars = true_sum[true_sum[0]==1]
	totals_dict = OrderedDict()
#	print passing_bars
	if passing_bars.shape[0]>0:
		passing_bars = passing_bars.drop(0,axis=1)
		temp = passing_bars.merge(dup_count_df,left_index=True,right_index=True)
		temp.fillna(0,inplace=True)
		dup_count_df = dup_count_df[~dup_count_df.index.isin(temp.index)]
		a = pd.DataFrame(temp.columns[np.argmax(temp.values,axis=1)])
		#a is now a list of the name of the column with the max value.
		temp.reset_index(inplace=True)
		temp.rename(columns={'index':'barcode'},inplace=True)
		
		dup_sum_df = pd.DataFrame(dup_count_df.sum(axis=1,skipna=True),columns=['count'])
		duplicated_bars.drop(['count'],axis=1,inplace=True)
		duplicated_bars = duplicated_bars.join(dup_sum_df,how='inner')
		for k in barcode_count_dict:
			print k
			if k in a[0].values:
				print '{0} in a.values'.format(k)
				temp_mut_ndxs = a[a[0]==k]
#				print temp_mut_ndxs
				temp_bar_count_df = temp.iloc[temp_mut_ndxs.index.values,[temp.columns.values.tolist().index(
					k),temp.columns.values.tolist().index('barcode')]]
				temp_bar_count_df.columns = ['count','barcode']
				print temp_bar_count_df.shape
				barcode_count_dict[k] = barcode_count_dict[k].append(
					temp_bar_count_df,ignore_index=True)
				print barcode_count_dict[k].shape
				print len(unique_barcodes[k])
				unique_barcodes[k].extend(temp_bar_count_df['barcode'].tolist())
				print len(unique_barcodes[k])
			else:
				pass
			totals_dict[k] = [barcode_count_dict[k].shape[0],np.sum(
				barcode_count_dict[k][barcode_count_dict[k].columns[0]])]
	else:
		for k in barcode_count_dict:
			totals_dict[k] = [barcode_count_dict[k].shape[0],np.sum(
				barcode_count_dict[k][barcode_count_dict[k].columns[0]])]
#	print barcode_count_dict
	return barcode_count_dict,totals_dict,unique_barcodes,dup_count_df,duplicated_bars
		
		
		
		
	
	
	

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='''
		Script to match barcodes and variants from miseq_spacer reads. 
		Requires two given constant regions surrounding the
		barcode and a fasta file of desired gene variants (Can be partial sequences).
		CSV inputs are just a list of all possible sequences
		'''
		)
	required = parser.add_argument_group('Required')
	required.add_argument('--input','-i',required=True,
		help='Input CCS-derived sequence file')
	required.add_argument('--fasta','-f',required=True,
		help='Fasta file to match sequences to')
	required.add_argument('--diagram','-d',required=True,
		help='Nucleotide/amino acid sequences with the location of a barcode and 1 or more variable region(s). Must include a barcode. Each region should be separated by a new line')
	parser.add_argument('--variables','-v',required=False,default='barcode,oligo',
		help='Comma-separated list of variable regions in the diagram file')
	parser.add_argument('--out','-o',required=False,default=None,
		help='Name of output files')
	parser.add_argument('--ignore_constant','-m',required=False,action='store_true',
		default=False,
		help='Call this flag to ingore the constant region and diagram files')
	parser.add_argument('--cutoff','-c',required=False,default=1,
		help='Specify minimum number of reads required to have a barcode attributed to a variant. This cutoff is applied after all filtering.')
	parser.add_argument('--strict','-s',required=False,action='store_true',
		default=False,
		help='Call this flag to force barcodes of ANY count to be unique to a single variant. Essentially the cutoff would be applied after the unique filter.')
	parser.add_argument('--ratio','-r',default=None,required=False,
		help='''
		Provide a float. If ratio of variant/max_variant for a shared barcode is LESS
		THAN the provided float for ALL variants, then the barcode will be accepted as
		specific for a single variant
		''')
		
	parsed = parser.parse_args()
	if parsed.out is None:
		outname = parsed.input.split('/')[-1].split('.')[0]
	else:
		outname = parsed.out
	fastas = read_fasta(parsed.fasta)
	seq_df = read_csv(parsed.input)
	diag_list = read_diagram(parsed.diagram)
	print diag_list
	indices = diag_indices(diag_list,parsed.variables.split(','))
	print indices
	fasta_seq_dict,final_out_dict,final_count_df = dna_match_fastas(
		seq_df,fastas,diag_list,indices,parsed.ignore_constant)
	barcode_dict,full_barcode_dict = get_barcodes(fasta_seq_dict,diag_list)
	if parsed.strict == False:
		print 'Strict filtering not imposed. Applying read count threshold before uniqueness'
		full_barcode_dict,barcode_dict = get_barcode_counts(
			full_barcode_dict,int(parsed.cutoff),
			barcode_dict,parsed.strict)
#		print barcode_dict
		print '------'
	else:
		pass
	unique_barcode_dict,duplicated_bars,unique_barcode_count_dict,dup_count_df = get_unique_barcodes(
		barcode_dict,full_barcode_dict)
	counted_bc_dict,counted_bc_totals,curated_unique_barcodes = get_barcode_counts(
		unique_barcode_count_dict,int(parsed.cutoff),unique_barcode_dict,True)
	if parsed.ratio is not None:
		print "Applying ratio cutoff"
		counted_bc_dict,counted_bc_totals,curated_unique_barcodes,dup_count_df,duplicated_bars = determine_ratio_curation(
		counted_bc_dict,curated_unique_barcodes,dup_count_df,np.float64(parsed.ratio),
		duplicated_bars)
	else:
		pass
#	Write output
	counted_bc_totals['duplicated_read_count_total']=[duplicated_bars.shape[0],np.sum(duplicated_bars['count'])]
	final_out_df = pd.DataFrame.from_dict(final_out_dict,orient='index')
	final_out_df.columns=['sequence','count']
	final_out_df.to_csv('{0}_curated_counts.csv'.format(outname))
	unique_barcode_df = pd.DataFrame.from_dict(curated_unique_barcodes,orient='index')
	unique_barcode_df.to_csv('{0}_barcodes.csv'.format(outname))
	counted_bc_totals_df = pd.DataFrame.from_dict(counted_bc_totals,orient='index')
	counted_bc_totals_df.columns=['num_barcodes','totals']
#	print duplicated_bars
	duplicated_bars.to_csv('{0}_duplicated_barcodes.csv'.format(outname))
	dup_count_df.to_csv('{0}_duplicated_barcode_variants.csv'.format(outname),na_rep='0')
	counted_bc_totals_df.to_csv('{0}_barcode_counts_totals.csv'.format(outname))
	for x in counted_bc_dict:
		counted_bc_dict[x].to_csv('{0}_{1}_barcode_counts.csv'.format(outname,x),
			index=False)
