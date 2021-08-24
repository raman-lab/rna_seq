'''
Script to match PacBio reads. Requires two given constant regions surrounding the barcode and a fasta file of desired gene variants (Can be partial sequences).
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
	x = pd.read_csv(input,header=None,index=None)
	x.columns = ['sequence']
	return x

def read_diagram(diagram_file):
	with open(diagram_file,'r') as f:
		diag_list = [line.rstrip() for line in f]
	return diag_list

def diag_indices(diag_list,variables):
	var_list = [diag_list.index(x) for x in variables]
	con_list = [x for x in range(0,len(diag_list)) if x not in var_list]
	return (var_list,con_list)

	
def dna_match_fastas(count_df,the_fastas,diag_list,indices,ignore_constant):
	out_dict=OrderedDict()
	fasta_list = []
	og_count_df = count_df.copy()
	count_df = count_df.dropna(how='any',axis=0)
	#initial curation for constant regions:
	if ignore_constant == True:
		pass
	else:
		for con in indices[1]:
			count_df = pd.DataFrame(count_df[count_df.iloc[:,1].str.contains(
				diag_list[con],regex=False)])
		no_con_df = og_count_df.drop(count_df.index)
#This for loop cycles through the fasta list and searches for the sequence of the fasta in the count_df dataframe. It also sorts the values of the df based on count and appends them to out_dict
	final_out_dict = OrderedDict()
	for rec in the_fastas:
		fasta_list.append(rec)
		out_df=pd.DataFrame(
			count_df[count_df.iloc[:,1].str.contains(
				str(rec.seq),regex=False)].copy())
		count_df = pd.DataFrame(
			count_df[~count_df.iloc[:,1].str.contains(
				str(rec.seq),regex=False)].copy())
		out_dict[rec.name]=out_df
		final_out_dict[rec.name]=[rec.seq,np.sum(out_df.iloc[:,0])]
#This for loop cycles through out dict, matches them to fasta_list, and then searches for the full 5'_leader+target_sequence+3'_ender
	total_reads = np.sum(og_count_df['sequence'])
	total_muts = np.sum([final_out_dict[k][1] for k in final_out_dict.keys()])
	final_out_dict['total_oligos'] = [None,total_muts]
	final_out_dict['total_reads']=[None,total_reads]
	final_out_dict['no_constants']=[None,no_con_df.shape[0]]
	final_out_dict['Unmatched']=[None,np.sum(count_df)-total_muts]
	return out_dict,final_out_dict,count_df

def get_barcodes(in_dict,diag_list):
	barcode_ndx = diag_list.index('barcode')
	constant_indices = [barcode_ndx-1,barcode_ndx+1]
	out_dict = OrderedDict()
	for x in in_dict:
		out_dict[x] = pd.DataFrame(in_dict[x]['sequence'].str.extract(
			r'{0}(.*){1}'.format(constant_indices[0],constant_indices[1]),expand=False)
	return out_dict

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='''
	Script to match PacBio reads. Requires two given constant regions surrounding the
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
	
	parsed = parser.parse_args()
	if parsed.out is None:
		outname = parsed.input.split('/')[-1].split('.')[0]
	else:
		outname = parsed.out
	fastas = read_fasta(parsed.fasta)
	seq_df = read_csv(parsed.input)
	diag_list = read_diagram(parsed.diagram)
	indices = diag_indices(diag_list,parsed.variables.split(','))
	fasta_seq_dict,final_out_dict,final_count_df = dna_match_fastas(
		seq_df,fastas,diag_list,indices)
	barcode_dict = get_barcodes(fasta_seq_dict,diag_list)
	final_out_df = pd.DataFrame.from_dict(final_out_dict,orient='index')
	final_out_df.columns=['sequence','count']
	final_out_df.to_csv('{0}_curated_counts.csv'.format(outname))
	barcode_df = pd.DataFrame.from_dict(barcode_dict,orient='index')
	barcode_df.to_csv('{0}_barcodes.csv'.format(outname))
	
		
		