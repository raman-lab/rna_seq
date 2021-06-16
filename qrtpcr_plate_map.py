import numpy as np
import pandas as pd
import argparse
import re


def read_input(sample_file):
	return pd.read_csv(sample_file,header=None,index_col=None)
	
def get_start_coords(start_well):
	cols=[str(x) for x in list(range(1,13))]
	rows = list('ABCDEFGH')
	start_well_list=re.split('(\d+)',start_well)
	start_well_coords = (rows.index(start_well_list[0]),cols.index(start_well_list[1]))
	return start_well_coords

def constrain_plate(xbyy,start_well_coords):
	xbyy=[int(x) for x in xbyy.split(',')]
	cols=[str(x) for x in list(range(1,13))]
	rows = list('ABCDEFGH')
	filled_rows = rows[start_well_coords[0]:start_well_coords[0]+xbyy[0]]
	filled_cols = cols[start_well_coords[1]:start_well_coords[1]+xbyy[1]]
	return filled_rows,filled_cols

def create_map(xbyy,ntrep,nbrep,vertical,random,samples,control_list,start_well):
	start_well_coords = get_start_coords(start_well)
	filled_rows,filled_cols=constrain_plate(xbyy,start_well_coords)
	total_samples = []
	total_wells = []
	for sample in samples.iloc[:,0]:
		controls = ['{0}||{1}'.format(sample,x) for x in control_list]
		for ctrl in controls:
			brep_samples = ['{0}||brep_{1}'.format(ctrl,int(x)) for x in range(1,nbrep+1)]
			for bsam in brep_samples:
				trep_samples = ['{0}||trep_{1}'.format(bsam,int(x)) for x in range(1,ntrep+1)]
				for x in trep_samples:
					total_samples.append(x)
	total_sam_arr = np.array(total_samples)
	if random == True:
		for row in filled_rows:
			for col in filled_cols:
				total_wells.append('{0}{1}'.format(row,col))
		np.random.shuffle(total_sam_arr)
		#print total_sam_arr
	else:
		if vertical==False:
			for row in filled_rows:
				for col in filled_cols:
					total_wells.append('{0}{1}'.format(row,col))
		else:
			for col in filled_cols:
				for row in filled_rows:
					total_wells.append('{0}{1}'.format(row,col))
	used_wells = [total_wells[x] for x in range(0,len(total_samples))]
	out_df = pd.DataFrame()
	out_df['Well']=used_wells
	out_df['Sample_ID']=total_sam_arr
	return out_df

if __name__ == '__main__':
	parser=argparse.ArgumentParser(description='Script to create plates for qRTPCR')
	required=parser.add_argument_group('Required')
	required.add_argument('--input','-i',required=True,
		help='Input file with sample names on separate lines')
	parser.add_argument('--xbyy','-x',required=False,default='8,12',
		help='row,col to be used. Default is entire plate (8,12)')
	parser.add_argument('--random','-r',required=False,action='store_true',default=False,
		help='Randomize wells')
	parser.add_argument('--vertical','-v',required=False,action='store_true',
		default=False,help='Call this flag to force VERTICAL orientation. Default horizontal')
	parser.add_argument('--ntrep','-t',required=False,default=3,
		help='Number of technical replicates')
	parser.add_argument('--nbrep','-b',required=False,default=3,
		help='Number of biological replicates')
	parser.add_argument('--start_well','-s',required=False,default='A1',
		help='Well in upper left corner of starting region')
	required.add_argument('--controls','-c',required=False,
		help=''''
			Comma-separated list of controls including actual sample.
			Example: 'gfp_full,gfp_noRT,rRNA_full,rRNA_noRT'
			'''
			)
	parsed=parser.parse_args()
	input_df = read_input(parsed.input)
	map = create_map(parsed.xbyy,int(parsed.ntrep),int(parsed.nbrep),parsed.vertical,
		parsed.random,input_df,[x for x in parsed.controls.split(',')],
		parsed.start_well)
	name = parsed.input.split('/')[-1].split('.')[0]
	map.to_csv('{0}_map.tsv'.format(name),sep='\t')
	
	
		
		
		
				
				
	
	
	