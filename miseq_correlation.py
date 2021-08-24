#Python script to map replicate NGS runs against each other.

import argparse
import pandas as pd
import numpy as np
import scipy.optimize
import matplotlib
matplotlib.use('agg')
from matplotlib import pyplot as plt
from collections import OrderedDict
import re
import traceback
import itertools

def read_input(input_file):
	variants = pd.read_csv(input_file,dtype='str',index_col=None)
	out_df = pd.read_csv(input_file,index_col=None)
	out_df['barcode']=variants['barcode']
	return out_df

def totals_cutoffs(in_df,cutoff):
	temp_df = in_df.copy()
	temp_df.set_index('barcode',inplace=True)
	temp_df = temp_df[temp_df>cutoff]
	temp_df.dropna(how='any',axis=0,inplace=True)
	perc_df = temp_df/temp_df.sum()
#	print perc_df
	perc_df.dropna(how='any',axis=0,inplace=True)
#	print perc_df
	perc_df.reset_index()
	temp_df.reset_index()
	return temp_df, perc_df

def get_combinations(in_df):
	temp = [x for x in itertools.combinations(in_df.columns.values.tolist()[1:],2)]
	return temp

def fit_line(x,m,b):
	return (m*x)+b

def rmscalc(yact,yfit):
	return np.sqrt(np.sum((yfit-yact)**2)/yfit.shape[0])

def rsquared_calc(yact,yfit):
	residuals = yact - yfit
	ss_res = np.sum(residuals**2)
	ss_tot = np.sum((yact-np.mean(yact))**2)
	return 1-(ss_res/ss_tot)

def fit_curve(df,df1_name,df2_name):
	fit_xvals = np.arange(0,max(df.iloc[:,1])*1.25,0.05)
	fig,ax = plt.subplots(1,1)
	fig.set_size_inches(8,6)
	try:
		popt,pcov = scipy.optimize.curve_fit(fit_line,df[df1_name],df[df2_name],
			method='lm')
	except:
		traceback.print_exc()
		quit()
	else:
		temp_list = [x for x in popt]
		for num in np.sqrt(np.diag(pcov)):
			temp_list.append(num)
		rmserror = rmscalc(df[df2_name],fit_line(df[df1_name],*popt))
		ax.errorbar(fit_xvals,fit_line(fit_xvals,*popt),color='red',label='fit: m = %5.3f\n b = %5.3f\n' % tuple(popt))
		temp_list.append(rmserror)
		temp_list.append(rsquared_calc(df[df2_name],fit_line(df[df1_name],*popt)))
		ax.errorbar(df[df1_name],df[df2_name],fmt='o',
			color='blue',capsize=5, elinewidth=1, markeredgewidth=1,
			label='data\n {0}'.format(rmserror),alpha=0.25
			)
		ax.set_xlabel(df1_name)
		ax.set_ylabel(df2_name)
#		ax.set_xscale('log')
		plt.legend(loc=2)
		plt.savefig('{0}_{1}_fit.pdf'.format(df1_name,df2_name),transparent=True)
		plt.close()
		out_df = pd.DataFrame([temp_list])
		out_df.columns = ['m','b','m_err','b_err','rmserror','rsquared']
		return out_df

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Script to plot two input files and calculate the best fit line between them')
	parser.add_argument('--input','-i',required=True,
		help='Input file containing barcode counts for multiple replicates')
	parser.add_argument('--cutoff','-c',required=False,default=1,
		help='Cutoff file for minimum sequence counts required for analysis')
	parser.add_argument('--percentage','-p',required=False,default=False,
		action='store_true',help='Calculate correlation of percentages instead of raw counts')
	parsed = parser.parse_args()
	in_df = read_input(parsed.input)
	cutoff_df,perc_df = totals_cutoffs(in_df,int(parsed.cutoff))
	combos = get_combinations(in_df)
	for x in combos:
		if parsed.percentage == True:
			out_df = fit_curve(perc_df,x[0],x[1])
		else:
			out_df = fit_curve(cutoff_df,x[0],x[1])
		out_df.to_csv('{0}_{1}_fit_metrics.csv'.format(
			x[0],x[1]),
			index=False)