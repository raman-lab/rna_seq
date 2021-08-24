#Script to compare 7b fold induction data to novaseq/pacbio data
import argparse
import numpy as np
import pandas as pd
import scipy.optimize
import matplotlib
matplotlib.use('agg')
from matplotlib import pyplot as plt
from collections import OrderedDict
import re
import traceback

def read_input(input_file):
	variants = pd.read_csv(input_file,dtype='str',index_col=None)
	out_df = pd.read_csv(input_file,index_col=None)
	out_df['variant']=variants['variant']
	return out_df

def re_index(df1,df2):
	df2.set_index(['variant'],inplace=True)
	indices = df2.index.values.tolist()
	df1.set_index(['variant'],inplace=True)
	df1 = df1[df1.index.isin(df2.index)]
	new_df1 = df1.reindex(indices)
	df2.reset_index(inplace=True)
	new_df1.reset_index(inplace=True)
	return new_df1,df2

def fit_line(x,m,b):
	return (m*x)+b

def rmscalc(yact,yfit):
	return np.sqrt(np.sum((yfit-yact)**2)/yfit.shape[0])

def rsquared_calc(yact,yfit):
	residuals = yact - yfit
	ss_res = np.sum(residuals**2)
	ss_tot = np.sum((yact-np.mean(yact))**2)
	return 1-(ss_res/ss_tot)
	
def fit_curve(df1,df2,df1_name,df2_name,log):
	fit_xvals = np.arange(0,max(df1.iloc[:,1])*1.25,0.05)
	fig,ax = plt.subplots(1,1)
	fig.set_size_inches(8,6)
	if log == True:
		df2.iloc[:,1]=np.log(df2.iloc[:,1])
		df2.iloc[:,2]=np.log(df2.iloc[:,2])
	try:
		popt,pcov = scipy.optimize.curve_fit(fit_line,df1.iloc[:,1],df2.iloc[:,1],
			method='lm')
	except:
		traceback.print_exc()
		quit()
	else:
		temp_list = [x for x in popt]
		for num in np.sqrt(np.diag(pcov)):
			temp_list.append(num)
		rmserror = rmscalc(df2.iloc[:,1],fit_line(df1.iloc[:,1],*popt))
		residuals = df2.iloc[:,1]-fit_line(df1.iloc[:,1],*popt)
		ax.errorbar(fit_xvals,fit_line(fit_xvals,*popt),color='red',label='fit: m = %5.3f\n b = %5.3f\n' % tuple(popt))
		temp_list.append(rmserror)
		temp_list.append(rsquared_calc(df2.iloc[:,1],fit_line(df1.iloc[:,1],*popt)))
		ax.errorbar(df1.iloc[:,1],df2.iloc[:,1],fmt='o',
			color='blue',capsize=5, elinewidth=1, markeredgewidth=1,
			label='data\n {0}'.format(rmserror)
			)
		ax.set_xlabel(df1_name)
		ax.set_ylabel(df2_name)
#		ax.set_xscale('log')
		plt.legend(loc=2)
		plt.savefig('{0}_{1}_fit.pdf'.format(df1_name,df2_name),transparent=True)
		plt.close()
		fig,ax = plt.subplots(1,1)
		ax.scatter(df1.iloc[:,1],residuals)
		ax.set_xlabel(df1_name)
		ax.set_ylabel('residuals')
		plt.savefig('{0}_{1}_residuals.pdf'.format(df1_name,df2_name),transparent=True)
		plt.close()
		out_df = pd.DataFrame([temp_list])
		out_df.columns = ['m','b','m_err','b_err','rmserror','rsquared']
		return out_df

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Script to plot two input files and calculate the best fit line between them')
	parser.add_argument('--input','-i',required=True,nargs=2,
		help='2 input files separated by a space')
	parser.add_argument('--log','-l',required=False,default=False,action='store_true',
		help='Log of the second input file')
	parsed = parser.parse_args()
	in_dict = OrderedDict()
	for file in parsed.input:
		base = file.split('/')[-1].split('.')[0]
		in_dict[base]=read_input(file)
	df2,df1 = re_index(in_dict[in_dict.keys()[1]],in_dict[in_dict.keys()[0]])
	print df1
	print df2
	out_df = fit_curve(df1,df2,in_dict.keys()[0],in_dict.keys()[1],parsed.log)
	out_df.to_csv('{0}_{1}_fit_metrics.csv'.format(in_dict.keys()[0],in_dict.keys()[1]),
		index=False)