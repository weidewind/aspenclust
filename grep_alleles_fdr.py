## collects global results
import argparse
import xlsxwriter
import pandas as pd
import numpy as np
import os
import datetime
import re

## take options: input folder, output tag
parser = argparse.ArgumentParser(description='Collect the results of global analysis into xlsx')
parser.add_argument("--input", dest='analysis_folder', type=str, help="full path to the folder containing analysis results")
parser.add_argument("--fdr", dest='fdr_folder', type=str, help="full path to the folder containing fdr results")
parser.add_argument("--state", dest='state', type=str, default="nsyn")
args = parser.parse_args()

proteins = ["h1", "n1", "n2", "h3", "h1pandemic", "n1pandemic"]

output = os.path.join(args.analysis_folder,'Alleles_with_FDR'+'.xlsx')
writer = pd.ExcelWriter(output, engine='xlsxwriter',datetime_format='YYYY-MM-DD HH:MM:SS',options={'strings_to_numbers': True}) 
# Add a format. Green fill with dark green text.
workbook = writer.book
green = workbook.add_format({'bg_color': '#C6EFCE',
							   'font_color': '#006100'})
wrapit = workbook.add_format({'text_wrap':True,'align': 'left', 'border':1,'bold':True})


for prot in proteins:
	fdrfile = os.path.join(args.fdr_folder, prot+'_'+args.state+'_mean_alleles_FDR')
	sitefile = os.path.join(args.analysis_folder,'Alleles.xlsx')
	if not os.path.isfile(fdrfile) or  not os.path.isfile(sitefile):
		continue
	fdrdf = pd.read_csv(fdrfile,  sep='	', names= ["site", "ancder", "par", "div", "diff", "pvalue"])
	sitesdf = pd.read_excel(sitefile, sheet_name=prot)
	if (fdrdf.shape[0] == 0 or sitesdf.shape[0] == 0):
		continue
	fdr = []
	for pval in sitesdf["pvalue_mean"]:
		f=(fdrdf[fdrdf["pvalue"]<=pval].shape[0]*sitesdf.shape[0])/(fdrdf.shape[0]*sitesdf[sitesdf["pvalue_mean"]<=pval].shape[0])
		fdr.append(f)

	sitesdf["FDR"] = fdr
	if 'pvalue_mean' in sitesdf.columns:
		pdf = sitesdf[(sitesdf.protein == prot)].sort_values(by = "pvalue_mean")
	elif 'pvalue_median' in sitesdf.columns:
                pdf = sitesdf[(sitesdf.protein == prot)].sort_values(by = "pvalue_median")
	pdf.to_excel(writer, sheet_name=prot,index = False, startrow=1,header = False) # Convert the dataframe to an XlsxWriter Excel object. Skip one row to print header later (otherwise user-defined format will be ignored)
	worksheet = writer.sheets[prot]
	range1 = 'G1:G'+str(len(sitesdf.index)+1)
	range2 = 'K1:K'+str(len(sitesdf.index)+1)
	worksheet.conditional_format(range1, {'type':     'cell',
                                        'criteria': '<=',
                                        'value':    0.01,
					'multi_range': range1+' '+range2,
                                        'format':   green})
	for colnum, value in enumerate(pdf.columns.values):
		worksheet.write(0,colnum,value,wrapit)
pd.Series({'input folder':args.analysis_folder, 'fdr input folder':args.fdr_folder, 'time':datetime.datetime.now()}).to_excel(writer, sheet_name="readme", index = False)
writer.save()
	


