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
parser.add_argument("--filter", dest='filter', type=str)
#parser.add_argument("--protein", dest='prot', type=str)
args = parser.parse_args()

proteins = ["h1", "n1", "n2", "h3", "h1pandemic", "n1pandemic"]
tofilter = []
if (args.filter):
	tofilter = args.filter.split(",")

output = os.path.join(args.analysis_folder,'Groups_with_FDR'+'.xlsx')
writer = pd.ExcelWriter(output, engine='xlsxwriter',datetime_format='YYYY-MM-DD HH:MM:SS',options={'strings_to_numbers': True}) 
# Add a format. Green fill with dark green text.
workbook = writer.book
green = workbook.add_format({'bg_color': '#C6EFCE',
							   'font_color': '#006100'})
wrapit = workbook.add_format({'text_wrap':True,'align': 'left', 'border':1,'bold':True})


for prot in proteins:
	fdrfile = os.path.join(args.fdr_folder, prot+'_nsyn_mean_group_FDR')
	if ( not os.path.exists(fdrfile)):
		continue
	fdrdf = pd.read_csv(fdrfile,  sep='	', names= ["fake", "group_name", "group_pvalue", "enrich_pvalue"])
	fdrdf = fdrdf[~fdrdf["group_name"].isin(tofilter)]
	groupsdf = pd.read_excel(os.path.join(args.analysis_folder,'Groups.xlsx'),  sheet_name=prot)
	groupsdf = groupsdf[~groupsdf["group"].isin(tofilter)]
	if (fdrdf.shape[0] == 0 or groupsdf.shape[0] == 0):
		continue
		
	groupfdr = []
	print (groupsdf["group_pvalue"])
	for pval in groupsdf["group_pvalue"]:
		f=(fdrdf[fdrdf["group_pvalue"]<=pval].shape[0]*groupsdf.shape[0])/(fdrdf.shape[0]*groupsdf[groupsdf["group_pvalue"]<=pval].shape[0])
		groupfdr.append(f)
	print (groupfdr)
	groupsdf["group_FDR"] = groupfdr
	
	enrichfdr = []
	print (groupsdf["enrich_pvalue"])
	for pval in groupsdf["enrich_pvalue"]:
		f=(fdrdf[fdrdf["enrich_pvalue"]<=pval].shape[0]*groupsdf.shape[0])/(fdrdf.shape[0]*groupsdf[groupsdf["enrich_pvalue"]<=pval].shape[0])
		enrichfdr.append(f)
	print (groupfdr)
	groupsdf["enrich_FDR"] = enrichfdr
	
	pdf = groupsdf[(groupsdf.protein == prot)].sort_values(by = "group_pvalue")
	
	pdf.to_excel(writer, sheet_name=prot,index = False, startrow=1,header = False) # Convert the dataframe to an XlsxWriter Excel object. Skip one row to print header later (otherwise user-defined format will be ignored)
	worksheet = writer.sheets[prot]
	range1 = 'K1:K'+str(len(groupsdf.index)+1)
	range2 = 'L1:L'+str(len(groupsdf.index)+1)
	worksheet.conditional_format(range1, {'type':     'cell',
                                        'criteria': '<=',
                                        'value':    0.05,
					'multi_range': range1+' '+range2,
                                        'format':   green})
	for colnum, value in enumerate(pdf.columns.values):
		worksheet.write(0,colnum,value,wrapit)
pd.Series({'input folder':args.analysis_folder, 'fdr input folder':args.fdr_folder, 'time':datetime.datetime.now()}).to_excel(writer, sheet_name="readme", index = False)
writer.save()
	


