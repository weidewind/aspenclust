## collects global results
import argparse
import xlsxwriter
import pandas as pd
import numpy as np
import os
import datetime
import re
import subprocess

## take options: input folder, output tag
parser = argparse.ArgumentParser(description='Collect the results of global analysis into xlsx')
parser.add_argument("--input", dest='folder', type=str, help="full path to the folder containing results")
args = parser.parse_args()

subprocess.run(["perl", "grep_groups_results.pl", "--input", args.folder])

## list through the directory and find all files matching the pattern $protein_$state_global_$statype_statistics
#df = pd.DataFrame(columns=['site','stattype','parallel','divergent','difference','pvalue','subs'])
#cols = ['protein','site','parallel_mean', 'divergent_mean', 'difference_mean', 'pvalue_mean','parallel_median', 'divergent_median', 'difference_median', 'pvalue_median']
#df = pd.DataFrame(columns=cols)
prots = []



take_nonempty = lambda s1, s2: s1 if s2.isna().all() else s2

def get_files(folder):
	files = []
	for r, d, f in os.walk(folder):
		for file in f:
			if "_groups_results_" in file:
				files.append(file)
	files.sort()
	proteins = np.unique([f.split("_")[0] for f in files])
	print("Got some files for ")
	print(proteins)
	return files, proteins

files, proteins = get_files(args.folder)
sdf = pd.DataFrame()
ldf = pd.DataFrame()
for filename in files:
	newdf = pd.read_csv(os.path.join(args.folder, filename),index_col=False)
	newdf = newdf.round(7)
	print("looking for "+newdf['protein'][0])
#	if df.empty:
#		df = newdf
#	elif newdf['protein'][0] in set(df['protein']):
#		print("Merging..")
#		df = df.merge(newdf, how='outer',on = ['protein', 'group', 'stat_type', 'group_count', 'group_same', 'group_diff', 'compl_count', 'compl_same', 'compl_diff','diffdiff'])
#	else:
#		print("Concatting..")
	if "_site" in filename:
		if sdf.empty:
			sdf = newdf
		else:
			sdf = pd.concat([sdf, newdf], sort=False)
	elif "_label" in filename:
		if ldf.empty:
			ldf = newdf
		else:
			ldf = pd.concat([ldf, newdf], sort=False)

df = ldf.merge(sdf, how='outer',on = ['protein', 'group', 'stat_type', 'group_count', 'group_same', 'group_diff', 'compl_count', 'compl_same', 'compl_diff','diffdiff'])

# Create a Pandas Excel writer using XlsxWriter as the engine.
output = os.path.join(args.folder,'Groups'+'.xlsx')
writer = pd.ExcelWriter(output, engine='xlsxwriter',datetime_format='YYYY-MM-DD HH:MM:SS',options={'strings_to_numbers': True}) 
# Add a format. Green fill with dark green text.
workbook = writer.book
green = workbook.add_format({'bg_color': '#C6EFCE',
                               'font_color': '#006100'})
wrapit = workbook.add_format({'text_wrap':True,'align': 'left', 'border':1,'bold':True})

range1 = 'K1:K'+str(len(df.index)+1)
range2 = 'L1:L'+str(len(df.index)+1)
range3 = 'N1:N'+str(len(df.index)+1)
range4 = 'R1:R'+str(len(df.index)+1)

for prot in proteins:
	pdf = df[(df.protein == prot)].sort_values(by = "group")
	pdf.to_excel(writer, sheet_name=prot,index = False, startrow=1,header = False) # Convert the dataframe to an XlsxWriter Excel object. Skip one row to print header later (otherwise user-defined format will be ignored)
	worksheet = writer.sheets[prot]
	worksheet.conditional_format(range1, {'type':     'cell',
                                        'criteria': '<=',
                                        'value':    0.05,
					'multi_range': range1+' '+range2+' '+range3+' '+range4,
                                        'format':   green})
	for colnum, value in enumerate(pdf.columns.values):
		worksheet.write(0,colnum,value,wrapit)
pd.Series({'input folder':args.folder, 'time':datetime.datetime.now()}).to_excel(writer, sheet_name="readme", index = False)
writer.save()
