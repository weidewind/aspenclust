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
parser.add_argument("--type", dest='type', default='alleles', type=str, help="alleles or ancestors")
parser.add_argument("--input", dest='folder', type=str, help="full path to the folder containing results")
args = parser.parse_args()

## list through the directory and find all files matching the pattern $protein_$state_global_$statype_statistics
prots = []

take_nonempty = lambda s1, s2: s1 if s2.isna().all() else s2

def get_files(folder, type):
	files = []
	for r, d, f in os.walk(folder):
		for file in f:
			if "_" + type + "_" in file:
				files.append(file)
	files.sort()
	proteins = np.unique([f.split("_")[0] for f in files])
	print("Got some files for ")
	print(proteins)
	return files, proteins

files, proteins = get_files(args.folder, args.type)

## foreach file - read the results, parse them and write to a dataframe
def parse_file(folder, filename,prots):
	with open(os.path.join(folder, filename)) as f:
		spl = filename.split("_") # $protein_$state_global_$statype_statistics
		subtype = "ancder" if args.type == "alleles" else "anc"
		newdf = pd.DataFrame(columns=['protein','site', subtype, 'parallel_'+spl[3], 'divergent_'+spl[3],'difference_'+spl[3], 'pvalue_'+spl[3]])
		lines = f.readlines()
		prots.append(spl[0])
		num = 0
		outstr = ""
		while num < len(lines):
			if re.match(r'^>site', lines[num]):
				res = lines[num+1].strip()
				values = res.split("\t")
				row = pd.Series({'protein':spl[0],'site':values[0], subtype:values[1], 'parallel_'+spl[3]:values[2], 'divergent_'+spl[3]:values[3], 'difference_'+spl[3]:values[4], 'pvalue_'+spl[3]:values[5]})
				newdf.loc[len(newdf.index)] = row
			num += 1
		return (newdf)

meandf = pd.DataFrame(columns=['protein','site'])
mediandf = pd.DataFrame(columns=['protein','site'])

for filename in files:
	df = parse_file(args.folder, filename,prots)
	if 'pvalue_mean' in df.columns:
		meandf = pd.concat([meandf, df], sort=False)
	elif 'pvalue_median' in df.columns:
                mediandf = pd.concat([mediandf, df], sort=False)
print(meandf)
print(mediandf)
if not meandf.empty and not mediandf.empty:
	df = meandf.merge(mediandf, on = ['protein', 'site', 'subs'])
else:
	if not meandf.empty:
		df = meandf
	else:
		df = mediandf
# Create a Pandas Excel writer using XlsxWriter as the engine.
output = os.path.join(args.folder, args.type.capitalize() + '.xlsx')
writer = pd.ExcelWriter(output, engine='xlsxwriter',datetime_format='YYYY-MM-DD HH:MM:SS',options={'strings_to_numbers': True}) 
# Add a format. Green fill with dark green text.
workbook = writer.book
green = workbook.add_format({'bg_color': '#C6EFCE',
                               'font_color': '#006100'})
wrapit = workbook.add_format({'text_wrap':True,'align': 'left', 'border':1,'bold':True})

range1 = 'G1:G'+str(len(df.index)+1)
range2 = 'K1:K'+str(len(df.index)+1)
for prot in proteins:
	pdf = df[(df.protein == prot)].sort_values(by = "pvalue_mean")
	pdf.to_excel(writer, sheet_name=prot,index = False, startrow=1,header = False) # Convert the dataframe to an XlsxWriter Excel object. Skip one row to print header later (otherwise user-defined format will be ignored)
	worksheet = writer.sheets[prot]
	worksheet.conditional_format(range1, {'type':     'cell',
                                        'criteria': '<=',
                                        'value':    0.01,
					'multi_range': range1+' '+range2,
                                        'format':   green})
	for colnum, value in enumerate(pdf.columns.values):
		worksheet.write(0,colnum,value,wrapit)
pd.Series({'input folder':args.folder, 'time':datetime.datetime.now()}).to_excel(writer, sheet_name="readme", index = False)
writer.save()
