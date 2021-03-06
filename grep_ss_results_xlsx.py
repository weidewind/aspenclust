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
parser.add_argument("--input", dest='folder', type=str, help="full path to the folder containing results")
args = parser.parse_args()

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
			if "_sites_" in file:
				files.append(file)
	files.sort()
	proteins = np.unique([f.split("_")[0] for f in files])
	print("Got some files for ")
	print(proteins)
	return files, proteins

files, proteins = get_files(args.folder)

## foreach file - read the results, parse them and write to a dataframe
def parse_file(folder, filename,prots):
	with open(os.path.join(folder, filename)) as f:
		spl = filename.split("_") # $protein_$state_global_$statype_statistics
		newdf = pd.DataFrame(columns=['protein','site', 'subs', 'parallel_'+spl[3], 'divergent_'+spl[3],'difference_'+spl[3], 'pvalue_'+spl[3]])
		lines = f.readlines()
		prots.append(spl[0])
		num = 0
		outstr = ""
		while num < len(lines):
			if re.match(r'^key', lines[num]):
				subs = {}
				num += 1
				while not re.match(r'^>', lines[num]) and not re.match(r'^Hist', lines[num]):
					splitter = lines[num].split()
					subs[re.findall(r'[a-zA-Z]+', splitter[0])[0]] = splitter[1]
					num += 1
				sorted_keys = sorted(subs.keys())
				outstr = ""
				sep = ","
				for i in range(len(sorted_keys)):
					if i+1 == len(sorted_keys):
						sep = ""
					elif sorted_keys[i][:3] != sorted_keys[i+1][:3]:
						sep = ";"
					else:
						sep = ","
					outstr = outstr+sorted_keys[i]+subs[sorted_keys[i]]+sep
			if re.match(r'^>site', lines[num]):
				res = lines[num+1].strip()
				values = res.split("\t")
				row = pd.Series({'protein':spl[0],'site':values[0], 'subs':outstr, 'parallel_'+spl[3]:values[1], 'divergent_'+spl[3]:values[2], 'difference_'+spl[3]:values[3], 'pvalue_'+spl[3]:values[4]})
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
output = os.path.join(args.folder,'Sites'+'.xlsx')
writer = pd.ExcelWriter(output, engine='xlsxwriter',datetime_format='YYYY-MM-DD HH:MM:SS',options={'strings_to_numbers': True}) 
# Add a format. Green fill with dark green text.
workbook = writer.book
green = workbook.add_format({'bg_color': '#C6EFCE',
                               'font_color': '#006100'})
wrapit = workbook.add_format({'text_wrap':True,'align': 'left', 'border':1,'bold':True})

range1 = 'G1:G'+str(len(df.index)+1)
range2 = 'K1:K'+str(len(df.index)+1)
for prot in proteins:
	if 'pvalue_mean' in df.columns:
		pdf = df[(df.protein == prot)].sort_values(by = "pvalue_mean")
	elif 'pvalue_median' in df.columns:
                pdf = df[(df.protein == prot)].sort_values(by = "pvalue_median")
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
