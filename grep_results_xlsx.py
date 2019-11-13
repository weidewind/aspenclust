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
df = pd.DataFrame(columns=['protein','state','stattype','parallel','divergent','difference','pvalue'])
def get_files(folder):
	files = []
	for r, d, f in os.walk(folder):
		for file in f:
			if "_global_" in file:
				files.append(file)
	files.sort()
	proteins = np.unique([f.split("_")[0] for f in files])
	print("Got some files for ")
	print(proteins)
	return files, proteins

files, proteins = get_files(args.folder)

## foreach file - read the results, parse them and write to a dataframe
def parse_file(folder, filename,df):
	with open(os.path.join(folder, filename)) as f:
		spl = filename.split("_") # $protein_$state_global_$statype_statistics
		res = f.readlines()[-1].strip()
		if not res:
			res = f.readlines()[-2].strip()
		values = res.split("\t")
		print(values)
		if len(values) != 5:
			print("Something is wrong with " + filename + " : its last line is " + res + ", expected 5 tab-delimited strings")
			return
		print(spl)
		print(values)
		df.loc[len(df.index)] = pd.Series({'protein':spl[0],'state':spl[1], 'stattype':spl[3], 'parallel':values[1], 'divergent':values[2], 'difference':values[3], 'pvalue':values[4]})

for filename in files:
	parse_file(args.folder, filename, df)

## Try to find some results for the other state
other_state_folder = re.sub("/nsyn/","/syn/",args.folder)
if args.folder == other_state_folder:
	other_state_folder = re.sub("/syn/","/nsyn/",args.folder)
print("Also searching in "+other_state_folder)
if os.path.exists(other_state_folder):
	files, proteins = get_files(other_state_folder)
	for filename in files:
	        parse_file(other_state_folder, filename, df)

# Create a Pandas Excel writer using XlsxWriter as the engine.
output = os.path.join(args.folder,'Global.xlsx')
shname = args.folder.split("/")[-1]
writer = pd.ExcelWriter(output, engine='xlsxwriter',datetime_format='YYYY-MM-DD HH:MM:SS',options={'strings_to_numbers': True}) 
df.to_excel(writer, sheet_name=shname, index = False) # Convert the dataframe to an XlsxWriter Excel object.
pd.Series({'input folder':args.folder, 'time':datetime.datetime.now()}).to_excel(writer, sheet_name="readme", index = False)

# Add a format. Green fill with dark green text.
workbook = writer.book
green = workbook.add_format({'bg_color': '#C6EFCE',
                               'font_color': '#006100'})
range = 'G2:G'+str(len(df.index)+1)
worksheet = writer.sheets[shname]
worksheet.conditional_format(range, {'type':     'cell',
                                        'criteria': '<=',
                                        'value':    0.05,
                                        'format':   green})
writer.save()

## write readme: input folder and date of collection
#workbook = xlsxwriter.Workbook(output,{'strings_to_numbers': True})
#worksheet = workbook.add_worksheet("readme")
#worksheet.write(0,1,args.folder)
#worksheet.write(1,1,datetime.datetime.now())
#workbook.close()
