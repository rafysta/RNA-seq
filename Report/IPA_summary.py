# -*- coding: utf-8 -*-
import argparse
import xlsxwriter
import codecs
import re
import numpy as np
import pandas as pd
from pandas import Series, DataFrame

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--ipa", help="IPA output file")
parser.add_argument("-o", "--out", help="output file")
args = parser.parse_args()


FILE_ipa = args.ipa
FILE_out = args.out


#==========================================
# load IPA
#==========================================
DATA_pathway = DataFrame()
DATA_regulator = DataFrame()
DATA_function = DataFrame()
DATA_network = DataFrame()
DATA_tox = DataFrame()

Flag_pathway = 0; NUM_pathway = 0
Flag_regulator = 0; NUM_regulator = 0
Flag_function = 0; NUM_function =0
Flag_network = 0; NUM_network = 0
Flag_tox = 0; NUM_tox = 0

fh_in = open(FILE_ipa, 'rt')
line = fh_in.readline()
while line:
    # pathway
    if Flag_pathway == 1 and line.rstrip('\r\n') == "":
        Flag_pathway = 0
    
    if Flag_pathway == 1:
        DATA_pathway.loc[NUM_pathway] =(map(str.strip, line.rstrip('\r\n').split('\t')))[:-1]
        NUM_pathway += 1
        
    if Flag_pathway == 2:
        DATA_pathway = DataFrame(columns=map(str.strip, line.rstrip('\r\n').split('\t')))
        Flag_pathway = 1
    
    if re.match("Canonical Pathways", line):
        Flag_pathway = 2


    # regulators
    if Flag_regulator == 1 and line.rstrip('\r\n') == "":
        Flag_regulator = 0
    
    if Flag_regulator == 1:
        DATA_regulator.loc[NUM_regulator] = map(str.strip, line.rstrip('\r\n').split('\t'))
        NUM_regulator += 1

    if Flag_regulator == 2:
        DATA_regulator = DataFrame(columns=map(str.strip, line.rstrip('\r\n').split('\t')))
        Flag_regulator = 1
    
    if re.match("Upstream Regulators", line):
        Flag_regulator = 2
    
    
    # disease and function
    if Flag_function == 1 and line.rstrip('\r\n') == "":
        Flag_function = 0
    
    if Flag_function == 1:
        DATA_function.loc[NUM_function] = map(str.strip, line.rstrip('\r\n').split('\t'))
        NUM_function += 1
        
    if Flag_function == 2:
        DATA_function = DataFrame(columns=map(str.strip, line.rstrip('\r\n').split('\t')))
        Flag_function = 1
    
    if re.match("Diseases and Bio Functions", line):
        Flag_function = 2
        
    
    # network
    if Flag_network == 1 and line.rstrip('\r\n') == "":
        Flag_network = 0
    
    if Flag_network == 1:
        DATA_network.loc[NUM_network] = map(str.strip, line.rstrip('\r\n').split('\t'))
        NUM_network += 1
        
    if Flag_network == 2:
        DATA_network = DataFrame(columns=map(str.strip, line.rstrip('\r\n').split('\t')))
        Flag_network = 1
    
    if re.match("Networks for My Projects", line):
        Flag_network = 2
    
    
    # tox
    if Flag_tox == 1 and line.rstrip('\r\n') == "":
        Flag_tox = 0
    
    if Flag_tox == 1:
        DATA_tox.loc[NUM_tox] =map(str.strip, line.rstrip('\r\n').split('\t'))
        NUM_tox += 1
        
    if Flag_tox == 2:
        DATA_tox = DataFrame(columns=map(str.strip, line.rstrip('\r\n').split('\t')))
        Flag_tox = 1
    
    if re.match("Tox Lists for", line):
        Flag_tox = 2
    
    line = fh_in.readline()
fh_in.close

workbook = xlsxwriter.Workbook(FILE_out)
worksheet_pathway = workbook.add_worksheet('Pathways')
worksheet_regulator = workbook.add_worksheet('Regulators')
worksheet_function = workbook.add_worksheet('Functions')
worksheet_network = workbook.add_worksheet('Networks')
worksheet_tox = workbook.add_worksheet('Tox')

format_sig = workbook.add_format({'bg_color': '#D1E68F', 'font_color': '#000000'})
format_blue = workbook.add_format({'bg_color': '#BDD7EE', 'font_color': '#0070C0'})
format_red = workbook.add_format({'bg_color': '#E8A4A4', 'font_color': '#C00000'})
format_bold = workbook.add_format({'bold': True})
format_normal = workbook.add_format({'bold': False})


#==========================================
# pathway
#==========================================
DATA_pathway = DATA_pathway.iloc[::-1]
DATA_pathway = DATA_pathway.reset_index(drop=True)
worksheet_pathway.set_column('A:A', 53)
worksheet_pathway.set_column('B:E', 8)
header = ["Canonical Pathways", "p-value", "Ratio", "number", "Molecules"]
string_col = [0,4]
string_name = ["Ingenuity Canonical Pathways", "Molecules"]
num_col = [2]
num_name = ["Ratio"]
worksheet_pathway.write_row(0, 0, header, format_bold)
for i in range(0,len(DATA_pathway)):
    row=i+1
    for col, name in zip(string_col, string_name):
        worksheet_pathway.write_string(row, col, unicode(DATA_pathway.loc[i, name], 'utf-8'))
    for col, name in zip(num_col, num_name):
        if DATA_pathway.loc[i, name] != u'NaN' and DATA_pathway.loc[i, name] != "":
            worksheet_pathway.write_number(row, col, float(DATA_pathway.loc[i, name]))
    
    pvalue = 10**(float(DATA_pathway.loc[i, "-log(p-value)"]) * -1)
    worksheet_pathway.write_number(row, 1, pvalue)
    
    molecule_count = len(DATA_pathway.loc[i, "Molecules"].split(","))
    worksheet_pathway.write_number(row, 3, int(molecule_count))

Range = "B2:B" + str(len(DATA_pathway)+1)
worksheet_pathway.conditional_format(Range, {'type': 'cell',
                                         'criteria': '<',
                                         'value': 0.05,
                                         'format': format_sig})


#==========================================
# Regulators
#==========================================
worksheet_regulator.set_column('A:A', 25)
worksheet_regulator.set_column('B:B', 20)
worksheet_regulator.set_column('C:E', 8)
header = ["Type", "Regulator", "p-value", "number", "Target molecules"]
string_col = [0,1,4]
string_name = ["Molecule Type", "Upstream Regulator", "Target molecules in dataset"]
num_col = [2]
num_name = ["p-value of overlap"]
worksheet_regulator.write_row(0, 0, header, format_bold)
for i in range(0,len(DATA_regulator)):
    row=i+1
    for col, name in zip(string_col, string_name):
        worksheet_regulator.write_string(row, col, unicode(DATA_regulator.loc[i, name], 'utf-8'))
    for col, name in zip(num_col, num_name):
        if DATA_regulator.loc[i, name] != u'NaN' and DATA_regulator.loc[i, name] != "":
            worksheet_regulator.write_number(row, col, float(DATA_regulator.loc[i, name]))
    
    molecule_count = len(DATA_regulator.loc[i, "Target molecules in dataset"].split(","))
    worksheet_regulator.write_number(row, 3, int(molecule_count))

Range = "C2:C" + str(len(DATA_regulator)+1)
worksheet_regulator.conditional_format(Range, {'type': 'cell',
                                         'criteria': '<',
                                         'value': 0.05,
                                         'format': format_sig})
                                         
                                         
                                         
#==========================================
# Functions
#==========================================
worksheet_function.set_column('A:A', 32)
worksheet_function.set_column('B:C', 44)
worksheet_function.set_column('D:F', 8)
header = ["Categories", "Function", "Diseases or Functions Annotation", "p-value","number", "Molecules"]
string_col = [0,1,2,5]
string_name = ["Categories", "Functions", "Diseases or Functions Annotation", "Molecules"]
num_col = [3,4]
num_name = ["p-Value", "# Molecules"]
worksheet_function.write_row(0, 0, header, format_bold)
for i in range(0,len(DATA_function)):
    row=i+1
    for col, name in zip(string_col, string_name):
        worksheet_function.write_string(row, col, unicode(DATA_function.loc[i, name], 'utf-8'))
    for col, name in zip(num_col, num_name):
        if DATA_function.loc[i, name] != u'NaN' and DATA_function.loc[i, name] != "":
            worksheet_function.write_number(row, col, float(DATA_function.loc[i, name]))

Range = "D2:D" + str(len(DATA_function)+1)
worksheet_function.conditional_format(Range, {'type': 'cell',
                                         'criteria': '<',
                                         'value': 0.05,
                                         'format': format_sig})


#==========================================
# Networks
#==========================================
worksheet_network.set_column('A:A', 5)
worksheet_network.set_column('B:B', 67)
worksheet_network.set_column('C:R', 8)
header = ["Rank", "Functions", "Score", "number", "Molecules"]
string_col = [1,4]
string_name = ["Top Diseases and Functions", "Molecules in Network"]
num_col = [0,2]
num_name = ["ID", "Score"]
worksheet_network.write_row(0, 0, header, format_bold)
for i in range(0,len(DATA_network)):
    row=i+1
    for col, name in zip(string_col, string_name):
        worksheet_network.write_string(row, col, unicode(DATA_network.loc[i, name], 'utf-8'))
    for col, name in zip(num_col, num_name):
        if DATA_network.loc[i, name] != u'NaN' and DATA_network.loc[i, name] != "":
            worksheet_network.write_number(row, col, float(DATA_network.loc[i, name]))

    molecule_count = len(DATA_network.loc[i, "Molecules in Network"].split(","))
    worksheet_network.write_number(row, 3, int(molecule_count))


#==========================================
# Tox
#==========================================
worksheet_tox.set_column('A:A', 53)
worksheet_tox.set_column('B:E', 8)
header = ["Toxicity list", "p-value", "Ratio", "number", "Molecules"]
string_col = [0,4]
string_name = ["Ingenuity Toxicity Lists", "Molecules"]
num_col = [2]
num_name = ["Ratio"]
worksheet_tox.write_row(0, 0, header, format_bold)
for i in range(0,len(DATA_tox)):
    row=i+1
    for col, name in zip(string_col, string_name):
        worksheet_tox.write_string(row, col, unicode(DATA_tox.loc[i, name], 'utf-8'))
    for col, name in zip(num_col, num_name):
        if DATA_tox.loc[i, name] != u'NaN' and DATA_tox.loc[i, name] != "":
            worksheet_tox.write_number(row, col, float(DATA_tox.loc[i, name]))
    
    pvalue = 10**(float(DATA_tox.loc[i, "-log(p-value)"]) * -1)
    worksheet_tox.write_number(row, 1, pvalue)
    
    molecule_count = len(DATA_tox.loc[i, "Molecules"].split(","))
    worksheet_tox.write_number(row, 3, int(molecule_count))

Range = "B2:B" + str(len(DATA_tox)+1)
worksheet_tox.conditional_format(Range, {'type': 'cell',
                                         'criteria': '<',
                                         'value': 0.05,
                                         'format': format_sig})
                                         
                                         
workbook.close()
