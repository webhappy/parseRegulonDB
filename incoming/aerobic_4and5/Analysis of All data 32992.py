# -*- coding: utf-8 -*-
"""
Created on Mon Jan  6 17:06:22 2014

@author: Toro
"""

    ### This script calculates the log ratios for all 32992 sgRNAs ###

import csv
import numpy as np
import os

#Set working Directory
os.chdir(r'/Users/Toro/Documents/Dropbox/Data/Cas9/Miseq/Analysis/')

# import the file with all the data
print 'Importing data'
Data_file = r'All data - 32992.csv'
ifile = open(Data_file)
reader = csv.reader(ifile)
firstline = True
Data = {}
variables = []
for i in range(1, 32993): # initialize a dictionary for each row
    Data['%s' % i] = {}

row_counter = 1
for row in reader:
    if firstline:
        variables.extend(row)
        firstline = False
    else:   
        var_counter = 0
        for var in variables:
            Data['%s' % row_counter][var] = row[var_counter]  # get data
            var_counter += 1 
        row_counter += 1
ifile.close()

  
# Calculate the log ratios
print 'Calculating Log Ratios'

for rec in Data.keys():
    Data[rec]['aerobic_t0_3_LR'] = 1
    Data[rec]['aerobic_t4_3_LR'] = np.log2(int(Data[rec]['count_t4_3'])+1) - np.log2(int(Data[rec]['count_t0_3'])+1)
    Data[rec]['aerobic_t8_3_LR'] = np.log2(int(Data[rec]['count_t8_3'])+1) - np.log2(int(Data[rec]['count_t0_3'])+1)
    Data[rec]['aerobic_t13_3_LR'] = np.log2(int(Data[rec]['count_t13_3'])+1) - np.log2(int(Data[rec]['count_t0_3'])+1)        
    Data[rec]['anaerobic_d24_1_LR'] = np.log2(int(Data[rec]['count_O2_d24_1'])+1) - np.log2(int(Data[rec]['count_O2_t0_1'])+1)
    Data[rec]['anaerobic_d24_2_LR'] = np.log2(int(Data[rec]['count_O2_d24_2'])+1) - np.log2(int(Data[rec]['count_O2_t0_2'])+1)
    Data[rec]['anaerobic_d24_3_LR'] = np.log2(int(Data[rec]['count_O2_d24_3'])+1) - np.log2(int(Data[rec]['count_O2_t0_3'])+1)
    Data[rec]['aerobic_t8_1_LR'] = np.log2(int(float(Data[rec]['count_t8_1']))+1) - np.log2(int(float(Data[rec]['count_t0_1']))+1) # the extra float is require because count_t8_1 is an average, so is float)
    Data[rec]['aerobic_t8_2_LR'] = np.log2(int(Data[rec]['count_t8_2'])+1) - np.log2(int(Data[rec]['count_t0_2'])+1)
    
 
#Calculate the median of LogRatios and subtract that from each entry to normalize
print 'Normalizing Log Ratios'
LRs_t4_3 = []
LRs_t8_3 = []
LRs_t13_3 = []
LRs_O2_1 = []
LRs_O2_2 = []
LRs_O2_3 = []
LRs_t8_1 = []
LRs_t8_2 = []

for rec in Data.keys():
    LRs_t4_3.append(Data[rec]['aerobic_t4_3_LR'])
    LRs_t8_3.append(Data[rec]['aerobic_t8_3_LR'])
    LRs_t13_3.append(Data[rec]['aerobic_t13_3_LR'])
    LRs_O2_1.append(Data[rec]['anaerobic_d24_1_LR'])
    LRs_O2_2.append(Data[rec]['anaerobic_d24_2_LR'])
    LRs_O2_3.append(Data[rec]['anaerobic_d24_3_LR'])
    LRs_t8_1.append(Data[rec]['aerobic_t8_1_LR'])
    LRs_t8_2.append(Data[rec]['aerobic_t8_2_LR'])

print 'Updating records with normalized values'
mLRs_t4_3 = np.median(LRs_t4_3)
mLRs_t8_3 = np.median(LRs_t8_3)
mLRs_t13_3 = np.median(LRs_t13_3)
mLRs_O2_1 = np.median(LRs_O2_1)
mLRs_O2_2 = np.median(LRs_O2_2)
mLRs_O2_3 = np.median(LRs_O2_3)
mLRs_t8_1 = np.median(LRs_t8_1)
mLRs_t8_2 = np.median(LRs_t8_2)

for rec in Data.keys():
    Data[rec]['aerobic_t4_3_LR'] = Data[rec]['aerobic_t4_3_LR'] - mLRs_t4_3
    Data[rec]['aerobic_t8_3_LR'] = Data[rec]['aerobic_t8_3_LR'] - mLRs_t8_3
    Data[rec]['aerobic_t13_3_LR'] = Data[rec]['aerobic_t13_3_LR'] - mLRs_t13_3
    Data[rec]['anaerobic_d24_1_LR'] = Data[rec]['anaerobic_d24_1_LR'] - mLRs_O2_1
    Data[rec]['anaerobic_d24_2_LR'] = Data[rec]['anaerobic_d24_2_LR'] - mLRs_O2_2
    Data[rec]['anaerobic_d24_3_LR'] = Data[rec]['anaerobic_d24_3_LR'] - mLRs_O2_3
    Data[rec]['aerobic_t8_1_LR'] = Data[rec]['aerobic_t8_1_LR'] - mLRs_t8_1
    Data[rec]['aerobic_t8_2_LR'] = Data[rec]['aerobic_t8_2_LR'] - mLRs_t8_2

#Writing function
def write_csv(filename, data_dict, variables):
    """Writes a dictionary to filename.csv.
    Dictionary must be in the form: filename[record][variable] = value
    Variable names are gotten from 'variables' and put into a dictionary entry.
    Then the script goes through filename.keys() line by line and adds them to the file.
    (NOTE: all commas in the values are deleted to avoid registering problems)
    """
    print 'Writing %s.csv file' % (filename)
    ofile = open('/Users/Toro/Documents/Dropbox/Data/Cas9/Miseq/Analysis/%s.csv' %(filename), 'w')
    ofile.write('%s\n' % ','.join(variables))  # Write the Headers
    for rec in data_dict.keys():
        line = []
        for variable in variables:
            line.append(str(data_dict[rec][variable]).replace(',',''))  # Create a list with all the data, forcing them as strings, and stripping off any commas
        ofile.write('%s\n' % ','.join(line))
    ofile.close()
    return


    ### Set up the variables that I want to write to each file ###

master_variables = ['bnum', 'product_name', 'gene_name', 'CategoryID', 'sgRNA_pos',
                  'sgRNA_strand', 'seq', 'RegulonDBID', 'PAMID', 'ref']
count_variables = ['count_O2_t0_1', 'count_O2_d24_1', 'count_O2_t0_2', 'count_O2_d24_2', 'count_O2_t0_3',
                   'count_O2_d24_3', 'count_t0_1', 'count_t4_1_minus', 'count_t4_1_1', 'count_t4_1_2', 'count_t4_1',
                  'count_t8_1_1', 'count_t8_1_2', 'count_t8_1', 'count_t0_2', 'count_t8_2', 'count_t0_3',
                  'count_t4_3', 'count_t8_3', 'count_t13_3']
LR_variables = ['aerobic_t0_3_LR', 'aerobic_t4_3_LR', 'aerobic_t8_3_LR', 'aerobic_t13_3_LR',
                'anaerobic_d24_1_LR', 'anaerobic_d24_2_LR', 'anaerobic_d24_3_LR',
                'aerobic_t8_1_LR', 'aerobic_t8_2_LR']


write_csv('All data - 32992 with LRs', Data, LR_variables + master_variables + count_variables)
write_csv('For David/Aerobic', Data, ['RegulonDBID', 'sgRNA_pos', 'seq', 'sgRNA_strand', 'aerobic_t8_1_LR', 'aerobic_t8_2_LR', 'aerobic_t8_3_LR'])
write_csv('For David/Anaerobic', Data, ['RegulonDBID', 'sgRNA_pos', 'seq', 'sgRNA_strand', 'anaerobic_d24_1_LR', 'anaerobic_d24_2_LR', 'anaerobic_d24_3_LR'])





















