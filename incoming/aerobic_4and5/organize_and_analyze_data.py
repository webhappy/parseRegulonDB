#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 16 21:03:28 2013

@author: Toro

Script to aggregate the data from my sequencing, starting from the mapping_summary1.csv file
"""

    ### Open each of the files and make a dictionary that contains all the data ###

import csv
import numpy as np
import os

#Set working Directory
os.chdir(r'/Users/Toro/Documents/Dropbox/Data/Cas9/Miseq/')

mapping_files = [r'2013-07-19 OLS in LB_1/mapping_summary_1.csv', 
r'2013-11-03 OLS in LB_2/mapping_summary_1.csv', 
r'2013-11-21 OLS in LB3 - Timecourse/mapping_summary_1.csv', 
r'2013-12-26 OLS no O2/mapping_summary_1.csv', 
r'2013-12-28 OLS no O2 2-3/mapping_summary_1.csv',
r'Analysis/OLS Master key - unique oligos 32992.csv']

#Writing function
def write_csv(filename, data_dict, variables):
    """Writes a dictionary to filename.csv.
    Dictionary must be in the form: filename[record][variable] = value
    Variable names are gotten from 'variables' and put into a dictionary entry.
    Then the script goes through filename.keys() line by line and adds them to the file.
    (NOTE: all commas in the values are deleted to avoid registering problems)
    """
    print 'Writing %s.csv file' % (filename)
    ofile = open('Analysis/%s.csv' %(filename), 'w')
    ofile.write('%s\n' % ','.join(variables))  # Write the Headers
    for rec in data_dict.keys():
        line = []
        for variable in variables:
            line.append(str(data_dict[rec][variable]).replace(',',''))  # Create a list with all the data, forcing them as strings, and stripping off any commas
        ofile.write('%s\n' % ','.join(line))
    ofile.close()
    return


    ### First, get all the data into a dictionary ###
    ### Dictionary looks like this: Data[record(unique to each sgRNA, goes from 1-32993)][variable] = 'value' ###

vars_done = 0 # keep track of how many variables have been entered, to avoid overwriting
firstfile = True 
Data = {} # Initialize the ultimate Data dictionary
print '\n'
Data['variables'] = []
for map_file in mapping_files:
    #Open file
    print 'Getting data from %s' % map_file  
    ifile = open('%s' % (map_file))
    reader = csv.reader(ifile)
    firstline = True
    row_counter = 1 
    for row in reader:
        var_counter = 0
        if firstline:    # get variable names
            Data['variables'].extend(row)
            firstline = False
            continue
        else:           # get Data values
            for var in Data['variables']: # iterate over variables
                if var_counter < vars_done: # if there are variables already in the file,
                    var_counter += 1        # skip them
                    continue
                else:
                    if firstfile:   # if this is the first file
                        for i in range(1, 32993): # initialize a dictionary for each row
                            Data['%s' % (i)] = {}
                        firstfile = False
                    Data['%s' % row_counter][var] = row[var_counter-vars_done]  # get data
                    var_counter += 1
        row_counter += 1
    vars_done = var_counter          
    ifile.close()

Vars = Data['variables']
del Data['variables']

# Calculate the average t4_1 and t8_1
for rec in Data.keys():
    Data[rec]['count_t4_1'] = sum([float(Data[rec]['count_t4_1_1']), float(Data[rec]['count_t4_1_2'])])
    Data[rec]['count_t8_1'] = sum([float(Data[rec]['count_t8_1_1']), float(Data[rec]['count_t8_1_2'])])
Vars.extend(['count_t4_1', 'count_t8_1'])

    ### testing module ###
Data_variables = ['bnum', 'product_name', 'gene_name', 'CategoryID', 'sgRNA_pos', 'sgRNA_strand', 'seq', 'RegulonDBID', 'PAMID', 'ref',
                  'count_O2_t0_1', 'count_O2_d24_1', 'count_O2_t0_2', 'count_O2_d24_2', 'count_O2_t0_3', 'count_O2_d24_3',
                  'count_t0_1', 'count_t4_1_minus', 'count_t4_1_1', 'count_t4_1_2', 'count_t4_1',
                  'count_t8_1_1', 'count_t8_1_2', 'count_t8_1', 'count_t0_2', 'count_t8_2',
                  'count_t0_3', 'count_t4_3', 'count_t8_3', 'count_t13_3']
write_csv('All data - 32992', Data, Data_variables)
    ###  ###

    ### Now aggregate the data so that there is only one record per feature ###

# First make a list of regDBIDs, with a list of recs for each
print '\nMaking a dictionary of regDBIDs, with their respective locations'
regDBIDs = {}
counter = 1
for rec in Data.keys():
    if Data[rec]['RegulonDBID'] not in regDBIDs.keys():
        regDBIDs['%s' % Data[rec]['RegulonDBID']] = [rec]
        counter +=1
        if counter%1000 == 0:
            print 'Calculated %i ECKs out of ~10500' % counter
    else:
        regDBIDs['%s' % Data[rec]['RegulonDBID']].append(rec)

# Now make the aggregated Data
print '\nAggregating the data'
Da = {}
for ECK in regDBIDs.keys():  # for each regulonDBID
    Da[ECK] = {}
    for var in Vars: # for every variable
        if 'count' in var:
            temp = 0
            for rec in regDBIDs[ECK]:  # for each of the sgRNAs that hit that ECK
                temp += int(Data[rec][var])  # sum the counts
            Da[ECK][var] = temp  # add them to Data aggregated (Da)
        elif 'freq' in var:
            continue  # Ignore the freqs
        else:
            Da[ECK][var] = Data[regDBIDs[ECK][0]][var]  # keep the first record
            
    ### Get rid of things that have too few counts and calculate the log ratios ###
print '\nCalculating log ratios and getting rid of low-count records'
Da_under60 = {}  # initialize the low count data dict
for rec in Da.keys():
    # calculate the mean log ratios, add one extra count to avoid dividing by zero
    Da[rec]['t8_1_LR'] = np.log2(Da[rec]['count_t8_1']+1) - np.log2(Da[rec]['count_t0_1']+1)
    Da[rec]['t8_2_LR'] = np.log2(Da[rec]['count_t8_2']+1) - np.log2(Da[rec]['count_t0_2']+1)
    Da[rec]['t8_3_LR'] = np.log2(Da[rec]['count_t8_3']+1) - np.log2(Da[rec]['count_t0_3']+1)
    Da[rec]['O2_1_LR'] = np.log2(Da[rec]['count_O2_d24_1']+1) - np.log2(Da[rec]['count_O2_t0_1']+1)
    Da[rec]['O2_2_LR'] = np.log2(Da[rec]['count_O2_d24_2']+1) - np.log2(Da[rec]['count_O2_t0_2']+1)
    Da[rec]['O2_3_LR'] = np.log2(Da[rec]['count_O2_d24_3']+1) - np.log2(Da[rec]['count_O2_t0_3']+1)
    # Pick out the records that have < 60 counts and throw them out
    Da[rec]['total_t0_count'] = sum([int(Da[rec]['count_t0_1']), int(Da[rec]['count_t0_2']), int(Da[rec]['count_t0_3'])])
    if Da[rec]['total_t0_count'] < 60:
        Da_under60[rec] = Da[rec]  # copy that entry to Da_under60
        del Da[rec]  # delete that entry

        
#Calculate the median of t8_LR and O2_LR and subtract that from each entry to normalize
print '\nNormalizing the Log ratios'
LRs_t8_1 = []
LRs_t8_2 = []
LRs_t8_3 = []
LRs_O2_1 = []
LRs_O2_2 = []
LRs_O2_3 = []

for rec in Da.keys():
    LRs_t8_1.append(Da[rec]['t8_1_LR'])
    LRs_t8_2.append(Da[rec]['t8_2_LR'])
    LRs_t8_3.append(Da[rec]['t8_3_LR'])
    LRs_O2_1.append(Da[rec]['O2_1_LR'])
    LRs_O2_2.append(Da[rec]['O2_2_LR'])
    LRs_O2_3.append(Da[rec]['O2_3_LR'])

for rec in Da.keys():
    Da[rec]['t8_1_LR'] = Da[rec]['t8_1_LR'] - np.median(LRs_t8_1)
    Da[rec]['t8_2_LR'] = Da[rec]['t8_2_LR'] - np.median(LRs_t8_2)
    Da[rec]['t8_3_LR'] = Da[rec]['t8_3_LR'] - np.median(LRs_t8_3)
    Da[rec]['O2_1_LR'] = Da[rec]['O2_1_LR'] - np.median(LRs_O2_1)
    Da[rec]['O2_2_LR'] = Da[rec]['O2_2_LR'] - np.median(LRs_O2_2)
    Da[rec]['O2_3_LR'] = Da[rec]['O2_3_LR'] - np.median(LRs_O2_3)

print 'There were %i ECKs with fewer than 60 counts' % len(Da_under60.keys())
print 'There were %i ECKs with more than 60 counts\n' % len(Da.keys())

        
        
    
    ### Output .csv files with the data ###

    ### Set up the variables that I want to write to each file ###

master_variables = ['bnum', 'product_name', 'gene_name', 'CategoryID', 'sgRNA_pos',
                  'sgRNA_strand', 'seq', 'RegulonDBID', 'PAMID', 'ref']
count_variables = ['count_O2_t0_1', 'count_O2_d24_1', 'count_O2_t0_2', 'count_O2_d24_2', 'count_O2_t0_3',
                   'count_O2_d24_3', 'count_t0_1', 'count_t4_1_minus', 'count_t4_1_1', 'count_t4_1_2', 'count_t4_1',
                  'count_t8_1_1', 'count_t8_1_2', 'count_t8_1', 'count_t0_2', 'count_t8_2', 'count_t0_3',
                  'count_t4_3', 'count_t8_3', 'count_t13_3']
LR_variables = ['t8_1_LR', 't8_2_LR', 't8_3_LR', 'O2_1_LR', 'O2_2_LR', 'O2_3_LR']

write_csv('All data - 32992', Data, master_variables + count_variables)  # Save the whole dataset
write_csv('Data aggregated - over 60 counts', Da, LR_variables + master_variables + count_variables)  # Save the aggregated data
write_csv('Data aggregated - under 60 counts', Da_under60, LR_variables + master_variables + count_variables)  # Save the records with <60 counts
write_csv('For Honglei/Data aggregated - over 60 counts', Da, LR_variables + master_variables + count_variables)
write_csv('For Honglei/Data aggregated - under 60 counts', Da_under60, LR_variables + master_variables + count_variables)

 




































