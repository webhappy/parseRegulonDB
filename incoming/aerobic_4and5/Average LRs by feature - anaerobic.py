# -*- coding: utf-8 -*-
"""
Created on Mon Jan  6 17:06:22 2014

@author: Toro
"""

    ### This script calculates the log ratios for all 32992 sgRNAs
    ### and averages them by ECK and replicate (only one value is kept for
    ### each ECK

import csv
import numpy as np
import os

#Set working Directory
os.chdir(r'/Users/Toro/Documents/Dropbox/Data/Cas9/Miseq/Analysis/For David')

# import the file with all the data
print 'Importing data'
Data_file = r'Anaerobic.csv'
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
counter = 0
for ECK in regDBIDs.keys():  # for each regulonDBID
    counter += 1
    if counter % 500 == 0:
        print 'Done with %i of 10560 ECKs' % counter
    Da[ECK] = {}
    temp = []
    for rec in Data.keys():
        first = True
        if Data[rec]['RegulonDBID'] == ECK:
            if first:
                Da[ECK]['sgRNA_strand'] = Data[rec]['sgRNA_strand']
                Da[ECK]['RegulonDBID'] = ECK
                Da[ECK]['sgRNA_pos'] = Data[rec]['sgRNA_pos']
                first = False
            temp.extend([float(Data[rec]['anaerobic_d24_1_LR']),
                        float(Data[rec]['anaerobic_d24_2_LR']),
                        float(Data[rec]['anaerobic_d24_3_LR'])])
    Da[ECK]['LR_avg'] = np.average(temp)


#Calculate the median of t8_LR and O2_LR and subtract that from each entry to normalize
print '\nNormalizing the Log ratios'
LRs = []

for rec in Da.keys():
    LRs.append(Da[rec]['LR_avg'])

LRs_median = np.median(LRs)

for rec in Da.keys():
    Da[rec]['LR_avg'] = Da[rec]['LR_avg'] - LRs_median

print 'The median LR was %f' % np.median(LRs)
print 'There were %i ECKs' % len(Da)
print 'There were %i LRs averaged' % len(LRs)


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


# Output .csv files with the data ###
write_csv(r'For David/Anaerobic_averages', Da, ['sgRNA_pos', 'RegulonDBID', 'sgRNA_strand', 'LR_avg'])



















