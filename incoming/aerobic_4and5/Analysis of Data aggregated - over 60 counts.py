# -*- coding: utf-8 -*-
"""
Created on Tue Jan  7 20:06:23 2014

@author: Toro
"""
import csv
import numpy as np
import os
#Set working Directory
os.chdir(r'/Users/Toro/Documents/Dropbox/Data/Cas9/Miseq/Analysis/')

# First, make a file with all the data from sRNAs ###

# Import file with aggregated data
def import_csv_into_dict(filename):
    """Imports a .csv file with the variables on the first row into a dictionary.
    The dictionary is in the form Data[row_number][variable] = value.
    The working directory must point to where the file is.
    The number of cases (rows) must be given as an argument"""
    print 'Importing %s' % filename
    ifile = open(filename)
    reader = csv.reader(ifile)
    num_rows = 0
    for row in reader:
        num_rows += 1 # count the number of rows

    firstline = True
    Data = {}
    variables = []
    ifile.seek(0)
    for i in range(1, (num_rows)): # initialize a dictionary for each row. 
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
    return Data

# Import data
over60 = import_csv_into_dict(r'Data aggregated - over 60 counts.csv')  #8612
sRNAs = import_csv_into_dict(r'OLS Master key - sRNA oligos - 344.csv') #344

# Add data from over60 into sRNAs
for query in range(1, len(sRNAs)+1):
    for rec in range(1, len(over60)+1):
        if sRNAs['%s' % query]['RegulonDBID'] == over60['%s' % rec]['RegulonDBID']:
            for var in over60['%s' % rec].keys():
                if var in sRNAs['%s' % query].keys():
                    continue
                else:
                    sRNAs['%s' % query][var] = over60['%s' % rec][var]
                    

# Save sRNAs as file

#Writing function
def write_csv(filename, data_dict):
    """Writes a dictionary to filename.csv.  
    Dictionary must be in the form: filename[record][variable] = value
    Variable names are gotten from data_dict.keys()[0] and put into a dictionary entry. 
    Then the script goes through filename.keys() line by line and adds them to the file.
    (NOTE: all commas in the values are deleted to avoid registering problems)
    """
    print 'Writing %s.csv file' % (filename)
    ofile = open('%s.csv' %(filename), 'w')
    ofile.write('%s\n' % ','.join(data_dict[data_dict.keys()[0]].keys()))  # Write the Headers
    for rec in data_dict.keys():
        line = []
        for variable in data_dict[rec].keys():
            line.append(str(data_dict[rec][variable]).replace(',',''))# Create a list with all the data, forcing them as strings, and stripping off any commas
        ofile.write('%s\n' % ','.join(line))
    ofile.close()
    return

write_csv('sRNAs with LRs', sRNAs)
