# -*- coding: utf-8 -*-
"""
Created on Thu Dec 26 15:04:46 2013

@author: Toro
"""

import csv
import numpy as np
import matplotlib.pyplot as plt
import os

    ### Import whole dataset ###

filename = r'/Users/Toro/Documents/Dropbox/Data/Cas9/Miseq/Analysis/All data - 32992 with LRs.csv'

# Set up the Dictionary that the data will go into
Data = {}
for i in range(1, 32993): # initialize a dictionary for each row
    Data['%s' % (i)] = {}
Data['variables'] = [] # initialize a dictionary for the variable names

print '\n'
print 'Getting data from %s' % filename

#Open file
ifile = open('%s' % (filename))
reader = csv.reader(ifile)
firstline = True
row_counter = 1
for row in reader:
    if firstline:    # get variable names
        Data['variables'].extend(row)
        firstline = False
        continue
    else:           # get Data values
        var_counter = 0
        for var in Data['variables']: # iterate over variables
            Data['%s' % (row_counter)][var] = row[var_counter]  # get data
            var_counter += 1
    row_counter += 1
ifile.close()

del Data['variables']

    ### Make the reproducibility plots

# Set up lists with the different variables
print 'Making lists of LRs'
aerobic_t0_3_LR = []
aerobic_t4_3_LR = []
aerobic_t8_3_LR = []
aerobic_t13_3_LR = []
anaerobic_d24_1_LR = []
anaerobic_d24_2_LR = []
anaerobic_d24_3_LR = []
aerobic_t8_1_LR = []
aerobic_t8_2_LR = []


for rec in Data.iterkeys():
    aerobic_t0_3_LR.append(Data[rec]['aerobic_t0_3_LR'])
    aerobic_t4_3_LR.append(Data[rec]['aerobic_t4_3_LR'])
    aerobic_t8_3_LR.append(Data[rec]['aerobic_t8_3_LR'])
    aerobic_t13_3_LR.append(Data[rec]['aerobic_t13_3_LR'])
    anaerobic_d24_1_LR.append(Data[rec]['anaerobic_d24_1_LR'])
    anaerobic_d24_2_LR.append(Data[rec]['anaerobic_d24_2_LR'])
    anaerobic_d24_3_LR.append(Data[rec]['anaerobic_d24_3_LR'])
    aerobic_t8_1_LR.append(Data[rec]['aerobic_t8_1_LR'])
    aerobic_t8_2_LR.append(Data[rec]['aerobic_t8_2_LR'])

print 'Done making lists'


# Scatterplots of sgRNA counts between technical and biological replicates

print 'Calculating plots'
plt.figure()
plt.subplot(231, share)
plt.scatter(aerobic_t8_1_LR, aerobic_t8_2_LR, s=0.1)
plt.ylabel('2')
plt.xlabel('1')
plt.subplot(232)
plt.scatter(aerobic_t8_1_LR, aerobic_t8_3_LR, s=0.1)
plt.ylabel('3')
plt.xlabel('1')
plt.subplot(233)
plt.scatter(aerobic_t8_2_LR, aerobic_t8_3_LR, s=0.1)
plt.ylabel('3')
plt.xlabel('2')
plt.subplot(234)
plt.scatter(anaerobic_d24_1_LR, anaerobic_d24_2_LR, s=0.1)
plt.ylabel('2')
plt.xlabel('1')
plt.subplot(235)
plt.scatter(anaerobic_d24_1_LR, anaerobic_d24_3_LR, s=0.1)
plt.ylabel('3')
plt.xlabel('1')
plt.subplot(236)
plt.scatter(anaerobic_d24_2_LR, anaerobic_d24_3_LR, s=0.1)
plt.ylabel('3')
plt.xlabel('2')
fig = plt.gcf()
fig.set_size_inches(8,8)
plt.savefig(r'/Users/Toro/Desktop/testplot.png',dpi=100)

# Consistency between sgRNAs on the same feature line graphs

# Consistency between sgRNAs SD histogram

# Plot the technical and biological reproducibility data

# Technical reproducibility

# Biological reproducibility

# Plot the sgRNA consistency data

# Histogram of SDs of sgRNAs

# Plots of individual sgRNAs


# Save the plots to a file