__author__ = 'davidc'
import csv, re

f = csv.reader(open('../promoter.txt', 'rb'), delimiter='\t')
t = 0
pattern = re.compile("[A-Z]")

for row in f:
    if row[0][0] == '#': #skip comments inside the txt file
        continue
    #print row
    promName = row[1];
    strand = row[2];
    TSS = row[3];
    seq = row[9];
    if seq == '' or TSS == '' :
        continue #This row is actually empty, why does regulonDB include it if we know so little
    m = pattern.search(seq)
    if m == None :
        print '%s has no TSS!' % seq
    t = t + 1__author__ = 'davidc'
