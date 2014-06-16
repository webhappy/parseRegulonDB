import csv
from Bio import SeqIO
from Bio.Seq import Seq
import cPickle

# Get E. coli gemone
refFasta = SeqIO.read(r'/Users/Toro/Documents/Dropbox/Data/Cas9/Miseq/Analysis/sgRNA strand finding/NC_000913_2.fas', 'fasta')

# Get sgRNA sequences
Data_file = r'/Users/Toro/Documents/Dropbox/Data/Cas9/Miseq/Analysis/OLS Master key - unique oligos 32992.csv'
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

# Put sgRNA seqs in a single list
sgRNAseqs=[]
for i in Data.keys():
    sgRNAseqs.append(Data[i]['seq'])

positions=[]
strand=[]

# Here I begin to calculate the strand from the sgRNAseqs
refSeq=str(refFasta.seq) # search in here
refSeqRC=str(refFasta.seq.reverse_complement()) # search in here if can't find in forward strand
monogamous_sgRNAs = 0
unmapped_sgRNAs = 0
degenerate_sgRNAs = 0
unmap_sgRNAs = []
deg_sgRNAs = []
count = 0
for s in sgRNAseqs:
    count+=1
    if count % 100 == 0:
        print count
    s=s.upper()
    #look in reference for this sequence
    l = refSeq.count(s)
    m = refSeqRC.count(s)
    hits = l + m
    if hits == 0:
        unmap_sgRNAs.append(s)
        unmapped_sgRNAs += 1
    elif hits == 1:
        monogamous_sgRNAs += 1
    elif hits > 1:
        degenerate_sgRNAs += 1
        deg_sgRNAs.append(s)

print '%i monogamous sgRNAs' % monogamous_sgRNAs
print '%i unmapped sgRNAs' % unmapped_sgRNAs
print '%i degenerate sgRNAs' % degenerate_sgRNAs
print '%i total sgRNAs checked, out of 32992 total' % (monogamous_sgRNAs+unmapped_sgRNAs+degenerate_sgRNAs)

# Write the degenerate sgRNAs to file
ofile = open(r'/Users/Toro/Documents/Dropbox/Data/Cas9/Miseq/Analysis/sgRNA strand finding/degenerate sgRNAs.csv', 'w')
ofile.write('seq\n')
for i in range(len(deg_sgRNAs)):
    ofile.write('%s\n' % deg_sgRNAs[i])
