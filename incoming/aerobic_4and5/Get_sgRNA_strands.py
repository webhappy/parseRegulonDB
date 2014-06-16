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
count=0
for s in sgRNAseqs:
    count+=1
    if count % 10 == 0:
        print count
    s=s.upper()
    #look in reference for this sequence
    l=refSeq.find(s)
    if l > -1:
        strand.append('+')
        positions.append(l+1)  #Notice I increment by 1 since strings are 0-indexed while NT coordinates are 1-indexed
    else:
        t=Seq(s)
        strand.append('-')
        l=refSeq.find(str(t.reverse_complement()))
        if l == -1:
            unmapped = True
        positions.append(l+1)  #Notice I increment by 1 since strings are 0-indexed while NT coordinates are 1-indexed

print '\n%i sgRNAseqs' % len(sgRNAseqs)
print '%i positions mapped' % len(positions)
print '%i strands found' % len(strand)
if unmapped:
    print '\nThere were some unmapped sgRNAs'
if ~unmapped:
    print '\nAll sgRNAs were perfectly mapped to at least one chromosomal location'

# Write the results to file

ofile = open(r'/Users/Toro/Documents/Dropbox/Data/Cas9/Miseq/Analysis/sgRNA strand finding/sgRNA positions.csv', 'w')
ofile.write('seq,sgRNA_pos,sgRNA_strand\n')
for i in range(len(sgRNAseqs)):
    ofile.write('%s,%i,%s\n' %(sgRNAseqs[i], positions[i],strand[i]))
