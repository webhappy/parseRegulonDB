__author__ = 'davidc'
import csv
from Bio import SeqIO
from Bio.Seq import Seq
import cPickle

refFasta = SeqIO.read("sequence.fasta", 'fasta')
#Dict mapping column header (must match exactly) to expt name (arbitrary)
exptCols={'t8_2_LR':'Aerobic #2','t8_3_LR':'Aerobic #3'}
correspondingCols={} #Store column header:Column #
exptResults={}
SEQ_COL = 2  # Hard-code for now that the sgRNA sequence is in column 3
            # Ideally, there would be a column with the position and I wouldn't need this column
sgRNAseqs=[]
positions=[]
strand=[]

with open('../All data - 32992 with LRs.csv','rb') as csvfile:
    reader=csv.reader(csvfile)
    header=reader.next()

    for k in range(len(header)):
        for j in exptCols.iterkeys():
            if str(header[k]) == j:
                #This column matches exptCols{j}
                correspondingCols[j]=k

    assert(len(correspondingCols.items())==len(exptCols.items()))
    for k in exptCols.iterkeys():
        exptResults[k]=[];  # Initialize to an empty list

    for row in reader:
        sgRNAseqs.append(row[SEQ_COL])
        for k in exptCols.iterkeys():
            exptResults[k].append(row[correspondingCols[k]])


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
            positions.append(l+1)  #Notice I increment by 1 since strings are 0-indexed while NT coordinates are 1-indexed

cPickle.dump((exptResults,positions,strand),open('expt_results.txt', 'wb'))