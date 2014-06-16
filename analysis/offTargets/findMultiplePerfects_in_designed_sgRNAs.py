import cPickle
import pandas.io.sql as psql
import pandas.io.parsers
import mysql.connector
import csv

results=cPickle.load(open('pickled.txt'))
designed=csv.reader(open('sgRNAs_to_design.csv','rU'))
print len(results)

sgRNAs_with_perfect_mismatches = []
numPerfect=0
for (seq,records) in results.iteritems():
    seq=seq.lower()
    numPerfect = 0
    for record in records:
        if record.identities == 20:
            numPerfect += 1
            if numPerfect == 2:
                sgRNAs_with_perfect_mismatches.append(seq)

print len(sgRNAs_with_perfect_mismatches),'have at least 2 20-bp matches'

designed.next()
for (name,dist,seq,logFC,counts,number,primer,plate1,plate2,primer) in designed:
    #print seq
    if seq in sgRNAs_with_perfect_mismatches:
        print seq,'has multiple perfect matches!'