import cPickle
import pandas.io.sql as psql
import pandas.io.parsers
import mysql.connector

results=cPickle.load(open('pickled.txt'))
designed=pandas.io.parsers.read_csv(open('sgRNAs_to_design.csv','rb'),index_col='seq')
print len(results)

numPerfect=0
for (seq,records) in results.iteritems():
    seq=seq.lower()
    allPerfect=True
    expData.loc[seq,'numOff']=len(records)
    expData.loc[seq,'hasOffTargets']=True
    numPerfect = 0
    for record in records:
        if record.identities < 20:
            numPerfect += 1
            allPerfect=False
    if allPerfect:
        expData.loc[seq,'onlyPerfectOffTargets']=True
        numPerfect+=1
        print expData.loc[seq]

print numPerfect,'have only matches of 20 bases'

expData.to_csv('offTargetResults.csv')