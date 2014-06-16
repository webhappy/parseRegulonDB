import pandas.io.sql as psql
import numpy as np
import csv, mysql.connector

con=mysql.connector.connect(user='root', host='127.0.0.1', database='t')
genes=psql.read_frame('select gene_name,gene_strand,gene_posleft,gene_posright,operon_pos,median as avg_lr from GENE',con)
print "Big join SQL query done!"
PECgenes=psql.read_frame('select gene_name from PEC',con)
PECgenes=PECgenes['gene_name'].values

tails=[]  # Append pointers to tails here
operonLengths = {}  # Dict mapping name of tail gene
for strand in ['forward','reverse']:
    last=None
    genesInDir = genes[(genes['gene_strand']==strand)]
    if strand == 'forward':
        genesInDir = genesInDir.sort('gene_posleft',ascending=True)
    else:
        genesInDir = genesInDir.sort('gene_posright',ascending=False)
    for (gene,strand,left,right,pos,avgLR) in genesInDir.itertuples(index=False):
        if gene =='' or gene=='None':
            continue
        if last is None:
            count=1
            last={'gene':gene,'val':avgLR,'prev':None}
        elif pos==1:  # we're starting a new operon here
            #print "%s has count %i"%(gene,count)
            tails.append(last)
            operonLengths[last['gene']] = count
            last={'gene':gene,'val':avgLR,'prev':None}
            count=1
        else:
            count+=1
            cur={'gene':gene,'prev':last,'val':avgLR}
            last=cur
    operonLengths[last['gene']] = count
    tails.append(last)

pointersToTails=tails
cutoff = -3
knownEssentials = PECgenes

outF = csv.writer(open('oddCases.csv','wb'))
outF.writerow(['Gene','logFC','Operon length',"This gene's position"])

determinedSick=[]
determinedNonSick=[]
undetermined=[]
unexpectedHealthy = []  # Genes that are healthy despite someone downstream being sick
for pointer in pointersToTails:
    cur=pointer
    lastGene = pointer['gene']
    canDetermine=True
    curPosition = 0
    while cur is not None:
        curPosition += 1
        if cur['val'] > cutoff:
            determinedNonSick.append(cur['gene'])
            if not canDetermine:
                outF.writerow([cur['gene'],cur['val'],operonLengths[lastGene],curPosition])
                unexpectedHealthy.append((cur['gene'],cur['val']))
        elif cur['val'] <= cutoff:
            if not canDetermine:
                undetermined.append(cur['gene'])
            else:
                determinedSick.append(cur['gene'])
                canDetermine=False
        else:
            print cur['gene'],'has no sgRNAs'
        cur=cur['prev']

TP = 0;  # True positives
missingEssentials = 0;  # False negative
undeterminedEssential = 0;

for gene in knownEssentials:
    if gene in determinedSick:
        TP+=1
    elif gene in determinedNonSick:
        missingEssentials+=1;
    elif gene in undetermined:
        undeterminedEssential+=1
    else:
        print gene,"not found in any of the arrays!"

# Want the true negative count
TN = 0
for gene in determinedNonSick:
    if gene not in knownEssentials:
        TN += 1

# How many actual negatives are there? Total # of genes - # in PEC
knownNumNonEssential = len(determinedNonSick) + len(determinedSick) + len(undetermined) - len(knownEssentials)

print "TP=%i, missing essentials=%i, undetermined essentials=%i, size of presumed essentials=%i, size of undetermined=%i"%\
      (TP,missingEssentials,undeterminedEssential,len(determinedSick),len(undetermined))
print len(unexpectedHealthy)

#for (gene,val) in unexpectedHealthy:
#    outF.writerow([gene,val])