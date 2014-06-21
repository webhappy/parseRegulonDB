import pandas.io.sql as psql
import numpy as np
import csv, mysql.connector

con=mysql.connector.connect(user='root', host='127.0.0.1', database='t')
genes=psql.read_frame('select gene_name,gene_strand,gene_posleft,gene_posright,operon_pos,extreme as AVG_LR from GENE',con)
#genes = psql.read_frame("select gene_name,gene_strand,gene_posleft,gene_posright,operon_pos, AVG(aerobic_t8_4_LR+aerobic_t8_5_LR) as avg_lr FROM gene LEFT JOIN aerobic_data ON aerobic_data.RegulonDBID=GENE.gene_id GROUP BY gene_id",con)
print "Big join SQL query done!"
PECgenes=psql.read_frame('select gene_name from PEC',con)
PECgenes=PECgenes['gene_name'].values

tails=[]  # Append pointers to tails here
for strand in ['forward','reverse']:
    last=None
    genesInDir = genes[(genes['gene_strand']==strand)]
    if strand == 'forward':
        genesInDir.sort('gene_posleft',ascending=True,inplace=True)
    else:
        genesInDir.sort('gene_posright',ascending=False,inplace=True)
    for (gene,strand,left,right,pos,avgLR) in genesInDir.itertuples(index=False):
        if gene =='' or gene=='None':
            continue
        if last is None:
            count=1
            last={'gene':gene,'val':avgLR,'prev':None}
        elif pos==1:
            #print "%s has count %i"%(gene,count)
            tails.append(last)
            last={'gene':gene,'val':avgLR,'prev':None}
            count=1
        else:
            count+=1
            cur={'gene':gene,'prev':last,'val':avgLR}
            last=cur
    tails.append(last)

def getPrecisionRecall (pointersToTails,cutoff,knownEssentials):
    determinedSick=[]
    determinedNonSick=[]
    undetermined=[]
    unexpectedHealthy = []  # Genes that are healthy despite someone downstream being sick
    for pointer in pointersToTails:
        cur=pointer
        canDetermine=True
        while cur is not None:
            if cur['val'] > cutoff:
                determinedNonSick.append(cur['gene'])
                if not canDetermine:
                    unexpectedHealthy.append(cur['gene'])
            else:
                if not canDetermine:
                    undetermined.append(cur['gene'])
                else:
                    determinedSick.append(cur['gene'])
                    canDetermine=False
            cur=cur['prev']

    TP=0;  # True positives
    missingEssentials=0;  # False negative
    undeterminedEssential=0;

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
    return TP*1.0/len(knownEssentials), TN*1.0/knownNumNonEssential, TP*1.0/len(determinedSick)


    # tp=fp=tn=fn=0
    # for gene in determinedSick:
    #     if gene in knownEssentials:
    #         tp+=1
    #     else:
    #         fp+=1
    # for gene in determinedNonSick:
    #     if gene in knownEssentials:
    #         fn += 1
    #     else:
    #         tn += 1
    # undeterminedEssentials=undeterminedNonessentials=0
    # for gene in undetermined:
    #     if gene in knownEssentials:
    #         undeterminedEssentials += 1
    #     else:
    #         undeterminedNonessentials += 1
    # print "TP=%i, FP=%i, TN=%i, FN=%i, undetermined essentials=%i, undetermined non-essentials=%i"%(tp,fp,tn,fn,undeterminedEssentials,undeterminedNonessentials)

outFile=csv.writer(open('essentialsPerformance.csv','wb'))
outFile.writerow(['Cutoff','Recall','Specificity','Precision'])
for cutoff in np.arange(0,-5,-.1):
    (recall,specificity,precision)=getPrecisionRecall(tails,cutoff,PECgenes)
    outFile.writerow([cutoff,recall,specificity,precision])