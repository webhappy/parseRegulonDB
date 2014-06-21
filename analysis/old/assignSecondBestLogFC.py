#For each gene, pick out the sgRNA with the second-strongest log FC effect
import pandas
import numpy as np
import mysql.connector
from GenericUtils import *
import csv

DIST_FROM_ATG_CUTOFF = 200  # Can't just use the gene label from All data since some sgRNAs refer to the gene with a different name

outFile = csv.writer(open('second_best_sgRNA.csv','wb'))
outFile.writerow(['Gene_name','Pval','LogFC','Position','Sequence','Num_sgRNAs'])

cnx = mysql.connector.connect(user='root', host='127.0.0.1',  database='t')
cursor = cnx.cursor(buffered=True)
cursor.execute('update gene set second_best=NULL')
cursor.execute('select gene_id, gene_name, gene_posleft, gene_posright, gene_strand from GENE')
cursorT=cnx.cursor(buffered=True)
genes={}
for (gene_id,name,left,right,strand) in cursor:
    cursorT=cnx.cursor(buffered=True)
    if strand == 'forward':
        # it's on the positive orientation, so want to let right-side by smaller of right and left+100
        right = min(right,DIST_FROM_ATG_CUTOFF + left)
    else:
        left = max(left,right - DIST_FROM_ATG_CUTOFF)

    flipped=convert_strand_to_minus_plus(flipStrand(strand))
    cursorT.execute("select s.pos, s.seq,p.pval,p.logFC, p.countsSum from sgRNA s LEFT JOIN aerobic_PVALS p on s.seq=p.seq"+
                    " where strand=%s and s.pos<=%s and s.pos>=%s",(flipped,right,left))

    items = []
    for (pos,seq,pval,logFC, countsSum) in cursorT:
        try:
            pval = float(pval)
            logFC = float(logFC)
            #dist=getDistanceForSgRNAFromStart(pos,strand,left,right)
            items.append((pos,seq,countsSum,pval,logFC))
        except:
            print seq,'has no pVal'

    allItems = items

    if len(items) == 0:
        print name,'has no sgRNAs'
        continue
    else:
        items = [x for x in allItems if x[-2] < .1]
        if len(items) > 0:
            if len(items) == 1:
                val = items[0][-1]
            else:
                ordered = sorted(items,key=lambda item: item[-1])
                val = ordered[1][-1]
        else:  # No sgRNA's pass the p-value test
            items = [x for x in allItems if x[2] > 15]  # filter for items with countsSum > 15
            if len(items) == 0:
                print name,'has no usable sgRNAs'
                continue
            else:
                ordered = sorted(items,key=lambda item: item[-2])  # sort by p-val if we're filtering by counts sum
                val = ordered[min(1,len(items)-1)][-1]  # if there's only one item that passes, take that; else, take second-best

    cursorN = cnx.cursor(buffered=True)
    cursorN.execute('update GENE set smart_second_best=%s where gene_id=%s',(val, gene_id) )