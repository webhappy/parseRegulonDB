#From genes that are essential, pick out the sgRNAs that are not working
import pandas
import numpy as np
import mysql.connector
from GenericUtils import *
import csv

outfile=csv.writer(open('defectiveSgRNAs.csv','wb'))
outfile.writerow(['seq','Gene_name','Distance from atg','logFC','p-Val','Gene mean','# sgRNAs'])
def getSgRNAByPosAndStrand (left, right, strand,name):
    """Flips the strand orientation inside here before getDistanceFromStart

    :param left:
    :param right:
    :param strand: The actual strand orientation of the gene
    """
    cnx = mysql.connector.connect(user='root',
                              host='127.0.0.1',
                              database='t')
    cursorT=cnx.cursor(buffered=True)
    flipped=convert_strand_to_minus_plus(flipStrand(strand))
    cursorT.execute("select s.pos, s.seq,p.pval,p.logFC from sgRNA s LEFT JOIN aerobic_PVALS p on s.seq=p.seq"+
                    " where strand=%s and s.pos<=%s and s.pos>=%s",(flipped,right,left))
    items=[]
    maxVal=0
    for (pos,seq,pval,logFC) in cursorT:
        try:
            pval = float(pval)
            logFC = float(logFC)
            #dist=getDistanceForSgRNAFromStart(pos,strand,left,right)
            items.append((pos,seq,pval,logFC))
            if abs(maxVal) < abs(logFC):
                maxVal = logFC
        except:
            print seq,'has no pVal'

    if len(items) == 0:
        return
    allFC = [x[-1] for x in items]  # Pull out the logFC in the final position of each item
    mean_of_list = meanOfList(allFC)
    if mean_of_list > 2:  # Look for an sgRNA that's weak
        for (pos,seq,pval,logFC) in items:
            if abs(logFC) < 1:
                outfile.writerow([seq,name,getDistanceForSgRNAFromStart(pos,strand,left,right),logFC,pval,mean_of_list,len(items)])

cnx = mysql.connector.connect(user='root', host='127.0.0.1',  database='t')
cursor = cnx.cursor(buffered=True)
cursor.execute('select gene_id, gene_name, gene_posleft, gene_posright, gene_strand from gene where avg_lr_aerobic > 1.5')
genes={}
for (gene_id,name,left,right,strand) in cursor:
    getSgRNAByPosAndStrand(left,right,strand,name)