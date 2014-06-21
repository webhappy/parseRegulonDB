#From genes that are essential, pick out the sgRNAs that are not working
import pandas
import numpy as np
import mysql.connector
from GenericUtils import *

def getSgRNAByPosAndStrand (left, right, strand):
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
            if pval < .05 :
                items.append(logFC)
                if abs(maxVal) < abs(logFC):
                    maxVal = logFC
        except:
            print seq,'has no pVal'

    if len(items) == 0:
        return 0,0,0
    else:
        return maxVal, medianOfList(items), meanOfList(items)

cnx = mysql.connector.connect(user='root', host='127.0.0.1',  database='t')
cursor = cnx.cursor(buffered=True)
cursor.execute('select gene_id, gene_name, gene_posleft, gene_posright, gene_strand from GENE')
cursorT=cnx.cursor(buffered=True)
genes={}
for (gene_id,name,left,right,strand) in cursor:
    val,median,mean = getSgRNAByPosAndStrand(left,right,strand)

    cursorN = cnx.cursor(buffered=True)
    cursorN.execute('update GENE set avg_lr_aerobic=%s where gene_id=%s',(mean,gene_id))