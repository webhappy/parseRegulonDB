import pandas.io.sql as psql
import numpy as np
import csv, mysql.connector
from GenericUtils import *

# This script looks for sgRNA's that I suspect are bad based on experimental results
# The proposed workflow is to take the output of this and feed into the BLAST searcher to see if there's a pattern
# Foreach gene
# If there are at least 3 sgRNA's and one is strong while the other two have weak effect
# count this one as a defective sgRNA

badCasesFile = csv.writer(open('badSgRNAs.csv','wb'))
badCasesFile.writerow(['Gene','Seq','Position','Strand','# sgRNAs on gene','pVal','log fold-change'])

def getSgRNAByPosAndStrand (left, right, strand,geneName):
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
    cursorT.execute("select s.pos, s.seq,p.pval,p.logFC, sum as countsSum from sgRNA s LEFT JOIN anaerobic_PVALS p on s.seq=p.seq"+
                    " where strand=%s and s.pos<=%s and s.pos>=%s",(flipped,right,left))
    items=[]
    maxVal=0
    maxSGRNA=()
    for (pos,seq,pval,logFC, countsSum) in cursorT:
        try:
            pval = float(pval)
            logFC = float(logFC)
            #dist=getDistanceForSgRNAFromStart(pos,strand,left,right)
            items.append((logFC,pval,countsSum))
            if pval < .05 :
                if logFC < maxVal:
                    maxVal = logFC
                    maxSGRNA = (seq,pos,pval,logFC)
        except:
            print seq,'has no pVal'

    DEFECTIVE_CUTOFF = -2.0
    if len(items) > 0:
        goodNeutralCount = 0
        growthDefectiveCount = 0
        for (logFC,pval,countsSum) in items:
            if logFC < DEFECTIVE_CUTOFF and pval < .05:
                growthDefectiveCount += 1
            if logFC > maxVal: # NOT the smallest fold change for this gene
                if countsSum > 10 and logFC > -1: #CUTOFF for an sgRNA supporting this gene being good
                    goodNeutralCount += 1

        if growthDefectiveCount==1 and maxVal < DEFECTIVE_CUTOFF and goodNeutralCount >= 2:
            print "Bad sgRNA found at",geneName
            badCasesFile.writerow([geneName,maxSGRNA[0],maxSGRNA[1],strand,len(items),maxSGRNA[2],maxSGRNA[3]])


cnx = mysql.connector.connect(user='root', host='127.0.0.1',  database='t')
cursor = cnx.cursor(buffered=True)
cursor.execute('select gene_id, gene_name, gene_posleft, gene_posright, gene_strand from GENE')
cursorT=cnx.cursor(buffered=True)
genes={}
for (gene_id,name,left,right,strand) in cursor:
    getSgRNAByPosAndStrand(left,right,strand,name)
