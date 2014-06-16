import mysql.connector
from GenericUtils import *
import csv

outfile=csv.writer(open('overlapping_sgRNAs.csv','wb'))
outfile.writerow(['Gene_name','Distance1','Distance2','delta distance','logFC','logFC2'])
def getSgRNAByPosAndStrand (left, right, strand,name):
    """Flips the strand orientation inside here before getDistanceFromStart

    :param left:
    :param right:
    :param strand: The actual strand orientation of the gene
    """
    cnx = mysql.connector.connect(user='root',host='127.0.0.1',  database='t')
    cursorT=cnx.cursor(buffered=True)
    flipped=convert_strand_to_minus_plus(flipStrand(strand))
    cursorT.execute("select s.pos, s.seq,p.pval,p.logFC from sgRNA s LEFT JOIN aerobic_PVALS p on s.seq=p.seq"+
                    " where strand=%s and s.pos<=%s and s.pos>=%s order by s.pos asc",(flipped,right,left))
    items=[]
    maxVal=0
    last = None
    lastDist = None
    for (pos,seq,pval,logFC) in cursorT:
        pval = float(pval)
        logFC = float(logFC)
        distFromATG = getDistanceForSgRNAFromStart(pos, strand, left, right)
        if last is not None:
            if abs(logFC-last) > 2 and distFromATG-lastDist <= 10:  # Print out this sgRNA and next's
                outfile.writerow([name, lastDist,distFromATG,distFromATG-lastDist,last,logFC])
        last=logFC
        lastDist = distFromATG

cnx = mysql.connector.connect(user='root', host='127.0.0.1',  database='t')
cursor = cnx.cursor(buffered=True)
cursor.execute('select gene_id, gene_name, gene_posleft, gene_posright, gene_strand from gene')
genes={}
for (gene_id,name,left,right,strand) in cursor:
    getSgRNAByPosAndStrand(left,right,strand,name)