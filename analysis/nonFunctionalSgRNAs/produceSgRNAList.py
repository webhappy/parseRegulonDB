import mysql.connector
from GenericUtils import *
import csv

# For each gene in class, print out relevant stuff
genesWithBadSgRNA = ['pheU','rrfG','cyaR','birA','yhhX']
overlappingSgRNA = ['icd','mtr','mepM','yigB','yghE','yqjH']
genesWithOfftarget = ['mdtf','ilvH','hdeD','caic','eutq']

#Print out gene name, position, sgRNA seq, logFC
outFile = csv.writer(open('sgRNAs_to_design.csv','wb'))
outFile.writerow(['Gene name','Distance from ATG','sgRNA seq','logFC','Sum of starting counts'])

def getSgRNAByPosAndStrand (left, right, strand,name):
    """Flips the strand orientation inside here before getDistanceFromStart

    :param left:
    :param right:
    :param strand: The actual strand orientation of the gene
    """
    cnx = mysql.connector.connect(user='root',host='127.0.0.1',  database='t')
    cursorT=cnx.cursor(buffered=True)
    flipped=convert_strand_to_minus_plus(flipStrand(strand))
    cursorT.execute("select s.pos, s.seq,p.pval,p.logFC, p.countsSum from sgRNA s LEFT JOIN aerobic_PVALS p on s.seq=p.seq"+
                    " where strand=%s and s.pos<=%s and s.pos>=%s order by s.pos "+("asc" if strand=='forward' else "desc"),(flipped,right,left))
    items=[]
    maxVal=0
    last = None
    lastDist = None
    for (pos,seq,pval,logFC,countsSum) in cursorT:
        pval = float(pval)
        logFC = float(logFC)
        distFromATG = getDistanceForSgRNAFromStart(pos, strand, left, right)
        outFile.writerow([name,distFromATG,seq,logFC,countsSum])

cnx = mysql.connector.connect(user='root', host='127.0.0.1',  database='t')
cursor = cnx.cursor(buffered=True)
allGenes = genesWithBadSgRNA
allGenes.extend(overlappingSgRNA)
allGenes.extend(genesWithOfftarget)
for gene_name in allGenes:
    cursor.execute('select gene_id, gene_name, gene_posleft, gene_posright, gene_strand from gene where gene_name=%s',(gene_name,))
    for (gene_id,name,left,right,strand) in cursor:
        getSgRNAByPosAndStrand(left,right,strand,name)