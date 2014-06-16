import csv
import mysql.connector
import cPickle
from Bio import SeqIO

cnx = mysql.connector.connect(user='root',
                              host='127.0.0.1',
                              database='t')
#sgRNA_pos,seq,sgRNA_strand,t8_1_LR,t8_2_LR,t8_3_LR
cursor=cnx.cursor(buffered=True)
cursor.execute('select gene_id, gene_name, avg_lr from GENE')
genes={}
for (regulonID,name,avgLR) in cursor:
    genes[regulonID]=(name,avgLR)

def findNearestGene (location, direction):
    """

    :param location:  starting chromosome location
    :param direction: right or left
    """
    sqlCursor=cnx.cursor(buffered=True)
    query='Select gene_name, avg_lr from GENE where '
    if direction == 'right':
        query=query+" gene_posleft > "+str(location)+" order by gene_posleft asc limit 1"
    else:
        query=query+" gene_posright < "+str(location)+" order by gene_posright desc limit 1"
    sqlCursor.execute(query)
    return sqlCursor.fetchone()

def isWithinGene (loc,buffer=10):
    sqlCursor=cnx.cursor(buffered=True)
    query='Select count(*) from GENE where gene_posleft<'+str(loc+10)+' and '+str(loc-10)+'<gene_posright'
    sqlCursor.execute(query)
    (t,)=sqlCursor.fetchone()
    return True if t>0 else False

def mean (a,b,c):
    return (a+b+c)/3.0

anaerobicFile=csv.reader(open('../../Anaerobic - filtered over 60.csv','rb'))
outFile=csv.writer(open('longRangeScores.csv','wb'))
outFile.writerow(['Regulon ID','position','strand','score','my LR','geneLeft','LRleft','geneRight','LRright'])
header=anaerobicFile.next()
bestScores={} # Assume scores can't be exactly the same
for (RegulonDBID,sgRNA_pos,seq,strand,d1,d2,d3) in anaerobicFile:

    #RegulonDBID,sgRNA_pos,seq,sgRNA_strand,anaerobic_d24_1_LR,anaerobic_d24_2_LR,anaerobic_d24_3_LR
    avg=mean(float(d1),float(d2),float(d3))
    sgRNA_pos=int(sgRNA_pos)
    if RegulonDBID in genes:
        continue # this is a gene
    if isWithinGene(sgRNA_pos):
        continue  # we're not in an intergenic region
    score=0
    genes[RegulonDBID]=(sgRNA_pos,strand,avg)
    #find nearest gene to left
    try:
        nameLeft,lrLeft=findNearestGene(sgRNA_pos,'left')
        nameRight,lrRight=findNearestGene(sgRNA_pos,'right')
        score=lrLeft+lrRight + -2*avg
        #score=(lrLeft-avg)**2 + (lrRight-avg)**2 - (lrLeft-lrRight)**2
        bestScores[score]=(RegulonDBID,sgRNA_pos)
        print RegulonDBID,' has score ',score
        outFile.writerow([RegulonDBID,sgRNA_pos,strand,score,avg,nameLeft,lrLeft,nameRight,lrRight])
    except:
        print RegulonDBID, ' gave an error, skipping'

