import csv
import mysql.connector
from numpy import polyfit
import cPickle
from Bio import SeqIO
# get genes from database with avg_lr under -3

cnx = mysql.connector.connect(user='root',
                              host='127.0.0.1',
                              database='t')

def getDistanceForSgRNAFromStart (sgRNApos, strand, gene_left, gene_right):
    if strand=='forward':
        return sgRNApos-gene_left
    else:
        return gene_right-(20+sgRNApos)

allPoints=[]
fits=[]
def getSgRNAByPosAndStrand (left, right, strand,geneName):
    """Flips the strand orientation inside here before getDistanceFromStart

    :param left:
    :param right:
    :param strand: The actual strand orientation the gene
    """
    cursorT=cnx.cursor(buffered=True)
    flipped=flipStrand(strand)
    cursorT.execute("select pos,avg(log2) from EXPT_RESULTS "+
                    "where strand=%s and pos<=%s and pos>=%s GROUP BY pos, strand",(flipped,right,left))
    items=[]
    minVal=0
    minDist=10^6
    maxDist=0
    for (pos,val) in cursorT:
        dist=getDistanceForSgRNAFromStart(pos,strand,left,right)
        items.append((pos,val,dist))
        if val < minVal:
            minVal=val
        minDist=min(minDist,dist)
        maxDist=max(maxDist,dist)

    if len(items) <3:
        print geneName,'does not have enough sgRNA\'s'
        return
    if maxDist-minDist>50:
        for (pos,val,dist) in items:
            allPoints.append((geneName,dist,val/minVal))
        (b,a)=polyfit([int(i[2]) for i in items],[i[1]/minVal for i in items],1)
       # print 'gene ',geneName,' has slope ',b,' intercept ',a
        fits.append((geneName,a,b))
       # print pos,'has',val,'w/ dist',getDistanceForSgRNAFromStart(pos,strand,left,right)

def flipStrand (strand):
    if strand=='forward':
        return 'reverse'
    elif strand=='reverse':
        return 'forward'
    else:
        raise Exception('strand doesn\'t match!')

outF=csv.writer(open('fitnessVsDist.csv','w'))
fitsFile=csv.writer(open('fittingPerGene.csv','w'))
cursor = cnx.cursor(buffered=True)
posStr=1;
cursor.execute('SELECT gene_name, gene_posleft, gene_posright, gene_strand from GENE where operon_pos=%s && avg_lr < -3'%posStr)
vals=[];
iterCount=0
#outFile=csv.writer(open('fitness_first_gene_in_operon.csv','wb'))
for (name,left,right,strand) in cursor:
    #get all sgRNA's that are within the range
    getSgRNAByPosAndStrand(left,right,strand,name)
    iterCount+=1

print len(allPoints)
print len(fits)
for (name,dist,val) in allPoints:
    outF.writerow([name,dist, val])

for (name,a,b) in fits:
    fitsFile.writerow([name,a,b])
