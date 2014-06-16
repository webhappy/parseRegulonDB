import pandas,mysql.connector,csv
from GenericUtils import *
# Foreach nucleotide sequence, generate 80-row binary indicator matrix

def genNTFeatures (seq):
    ret=[]
    seq=seq.lower()
    NT = ['a','c','g','t']
    for j in xrange(len(seq)):
        for n in NT:
            if seq[j]==n:
                ret.append(1)
            else:
                ret.append(0)
    return ret

def convertNTtoInt (nt):
    if nt=='A':
        return 0
    elif nt=='C':
        return 1
    elif nt=='G':
        return 2
    elif nt=='T':
        return 3

def median (items):
    items=sorted(items)
    length=len(items)
    if (length % 2) ==0: #even number
        return (items[length/2]+items[length/2-1])/2.0
    else:
        return items[length/2]

def mean (items):
    sum=0.0
    for j in items:
        sum+=j
    return sum / len(items)

cnx = mysql.connector.connect(user='root',host='127.0.0.1',database='t')
data=pandas.io.parsers.read_csv('../Anaerobic - filtered over 60.csv')

outFile=csv.writer(open('heatmap_mean.csv','wb'))

allData=[[] for n in xrange(4)]
badData=[[] for n in xrange(4)]
for n in xrange(4):
    allData[n]=[[] for j in xrange(20)]
    badData[n]=[[] for j in xrange(20)]
cursor = cnx.cursor(buffered=True)
posStr=1;
seen={}
cursor.execute('SELECT gene_name, gene_posleft, gene_posright, gene_strand,operon_pos from GENE where avg_lr < -3')
genesCount=0
for (name,left,right,strand,operon_pos) in cursor:
    genesCount+=1
    items,minVal=getSgRNAByPosAndStrand(left,right,strand)
    for (val,dist,seq) in items:
        if seq in seen: #Notice some sgRNA's are hit by multiple genes
            continue
        else:
            seen[seq]=True
        seq=seq.upper()
        # if val/minVal < .3:
        #     storeIn=badData
        # else:
        #     storeIn=goodData

        for k in xrange(20):
            allData[convertNTtoInt(seq[k])][k].append(float(val))
print genesCount,'genes found'
print len(seen),'sgRNAs found'

for n in xrange(4):
    outrow=[]
    for k in xrange(20):
        outrow.append(mean(allData[n][k]))
    outFile.writerow(outrow)


