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
    else:
        return 3

def median (items):
    items=sorted(items)
    length=len(items)
    if length == 0:
        return 0
    if (length % 2) ==0: #even number
        return (items[length/2]+items[length/2-1])/2.0
    else:
        return items[length/2]

cnx = mysql.connector.connect(user='root',host='127.0.0.1',database='t')
data=pandas.io.parsers.read_csv('../Anaerobic - filtered over 60.csv')

outFile=csv.writer(open('heatmap_test.csv','wb'))

allData=[[] for n in xrange(4)]
badData=[[] for n in xrange(4)]
for n in xrange(4):
    allData[n]=[[] for j in xrange(20)]
    badData[n]=[[] for j in xrange(20)]
cursor = cnx.cursor(buffered=True)
posStr=1;
seen={}
base='A'*20

items=[('AAAAAAAAAAAAAAAAAAAA',-10),('AAAAAAAAAAAAAAAAAAAG',-1),('AAAAAAAAAAAAAAAAAAAC',-6),('CAAAAAAAAAAAAAAAAAAA',-100)]
for (seq,val) in items:
    if seq in seen: #Notice some sgRNA's are hit by multiple genes
        continue
    else:
        seen[seq]=True
    seq=seq.upper()

    for k in xrange(20):
        allData[convertNTtoInt(seq[k])][k].append(float(val))

for n in xrange(4):
    outrow=[]
    for k in xrange(20):
        outrow.append(median(allData[n][k]))
    outFile.writerow(outrow)


