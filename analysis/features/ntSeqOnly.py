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

cnx = mysql.connector.connect(user='root',host='127.0.0.1',database='t')
data=pandas.io.parsers.read_csv('../Anaerobic - filtered over 60.csv')

outFile=csv.writer(open('features.csv','wb'))
#Prep header row to print
outHeader=['GeneName','Target','entireSeq']
for j in xrange(1,21):
    for n in 'ACGT':
        outHeader.append(n+str(j))
outFile.writerow(outHeader)

cursor = cnx.cursor(buffered=True)
posStr=1;
seen={}
cursor.execute('SELECT gene_name, gene_posleft, gene_posright, gene_strand,operon_pos from GENE where avg_lr < -3')
for (name,left,right,strand,operon_pos) in cursor:
    items,minVal=getSgRNAByPosAndStrand(left,right,strand)
    for (val,dist,seq) in items:
        if seq in seen: #Notice some sgRNA's are hit by multiple genes
            continue
        else:
            seen[seq]=True
        if val/minVal < .3:
            target=0  # target = 0 denotes a bad sgRNA
        else:
            target=1
        row=[name,target,seq]
        row.extend(genNTFeatures(seq))
        outFile.writerow(row)

