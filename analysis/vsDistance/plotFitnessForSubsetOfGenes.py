# Get all genes w/ specified operon position
# Calculate fitness for each gene by averaging all sgRNA's in opposite strand within the first 100 bp (check it's not further than gene size)
import csv
import mysql.connector
import cPickle
from Bio import SeqIO
import numpy
import matplotlib

anaerobicFile=csv.reader(open('../../../Anaerobic.csv','rb'))
header=anaerobicFile.next()
#sgRNA_pos,seq,sgRNA_strand,t8_1_LR,t8_2_LR,t8_3_LR
anaerobicExptData={'pos':[],'strand':[],'replicate1':[],'replicate2':[],'replicate3':[],'all':[]}
for row in anaerobicFile:
    anaerobicExptData['pos'].append(int(row[0]))
    strandProper='forward' if row[2]=='+' else 'reverse'
    anaerobicExptData['strand'].append(strandProper)
    anaerobicExptData['replicate1'].append(float(row[3]))
    anaerobicExptData['replicate2'].append(float(row[4]))
    anaerobicExptData['replicate3'].append(float(row[5]))
    anaerobicExptData['all'].append([float(row[3]),float(row[4]),float(row[5])])

def getSgRNAByPosAndStrand (data,left, right, strand):
    ret=[]
    for j in range(len(data['pos'])):
        if data['strand'][j] == strand and left <= data['pos'][j] <= right:
            #this sgRNA should be kept
            ret.append((data['pos'][j],data['all'][j]))
    return ret

refFasta = SeqIO.read("../sequence.fasta", 'fasta')

cnx = mysql.connector.connect(user='root',
                              host='127.0.0.1',
                              database='t')
cursor = cnx.cursor(buffered=False)
posStr=1;
cursor.execute('SELECT gene_name, gene_posleft, gene_posright, gene_strand from GENE where operon_pos=%s && avg_lr < -3'%posStr)
vals=[];
#outFile=csv.writer(open('fitness_first_gene_in_operon.csv','wb'))
for (name,left,right,strand) in cursor:
    #get all sgRNA's that are within the range
    candidates=getSgRNAByPosAndStrand(anaerobicExptData,left,right,strand)
    if len(candidates) >0:
        allArray=[]
        for (pos,allVals) in candidates:
            allArray.extend(allVals)
        meanVal=sum(allArray) / float(len(allArray))
        vals.append(meanVal)
        outFile.writeRow([name,left,right,strand,meanVal])
        print '%s has mean %f'%(name,meanVal)