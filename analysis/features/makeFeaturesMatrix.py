import csv
import cPickle

import mysql.connector
import pandas.io.parsers

import runNUPACK
from analysis.features.calcBindingProperties import *
from GenericUtils import *




#parse fitnessVsDist_badcases.csv and find rows with 4th entry
cnx = mysql.connector.connect(user='root',host='127.0.0.1',database='t')



def getSgRNASeq (data, coord, strand):
    if strand.lower()=='forward':
        strand='+'
    else:
        strand='-'
    t=data[ (data.sgRNA_pos==coord) & (data.sgRNA_strand==strand) ]
    assert len(t)==1

    return t.iloc[0]['seq']

def calcBindingProbabilities(seq):
    probs = runNUPACK.callNuPack(seq)
    sumAll = 0
    for i in xrange(20):
        for j in xrange(20, len(seq)):
            #Add up probability of interacting for each base to a base outside the first 20
            sumAll += probs[i][j]
    sumSeed = 0  # Sum of probability of seed region interacting with any nucleotide outside the seed region
    restrictedRange = range(15, 20)
    for i in restrictedRange:
        for j in xrange(len(seq)):
            if j in restrictedRange:
                continue
            sumSeed += probs[i][j]
    return sumAll / 20.0, sumSeed /5.0

posWeights=[0,0,.014,0,0,0.395,0.317,0,0.389,0.079,0.445,0.508,0.613,0.851,0.732,0.828,0.615,0.804,0.685,0.583]
def getNucleotideWeight (seq, nt):
    seq=seq.lower()
    nt=nt.lower()
    sum=0.0
    for j in xrange(len(seq)):
        if seq[j]==nt:
            sum+=posWeights[j];
    return sum

offTargetsData=cPickle.load(open('../analysis/offTargets/pickled.txt'))

cursor = cnx.cursor(buffered=True)
badSgRNAAutomatic=[]
goodSgRNA=[]
lastGene=''
lastDist=10^6
lastVal=10^6
outFile=csv.writer(open('features.csv','wb'))
outHeader=['GeneName','Target','entireSeq','numOffTargets','GCperc','DistanceFromATG','probAll','probSeed','bindingEnergy']
outHeader.extend(['wt'+j for j in ['A','C','G','T']])
outFile.writerow(outHeader)

seen={}
data=pandas.io.parsers.read_csv('../Anaerobic - filtered over 60.csv'   )
cursor.execute('SELECT gene_name, gene_posleft, gene_posright, gene_strand,operon_pos from GENE where avg_lr < -3')
genesCount=0
for (name,left,right,strand,operon_pos) in cursor:
    genesCount+=1
    items,minVal=getSgRNAByPosAndStrand(left,right,strand)
    for (val,dist,seq) in items: # Foreach sgRNA
        if seq in seen: #Notice some sgRNA's are hit by multiple genes
            continue
        else:
            seen[seq]=True

        if strand == 'forward':
            coord=dist+int(left)
        else:
            coord=int(right)-dist-20
        sgRNAseq=getSgRNASeq(data,coord,flipStrand(strand))
        probAll,probSeed=calcBindingProbabilities(sgRNAseq+"GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGT")
        outLine=[name,val,sgRNAseq+ 'GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGT',len(offTargetsData[sgRNAseq]) if sgRNAseq in offTargetsData else 0,getGCPerc(sgRNAseq),dist,probAll,probSeed]
        outLine.append(getEnergyForSgRNA(sgRNAseq))
        outLine.extend([getNucleotideWeight(sgRNAseq,j) for j in ['A','C','G','T']])
        outFile.writerow(outLine)
        print name,dist,val
