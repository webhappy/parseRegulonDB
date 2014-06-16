import csv
import mysql.connector
import cPickle
from Bio import SeqIO
import numpy
import matplotlib,pandas
from GenericUtils import *

expData=pandas.io.parsers.read_csv(open('../../Anaerobic - filtered over 60.csv','rb'))
expData['mean']=expData[[elem for elem in expData.columns if elem[-3:]=='_LR']].mean(axis=1)


cnx = mysql.connector.connect(user='root',
                              host='127.0.0.1',
                              database='t')
cursor = cnx.cursor(buffered=True)
posStr=1;
cursor.execute('SELECT gene_name, gene_posleft, gene_posright, gene_strand,operon_pos from GENE where avg_lr < -3')
vals=[];

outFile=csv.writer(open('sgRNAs_in_essential_genes.csv','wb'))
outFile.writerow(['avg_lr','normalized_lr','distance_ATG','operon_pos','geneName'])
geneCount = 0
genes = {}
for (name,left,right,strand,operon_pos) in cursor:
    items,minVal=getSgRNAByPosAndStrand(left,right,strand)
    for (val,dist,seq) in items:
        outFile.writerow([val,val/minVal,dist,operon_pos,name])
    geneCount += 1
    name=name.lower()
    if name in genes:
        print name,'was already seen before!'
    genes[name] = True

print geneCount,'genes found'
