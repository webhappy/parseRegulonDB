import mysql.connector
from Bio import SeqIO
import pickle
import csv

cnx = mysql.connector.connect(user='root',
                              host='127.0.0.1',
                              database='t')
cursor = cnx.cursor(buffered=True)
anaerobicFile=csv.reader(open('../Anaerobic - filtered over 60.csv','rb'))
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

    cursor.execute('INSERT INTO EXPT_RESULTS (condition_name,replicate,pos,strand,log2) VALUES (%s,%s,%s,%s,%s)',('anaerobic',1,int(row[0]),strandProper,float(row[3])))
    cursor.execute('INSERT INTO EXPT_RESULTS (condition_name,replicate,pos,strand,log2) VALUES (%s,%s,%s,%s,%s)',('anaerobic',2,int(row[0]),strandProper,float(row[4])))
    cursor.execute('INSERT INTO EXPT_RESULTS (condition_name,replicate,pos,strand,log2) VALUES (%s,%s,%s,%s,%s)',('anaerobic',3,int(row[0]),strandProper,float(row[5])))
    cnx.commit();