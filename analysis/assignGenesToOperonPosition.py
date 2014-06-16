import mysql.connector
from Bio import SeqIO
import pandas.io.sql as psql
import pickle

refFasta = SeqIO.read("../sequence.fasta", 'fasta')

cnx = mysql.connector.connect(user='root',
                              host='127.0.0.1',
                              database='CRISPR')
cursor = cnx.cursor(buffered=True)
cursor.execute('SELECT operon_id, operon_name, firstgeneposleft,lastgeneposright,operon_strand from OPERON')
promoters = psql.read_frame('select transcription_unit_name, operon_id, pos_1 from transcription_unit LEFT JOIN promoter on transcription_unit.promoter_id=promoter.promoter_id',cnx)

ret = []
genePos={}  # Store found values to here
for (operon_id, name, left, right,strand) in cursor:
    geneCursor=cnx.cursor(buffered=True)
    if strand=='forward':
        orderBy='asc'
    else:
        orderBy='desc'
    promoters_in_operon = promoters[promoters.operon_id==operon_id].dropna()
    geneCursor.execute(
        "SELECT gene_name from GENE where gene_posleft>=%s and gene_posright <=%s and gene_strand='%s' ORDER BY gene_posleft %s"%(left,right,strand,orderBy))
    curPos = 1
    promoters_passed = 0
    for geneName in geneCursor:
        if strand == 'forward':
            promoters_upstream = promoters_in_operon.pos_1 < left
        else:
            promoters_upstream = promoters_in_operon.pos_1 > right
        promoters_upstream = promoters_upstream.shape[0]
        if promoters_upstream > promoters_passed:
            curPos = 1
            promoters_passed = promoters_upstream
        geneName=str(geneName[0])
        if geneName in genePos and genePos[geneName]!=0:
            if genePos[geneName] != curPos:
                print geneName,' has conflict! setting to 0'
                genePos[geneName]=0
        else:
            genePos[geneName]=curPos
            #print geneName,'to',curPos
        curPos+=1

for (k,v) in genePos.iteritems():
    #query="UPDATE GENE SET operon_pos=%s where gene_name='%s'"%(v,str(k))
    #print query
    cursor.execute("UPDATE GENE SET operon_pos=%s where gene_name=%s",(v,k))
    cnx.commit()

'''READ ME BELOW'''
# Need to manually set remaining nulls to 1