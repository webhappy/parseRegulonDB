import mysql.connector
from Bio import SeqIO
import pickle

refFasta = SeqIO.read("sequence.fasta", 'fasta')

cnx = mysql.connector.connect(user='root',
                              host='127.0.0.1',
                              database='t')
PROM_BOUNDARY = 35;
cursor = cnx.cursor()
cursor.execute('SELECT promoter.promoter_id, promoter_name,promoter_strand,pos_1,promoter_sequence,box_35_left,box_35_right '
               'from `PROMOTER` LEFT JOIN `PROMOTER_FEATURE` ON PROMOTER.`promoter_id`=PROMOTER_FEATURE.`promoter_id`'
               ' where pos_1 != \'\'')
ret = []
for (id,name, strand, TSS, seq, box35_left, box35_right) in cursor:
    TSS=int(TSS)
    if box35_left and box35_right:
        box35_right=int(box35_right)
        box35_left=int(box35_left)
    left_boundary = right_boundary = 0
    if strand == 'forward':
        if not box35_left:  # if empty
            right_boundary = TSS+10
            #if TSS < PROM_BOUNDARY:  #This complicated logic was to handle cases where the TSS was near 0
            #    left_boundary = len(refFasta) - 35
            #else:
            left_boundary = TSS - 60
        else:  # box35 is defined
            right_boundary = TSS
            left_boundary = box35_left
    else:  # strand is reverse
        if not box35_right:  # if empty
            left_boundary = TSS - 10
            right_boundary = 60 + TSS
           # if right_boundary > len(refFasta):
           #     right_boundary -= len(refFasta)
        else:
            left_boundary = TSS
            right_boundary = box35_right

    if right_boundary > len(refFasta):
        right_boundary=len(refFasta)
    if left_boundary < 1:
        left_boundary=1
#    assert(left_boundary<=right_boundary)  # These errors are in regulonDB's fault!
    assert (len(refFasta) >= left_boundary >= 0 and 0 <= right_boundary <= len(refFasta))
    ret.append((id, name, strand, TSS, seq, left_boundary, right_boundary))


pickle.dump(ret, open('../convertToPickled/promoters_pickled.txt', 'wb'))
cnx.close()