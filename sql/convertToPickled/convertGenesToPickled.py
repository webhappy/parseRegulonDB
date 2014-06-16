import mysql.connector
from Bio import SeqIO
import pickle

refFasta = SeqIO.read("sequence.fasta", 'fasta')

cnx = mysql.connector.connect(user='root',
                              host='127.0.0.1',
                              database='t')
PROM_BOUNDARY = 35;
cursor = cnx.cursor()
cursor.execute('SELECT gene_id,gene_name, gene_posleft, gene_posright, gene_strand from GENE')
ret = []
for (id,name, left, right,strand) in cursor:
    if left and right: #regulonDB includes some rows that are essentially empty
        id=str(id)
        name=str(name)
        left=int(left)
        right=int(right)
        strand=str(strand)
        left_boundary = left
        right_boundary = right

        assert (len(refFasta) >= left_boundary >= 0 and 0 <= right_boundary <= len(refFasta))
        ret.append((id,name, strand, left_boundary, right_boundary))


pickle.dump(ret, open('../convertToPickled/genes_pickled.txt', 'wb'))
cnx.close()