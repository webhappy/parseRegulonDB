import mysql.connector
from Bio import SeqIO
import pickle

refFasta = SeqIO.read("sequence.fasta", 'fasta')

cnx = mysql.connector.connect(user='root',
                              host='127.0.0.1',
                              database='t')

sql='select srna_id,srna_posleft,srna_posright, srna_sequence, gene_name,gene_strand from `SRNA_INTERACTION` left join GENE on `SRNA_INTERACTION`.`srna_gene_id`=GENE.`gene_id` where srna_posleft > 0 ORDER BY srna_posleft asc'
cursor = cnx.cursor()
cursor.execute(sql)
ret = []
for (id, left, right,seq,name,strand) in cursor:
    if left and right: #regulonDB includes some rows that are essentially empty
        name=str(name)
        left=int(left)
        right=int(right)
        strand=str(strand)
        left_boundary = left
        right_boundary = right

        assert (left_boundary<=right_boundary and len(refFasta) >= left_boundary >= 0 and 0 <= right_boundary <= len(refFasta))
        ret.append((id,name, strand, left_boundary, right_boundary))


pickle.dump(ret, open('../exptServer/../../exptServer/sRNA_pickled.txt', 'wb'))
cnx.close()