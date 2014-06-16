import mysql.connector
import csv
from pandas.io.parsers import *

sgRNAs = read_csv(open('../All data - 32992 with LRs.csv'))
cnx = mysql.connector.connect(user='davidc', password='mysql_password', host='127.0.0.1',  database='CRISPR')
cursor = cnx.cursor()

for j in xrange(sgRNAs.shape[0]):  # for each row
    row = sgRNAs.iloc[j,:]
    cursor.execute('update sgRNAs set category=%s, gene_name=%s, product_name=%s where seq=%s',
        (row['CategoryID'], row['gene_name'], row['product_name'], row['seq']))
    #print cursor.statement
cnx.commit()