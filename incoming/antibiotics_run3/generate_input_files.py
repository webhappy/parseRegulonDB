from pandas.io.parsers import *
import numpy as np
import mysql.connector
from pandas.io.sql import *

cnx = mysql.connector.connect(user='davidc', password='mysql_password', host='127.0.0.1',  database='CRISPR')
cursor = cnx.cursor()

#mappings: id_R1, id_R2, Description
mappings_file = open('mapping.csv', 'w')
mappings_file.write('id_R1, id_R2, Description\n')

#barcodes: orientation,sequence,id
barcodes_file = open('barcodes.csv', 'w')
barcodes_file.write('orientation,sequence,id\n')

cursor.execute('SELECT sampleID, description, barcode from samples where runID=4')  # NOTICE I've hardcoded in 4 here!
for (sampleID, description, barcode) in cursor:
    sanitized_description = str(description).replace(' ','_')
    mappings_file.write('%s,0,%s\n'%(sampleID, sanitized_description))
    barcodes_file.write('r1,%s,%s\n'%(barcode,sampleID))