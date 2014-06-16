import mysql.connector
import csv
import pandas
from parameters import *

inFile = csv.reader(open('counts.csv'))
headers = inFile.next()  # ignore headers
cols = headers[1:-1]  # first column is sgRNA_id and last column is sequence

cnx = mysql.connector.connect(user='davidc', password='mysql_password', host='127.0.0.1',  database='CRISPR')
cursor = cnx.cursor()

sample_names = []
sample_ids = {}  # map description to ID
for col in cols:
    if col != 'seq':
        sample_names.append(col)

print len(sample_names)
print sample_names

# Check that each sample_name actually corresponds
for sample in sample_names:
    cursor.execute("select sampleID from samples where runID=%s and description=%s",(RUN_ID,sample))
    for (sampleID,) in cursor:
        if sample in sample_ids:
            raise '%s is already assigned to sample_ID %i'% (sample, sample_ids[sample])
        sample_ids[sample] = sampleID
    if sample not in sample_ids:
        raise Exception('%s failed to get assigned a sample_ID using SQL=%s'%(sample,cursor.statement))

# Now insert the counts now that we've verified
temp_cursor = cnx.cursor()
row_count = 0
for row in inFile:
    row_count += 1
    if row_count % 1000 == 0:
        print row_count, 'rows done'
    seq = row[-1]  # NOTICE! We assume first column is the ref int for the sgRNA and all next columns are counts
    for i in range(len(sample_names)):
        sample = sample_names[i]
        temp_cursor.execute('INSERT into counts SET sampleID=%s, sgRNA_seq=%s, count=%s',(sample_ids[sample], seq, row[i+1]))
        #print temp_cursor.statement
cnx.commit()
