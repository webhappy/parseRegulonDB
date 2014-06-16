import mysql.connector
import csv
import pandas

RUN_ID = 1

inFile = csv.reader(open('mapping_summary_1.csv'))#pandas.read_csv('mapping_summary_1.csv',header=1)
headers = inFile.next()  # ignore headers
cols = headers[1:1+7]  # 7 samples to parse

cnx = mysql.connector.connect(user='davidc', password='mysql_password', host='127.0.0.1',  database='CRISPR')
cursor = cnx.cursor()

db_fasta = open('../../OLS Master key - unique oligos 32992 with 5bp context.fas')
refs = {}  # Map of key to sgRNA + flanking nucleotides
refs_to_trimmed_sgRNA = {}  # Map of key to sgRNA 20-mer only
for l in db_fasta:
    if ">" in l:
        name = l[1:].strip()
    if not ">" in l:
        refs_to_trimmed_sgRNA[name] = l[6:26]
        refs[name] = l.strip()
        #print name, refs[name]
db_fasta.close()

# cols = 'count_OLS antibiotics 2 - aerobic t180 6_nor,count_OLS antibiotics 1 - aerobic t180 2_chlor,count_OLS antibiotics 1 - aerobic t60 6_nor,count_OLS antibiotics 1 - aerobic t60 2_chlor,count_OLS antibiotics 1 - aerobic t30 6_nor,count_OLS antibiotics 1 - aerobic t30 2_chlor,count_OLS antibiotics 2 - aerobic t180 2_chlor'
# cols = cols.split(',')
sample_names = []
sample_ids = {}  # map description to ID
for col in cols:
    sample_names.append(col[6:])

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
        raise '%s failed to get assigned a sample_ID',sample


# Now insert the counts now that we've verified
temp_cursor = cnx.cursor()
row_count = 0
for row in inFile:
    row_count += 1
    if row_count % 1000 == 0:
        print row_count, 'rows done'
    ref = row[0]  # NOTICE! We assume first column is the ref int for the sgRNA and all next columns are counts
    for i in xrange(len(sample_names)):
        sample = sample_names[i]
        temp_cursor.execute('INSERT into counts SET sampleID=%s, sgRNA_seq=%s, count=%s',(sample_ids[sample], refs_to_trimmed_sgRNA[ref], row[i+1]))
        #print temp_cursor.statement
        cnx.commit()
