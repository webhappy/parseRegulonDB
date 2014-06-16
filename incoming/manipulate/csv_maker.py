import mysql.connector
import csv

cnx = mysql.connector.connect(user='davidc', password='mysql_password', host='127.0.0.1',  database='CRISPR')

# input is a list of sampleIDs
sampleIDs = []

cursor = cnx.cursor()
cursor.execute('SELECT sampleID, samples.description, date FROM samples, runs where samples.runID=runs.runID and samples.sampleID IN (27,28,29,30)')
#cursor.execute('SELECT sampleID, samples.description, date FROM samples, runs where samples.runID=runs.runID and runs.runID IN (1,2,3,4)')
for items in cursor:
    sampleIDs.append(items[0])


# For each sample in samples, go grab the counts
sample_order = {}  # Map sampleID to 0-indexed position
for x in xrange(len(sampleIDs)):
    sample_order[ sampleIDs[x] ] = x

sample_names = []  # List of strings in same order as samples
for sampleID in sampleIDs:
    cursor.execute('select description from samples where sampleID=%s',(sampleID,))
    sample_names.append(str(cursor.next()[0]))

outFile = csv.writer(open('../../cluster/selected_aerobic.csv','wb'))
temp = ['seq']
temp.extend(sample_names)
outFile.writerow(temp)

sgRNA_cursor = cnx.cursor(buffered=True)
sgRNA_cursor.execute('select seq from sgRNAs')
for (seq,) in sgRNA_cursor:
    cursor.execute('SELECT sampleID,count FROM counts where sgRNA_seq=\'%s\' and sampleID in (%s)'%(seq,','.join(str(x) for x in sampleIDs)))
    countsMap = {}
    for sampleID, count in cursor:
        countsMap[sampleID] = count

    row = [seq]
    for sampleID in sampleIDs:
        row.append(countsMap[sampleID])
    outFile.writerow(row)