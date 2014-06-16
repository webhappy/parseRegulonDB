import mysql.connector
import csv


RUN_ID = 2  # Update this to correct run in

#####


cnx = mysql.connector.connect(user='davidc', password='mysql_password', host='127.0.0.1',  database='CRISPR')
cursor = cnx.cursor()

inFile = csv.reader(open('mapping.csv'))
inFile.next()  # skip the header row

for (id_R1, id_R2, description ) in inFile:
    cursor.execute('insert into samples SET description=%s, runID=%s',(description,RUN_ID))
    print cursor.statement
    cnx.commit()