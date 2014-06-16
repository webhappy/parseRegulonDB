import mysql.connector
import pandas as pd
import pandas.io.parsers

cnx = mysql.connector.connect(user='davidc', password='mysql_password', host='127.0.0.1',  database='CRISPR')
df = pd.io.parsers.read_csv('../All data - 32992 with LRs.csv', index_col='seq')

print df

# first move in the two aerobic ones