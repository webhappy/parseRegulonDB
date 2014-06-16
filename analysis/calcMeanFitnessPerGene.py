import pandas
import numpy as np
import mysql.connector

# Read filtered data
# expData = pandas.read_csv('../Anaerobic - filtered over 60.csv')
expData = pandas.read_csv('../aerobic.csv')
expData['mean_lr'] = expData[[elem for elem in expData.columns if elem[-3:]=='_LR']].mean(axis=1)
expDataGrouped = expData.groupby('RegulonDBID').agg({'mean_lr':np.mean})
print expDataGrouped

cnx = mysql.connector.connect(user='root', host='127.0.0.1',  database='t')
cursor = cnx.cursor(buffered=True)
cursor.execute('select gene_id, gene_name from GENE')
cursorT=cnx.cursor(buffered=True)
genes={}
for (id,name) in cursor:
    try:
        val=float(expDataGrouped.loc[id][0])
    except:
        print name,id,'has no records'
        continue

    cursorT.execute('update GENE set avg_lr_aerobic=%s where gene_id=%s',[val,id])