from pandas.io.parsers import *
import numpy as np
import mysql.connector
from pandas.io.sql import *

df = read_csv(open('selected.csv'))
df.index = df['seq']
cnx = mysql.connector.connect(user='davidc', password='mysql_password', host='127.0.0.1',  database='CRISPR')
NUM_SGRNA = len(df.iloc[:,0])
proportions = (df.iloc[:, 1:] + 1).iloc[:,:].apply(lambda x:x/sum(x),1)  # Normalize to proportions

sgRNAs_with_gene = read_frame('select * from sgRNAs where gene_name!=" "', cnx, index_col='seq')
df = df.merge(sgRNAs_with_gene, 'left', left_index=True, right_index=True)

column_names = {'control':('1_aerobic_t30_1_control', '1_aerobic_t60_1_control', '1_aerobic_t180_1_control'),
              'chlor':('1_aerobic_t30_2_chlor', '1_aerobic_t60_2_chlor', '1_aerobic_t180_2_chlor'),
              'nor':('1_aerobic_t30_6_nor', '1_aerobic_t60_6_nor', '1_aerobic_t180_6_Nor')  }
conditions = column_names.keys()[1:]  # assume first entry at index 0 is the control

outFile = open('consistency_of_genes.csv','w')
outFile.write('Gene')
for c in conditions:
        for k in xrange(3):
            outFile.write(','+'mean_'+column_names[c][k]+', ' +'std_'+column_names[c][k] )
outFile.write('\n')

genes = df.gene_name.dropna().unique()
for gene in genes:
    outFile.write(gene)
    rows = df[df.gene_name==gene]
    for c in conditions:
        for k in xrange(3):
            cur = np.log2(rows[column_names[c][k]] / rows[column_names['control'][k]])
            outFile.write(", %.3f, %.3f"%(cur.mean(),cur.std()))
    outFile.write("\n")
outFile.close()