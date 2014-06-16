import pandas.io.sql as psql
import numpy as np
import csv, mysql.connector

con=mysql.connector.connect(user='root', host='127.0.0.1', database='t')
genes=psql.read_frame('select gene_name,gene_strand,gene_posleft,gene_posright,operon_pos,median as avg_lr from GENE',con)

df = np.io.parsers.read_csv(open('selected.csv'))
df.index = df['seq']
cnx = mysql.connector.connect(user='davidc', password='mysql_password', host='127.0.0.1',  database='CRISPR')
NUM_SGRNA = len(df.iloc[:,0])
proportions = (df.iloc[:, 1:] + 1).iloc[:,:].apply(lambda x:x/sum(x),1)  # Normalize to proportions

sgRNAs_with_gene = psql.read_frame('select * from sgRNAs where gene_name!=" "', cnx, index_col='seq')
df = df.merge(sgRNAs_with_gene, 'left', left_index=True, right_index=True)

column_names = {'control':('1_aerobic_t30_1_control', '1_aerobic_t60_1_control', '1_aerobic_t180_1_control'),
              'chlor':('1_aerobic_t30_2_chlor', '1_aerobic_t60_2_chlor', '1_aerobic_t180_2_chlor'),
              'nor':('1_aerobic_t30_6_nor', '1_aerobic_t60_6_nor', '1_aerobic_t180_6_Nor')  }
conditions = column_names.keys()[1:]  # assume first entry at index 0 is the control



outFile = open('operon_consistency.csv', 'w')
temp = ['pos'+str(x) for x in range(1,17)]
outFile.write('first_gene,'+','.join(temp))
outFile.write('\n')

tails=[]  # Append pointers to tails here
operonLengths = {}  # Dict mapping name of tail gene


def write_cur_vals():
    line = [first]
    padded_vals = vals
    padded_vals.extend([' ']*(16-len(vals)))
    line.extend(padded_vals)
    outFile.write(','.join(line) + '\n')


for strand in ['forward','reverse']:
    first=None
    genesInDir = genes[(genes['gene_strand']==strand)]
    if strand == 'forward':
        genesInDir = genesInDir.sort('gene_posleft',ascending=True)
    else:
        genesInDir = genesInDir.sort('gene_posright',ascending=False)
    for (gene,strand,left,right,pos,avgLR) in genesInDir.itertuples(index=False):
        if gene =='' or gene=='None':
            continue
        if first is None:
            vals = [str(avgLR)]
            first = gene
        elif pos==1:  # we're starting a new operon here
            write_cur_vals()
            vals = [str(avgLR)]
            first = gene
        else:
            vals.append(str(avgLR))
    write_cur_vals()