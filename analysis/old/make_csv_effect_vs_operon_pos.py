import pandas.io.sql as psql
import numpy as np
import csv, mysql.connector

con=mysql.connector.connect(user='root', host='127.0.0.1', database='t')
genes=psql.read_frame('select gene_name,gene_strand,gene_posleft,gene_posright,operon_pos,median as avg_lr from GENE',con)

outFile = open('genes_vs_operon_pos.csv', 'w')
temp = ['pos'+str(x) for x in range(1,17)]
outFile.write('first_gene,'+','.join(temp))
outFile.write('\n')

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