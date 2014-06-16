import math
import pandas.io.sql as psql
import numpy as np
import csv, mysql.connector

con=mysql.connector.connect(user='root', host='127.0.0.1', database='t')
genes=psql.read_frame('select gene_name,gene_strand,gene_posleft,gene_posright,operon_pos,second_best as AVG_LR from GENE',con)
#genes = psql.read_frame("select gene_name,gene_strand,gene_posleft,gene_posright,operon_pos, AVG(aerobic_t8_4_LR+aerobic_t8_5_LR) as avg_lr FROM gene LEFT JOIN aerobic_data ON aerobic_data.RegulonDBID=GENE.gene_id GROUP BY gene_id",con)
print "Big join SQL query done!"
PECgenes=psql.read_frame('select gene_name from PEC',con)
PECgenes=PECgenes['gene_name'].values

tails=[]  # Append pointers to tails here
for strand in ['forward','reverse']:
    last=None
    genesInDir = genes[(genes['gene_strand']==strand)]
    if strand == 'forward':
        genesInDir.sort('gene_posleft',ascending=True,inplace=True)
    else:
        genesInDir.sort('gene_posright',ascending=False,inplace=True)
    for (gene,strand,left,right,pos,avgLR) in genesInDir.itertuples(index=False):
        if gene =='' or gene=='None':
            continue
        if last is None:
            count=1
            last={'gene':gene,'val':avgLR,'prev':None}
        elif pos==1:
            #print "%s has count %i"%(gene,count)
            tails.append(last)
            last={'gene':gene,'val':avgLR,'prev':None}
            count=1
        else:
            count+=1
            cur={'gene':gene,'prev':last,'val':avgLR}
            last=cur
    tails.append(last)


essentials_with_essential_later = 0
first_essentials = 0
genes_before_essential = 0

for cur in tails:
    seen_essential = False
    while cur is not None:
        if seen_essential:
            genes_before_essential += 1
        if cur['gene'] in PECgenes:
            if seen_essential:
                essentials_with_essential_later += 1
            else:
                first_essentials += 1
            seen_essential = True
        cur = cur['prev']

print "# essentials at end of operon",first_essentials,' and # before another essential',essentials_with_essential_later
print 'Total # of genes (essential+nonessential) before an essential gene:',genes_before_essential