# Count how often I get a different result if I opt to sort by p-value versus LR

import mysql.connector
from GenericUtils import *
import csv

DIST_FROM_ATG_CUTOFF = 200  # Can't just use the gene label from All data since some sgRNAs refer to the gene with a different name

outfile=csv.writer(open('bad_sgRNAs.csv','wb'))
outfile.writerow(['Gene name','# sgRNA','Mean LR','Median LR','sgRNA LR','sgRNA p-val','position', 'sgRNA counts sum'])

genes_with_problem = []  # append names of genes to here
cnx = mysql.connector.connect(user='root', host='127.0.0.1',  database='t')
cursor = cnx.cursor(buffered=True)
cursor.execute('select gene_id, gene_name, gene_posleft, gene_posright, gene_strand from gene')
for (gene_id,name,left,right,strand) in cursor:
    gene_has_problem = False
    cursorT=cnx.cursor(buffered=True)
    if strand == 'forward':
        # it's on the positive orientation, so want to let right-side by smaller of right and left+100
        right = min(right,DIST_FROM_ATG_CUTOFF + left)
    else:
        left = max(left,right - DIST_FROM_ATG_CUTOFF)

    flipped=convert_strand_to_minus_plus(flipStrand(strand))
    cursorT.execute("select s.pos, s.seq,p.pval,p.logFC, p.countsSum from sgRNA s LEFT JOIN aerobic_PVALS p on s.seq=p.seq"+
                    " where strand=%s and s.pos<=%s and s.pos>=%s",(flipped,right,left))

    items = []
    for (pos,seq,pval,logFC, countsSum) in cursorT:
        try:
            pval = float(pval)
            logFC = float(logFC)
            #dist=getDistanceForSgRNAFromStart(pos,strand,left,right)
            items.append((pos,seq,countsSum,pval,logFC))
        except:
            print seq,'has no pVal'

    if len(items) == 0:
        print name,'has no sgRNAs'
        continue

    if len(items) <= 2:
        print name,'has only',len(items),'sgRNAs, skipping'
        continue

    sgRNA_logFC = [ x[-1] for x in items ]  # extract only logFC
    gene_mean = meanOfList(sgRNA_logFC)
    gene_median = medianOfList(sgRNA_logFC)

    sorted_by_p_val = sorted(items, key=lambda x: x[-2])  # sorted by pvalue
    sorted_by_logFC = sorted(items, key=lambda x:x[-1])  # sorted by logFC
    if sorted_by_logFC[1] != sorted_by_p_val[1]:
        print name,'is tricky!'

    if gene_has_problem:
        genes_with_problem.append(name)

print len(genes_with_problem),'genes have at least one problematic sgRNA!'