import mysql.connector
from GenericUtils import *
import csv

DIST_FROM_ATG_CUTOFF = 200  # Can't just use the gene label from All data since some sgRNAs refer to the gene with a different name

outfile=csv.writer(open('sgRNAs_per_gene.csv','wb'))
outfile.writerow(['Gene name','# sgRNA','Mean LR','Median LR'])

genes_with_problem = []  # append names of genes to here
cnx = mysql.connector.connect(user='root', host='127.0.0.1',  database='t')
cursor = cnx.cursor(buffered=True)
cursor.execute('select gene_id, gene_name, gene_posleft, gene_posright, gene_strand from gene')
for (gene_id,name,left,right,strand) in cursor:
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

    sgRNA_logFC = [ x[-1] for x in items ]  # extract only logFC
    gene_mean = meanOfList(sgRNA_logFC)
    gene_median = medianOfList(sgRNA_logFC)

    outfile.writerow([name,len(items),gene_mean,gene_median])
