import mysql.connector
from GenericUtils import *
import csv

DIST_FROM_ATG_CUTOFF = 200  # Can't just use the gene label from All data since some sgRNAs refer to the gene with a different name

outfile_should_be_sick=csv.writer(open('sgRNAs_that_should_be_sick.csv','wb'))
outfile_should_be_sick.writerow(['Gene name','# sgRNA','Mean LR','Median LR','sgRNA LR','sgRNA p-val','position', 'sgRNA counts sum'])
outfile_should_be_healthy=csv.writer(open('sgRNAs_that_should_be_healthy.csv','wb'))
outfile_should_be_healthy.writerow(['Gene name','# sgRNA','Mean LR','Median LR','sgRNA LR','sgRNA p-val','position', 'sgRNA counts sum'])

genes_with_problem = []  # append names of genes to here
sgRNAs_that_should_be_sick = []
sgRNAs_that_should_be_healthy = []
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

    sgRNA_logFC = [ x[-1] for x in items ]  # extract only logFC
    gene_mean = meanOfList(sgRNA_logFC)
    gene_median = medianOfList(sgRNA_logFC)

    for (pos,seq,countsSum,pval,logFC) in items:
        if gene_median < -2:  # each logFC should be sick
            if logFC > gene_median+2 and countsSum > 15:
                gene_has_problem = True
                sgRNAs_that_should_be_sick.append((pos,seq,countsSum,pval,logFC))
                outfile_should_be_sick.writerow([name,len(items),gene_mean,gene_median,logFC,pval,pos, countsSum])
        else:  # this gene is not sick
            if logFC < gene_median-2 and pval < .01:
                gene_has_problem = True
                sgRNAs_that_should_be_healthy.append((pos,seq,countsSum,pval,logFC))
                outfile_should_be_healthy.writerow([name,len(items),gene_mean,gene_median,logFC,pval,pos, countsSum])

    if gene_has_problem:
        genes_with_problem.append(name)

print len(sgRNAs_that_should_be_healthy),'sgRNAs should have been healthy!'
print len(sgRNAs_that_should_be_sick),'sgRNAs should have been sick'
print len(genes_with_problem),'genes have at least one problematic sgRNA!'