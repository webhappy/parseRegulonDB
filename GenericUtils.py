import mysql.connector
import pandas

def getGCPerc (seq):
    seq=seq.lower()
    total=0
    gc=0
    for j in xrange(len(seq)):
        total+=1
        if seq[j]=='g' or seq[j]=='c':
            gc+=1
    return 1.0*gc/total  # Cast to float

def meanOfList(inList):
    sum = 0.0
    for v in inList:
        sum += v
    return sum/len(inList)

def medianOfList(inList):
    sorts = sorted(inList)
    length = len(sorts)
    if not length % 2:
        return (sorts[length / 2] + sorts[length / 2 - 1]) / 2.0
    return sorts[length / 2]

def flipStrand (strand):
    if strand=='forward':
        return 'reverse'
    elif strand=='reverse':
        return 'forward'
    elif strand == '+':
        return '-'
    elif strand == '-':
        return '+'
    raise Exception('strand doesn\'t match!')

def convert_strand_to_forward_reverse (strand):
    if strand.lower() =='forward' or strand=='+':
        return 'forward'
    elif strand.lower()=='reverse' or strand=='-':
        return 'reverse'
    else:
        raise Exception("Error in strand format")

def convert_strand_to_minus_plus (strand):
    if strand.lower() =='forward' or strand=='+':
        return '+'
    elif strand.lower()=='reverse' or strand=='-':
        return '-'
    else:
        raise Exception("Error in strand format")

def getDistanceForSgRNAFromStart (sgRNApos, strand, gene_left, gene_right):
    """

    :param sgRNApos: Corresponds to left of sgRNA
    :param strand: strand for gene (opposite of strand for sgRNA)
    :param gene_left:
    :param gene_right:
    :return:
    """
    if strand=='forward':
        return sgRNApos-gene_left
    else:
        return gene_right-(20+sgRNApos)


def getSgRNAByPosAndStrand (left, right, strand):
    """Flips the strand orientation inside here before getDistanceFromStart

    :param left:
    :param right:
    :param strand: The actual strand orientation of the gene
    """
    cnx = mysql.connector.connect(user='root',
                              host='127.0.0.1',
                              database='t')
    cursorT=cnx.cursor(buffered=True)
    flipped=convert_strand_to_minus_plus(flipStrand(strand))
    cursorT.execute("select pos,avg(replicate1+replicate2+replicate3) as mean,seq from EXPT_RESULTS "+
                    "where strand=%s and pos<=%s and pos>=%s group by pos,strand",(flipped,right,left))
    items=[]
    minVal=10^6
    for (pos,val,seq) in cursorT:
        dist=getDistanceForSgRNAFromStart(pos,strand,left,right)
        minVal=min(minVal,val)
        items.append((val, dist,seq))

    return items,minVal

# expData = pandas.read_csv('Anaerobic - filtered over 60.csv')
# expData['mean_lr'] = expData[[elem for elem in expData.columns if elem[-3:]=='_LR']].mean(axis=1)
def getSgRNAByPosAndStrandPandas ( left, right, strand):
    """Flips the strand orientation inside here before getDistanceFromStart

    :param left:
    :param right:
    :param strand: The actual strand orientation of the gene
    """
    flipped=convert_strand_to_minus_plus(flipStrand(strand))
    sgRNAs = expData[(expData['sgRNA_pos']<=right)
                     & (expData['sgRNA_pos']>=left) & (expData['sgRNA_strand']==flipped)][['sgRNA_pos','mean']]
   # cursorT.execute("select pos,avg(log2) from EXPT_RESULTS "+
   #                 "where strand=%s and pos<=%s and pos>=%s GROUP BY pos, strand",(flipped,right,left))
    items=[]
    minVal=10^6
    for (pos,val) in sgRNAs.itertuples(index=False):
        dist=getDistanceForSgRNAFromStart(pos,strand,left,right)
        minVal=min(minVal,val)
        items.append((val, dist))

    return items,minVal

all_sgRNAs = {}
def get_genes_hit_by_sgRNA (sgRNA_seq):
    if len(all_sgRNAs) == 0:
        cnx = mysql.connector.connect(user='root', host='127.0.0.1', database='t')
        cursorT = cnx.cursor(buffered=True)
        cursorT.execute('select pos, strand, seq from sgRNA')
        for (pos,strand,seq) in cursorT:
            pos = int(pos)
            all_sgRNAs[seq] = (pos, strand)  # strand is in -/+ form
        print "Finished initializing hash table of all_sgRNAs"

    sgRNA_seq = sgRNA_seq.lower()
    pos = all_sgRNAs[sgRNA_seq][0]
    strand = convert_strand_to_forward_reverse(flipStrand(all_sgRNAs[sgRNA_seq][1]))  # In format to match to gene orientation
    cursorT = cnx.cursor(buffered=True)
    ret = []
    cursorT.execute('select gene_name from GENE where gene_posleft < %s and gene_posright>%s and gene_strand=%s',(pos,pos+20,strand))
    for (gene_name,) in cursorT:
        ret.append(str(gene_name) )
    return ret



if __name__ == '__main__':
    print get_genes_hit_by_sgRNA('actcatggttttctcctgtc')
    # all_data = pandas.io.parsers.read_csv("All data - 32992 with LRs.csv")
    # def getSgRNAsForGene (gene_name):
    #     subset = all_data[all_data['gene_name'] == gene_name]
    #     return subset
    #
    # cnx = mysql.connector.connect(user='root', host='127.0.0.1',  database='t')
    # cursor = cnx.cursor(buffered=False)
    # cursor.execute('select gene_id, gene_name, gene_posleft, gene_posright, gene_strand from gene')
    # genes = []
    # for (gene_id,name,left,right,strand) in cursor:
    #     genes.append(name)
    #
    # print genes
    # for name in genes:
    #     sgRNA_seqs = getSgRNAsForGene(name)['seq'].values
    #     for s in sgRNA_seqs:
    #         print 'Update',s,'to',name
    #         cursor = cnx.cursor(buffered=False)
    #         cursor.execute('update sgRNA set gene=%s where seq=%s',(name,s))
    #     cnx.commit()
