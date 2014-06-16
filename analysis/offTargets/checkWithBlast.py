from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
from Bio import SeqIO
from subprocess import call
from Bio.Seq import Seq
import pandas
import cPickle


refFasta = SeqIO.read("../../sequence.fasta", 'fasta')
refFasta='_'+str(refFasta.seq).upper() # Everything coming out of BLAST is 1-indexed
print len(refFasta)

def findOffTargets (refSeq,sgRNAseq):
    candidates=[]  # Return this list of candidates
    f = open('temp.fasta','wb')
    f.write(sgRNAseq+'\n')
    f.close()

    cline = NcbiblastnCommandline(query="temp.fasta", db="testdb",outfmt=5, out="temp.xml",task='blastn-short')
    cline()
    result=open('temp.xml','r')
    records = NCBIXML.read(result)
    if len(records.alignments) == 0 :
        return candidates
    records=records.alignments[0].hsps

    for record in records:
        if record.identities < 17 or record.query_end < 20:  # Require at least 14 direct matches and ends at the seed
            continue
        if record.sbjct_end > record.sbjct_start:
            end=record.sbjct_end
            # on the + strand, sequence is from [start,end]
            if refSeq[end+2:end+4]=='GG':
                candidates.append(record)
        else:  # On the - strand
            end=record.sbjct_end
            if refSeq[end-3:end-1] == 'CC':
                candidates.append(record)
    return candidates

results={}
outFile=open('results.txt','wb')
expData=pandas.io.parsers.read_csv(open('../../anaerobic.csv','rb'))
i=0
for seq in expData['seq']:
    c=findOffTargets(refFasta,seq)
    if len(c) > 1:
        results[seq]=c
        print seq,'has multiple targets',len(c)
        outFile.write(seq+'\n')
        for record in c:
            if record.sbjct_end > record.sbjct_start:
                querySeq=record.query
                seq='   '+refFasta[record.sbjct_start:record.sbjct_end+4] # positive strand
            else:  # negative strand
                querySeq=Seq(record.query).reverse_complement().tostring()
                seq=refFasta[min(record.sbjct_start,record.sbjct_end)-3:record.sbjct_start+1]+'   '
                #need to reverseComplement

            outFile.write("\t"+"   "+querySeq+"   \t"+str(record.identities)+'\n')
            outFile.write("\t"+"   "+record.match+"   "+'\n')
            outFile.write("\t"+seq+"\t"+str(record.sbjct_end)+'\n')
            outFile.write('\n')
        outFile.write('-----------------\n')
    i+=1
    if i % 100 == 0:
        outFile.flush()
        print i
outFile.close()
cPickle.dump(results,open('pickled.txt','wb'))