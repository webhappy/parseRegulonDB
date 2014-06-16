'''
Created on Dec 15, 2012

@author: cambray

Find barcodes/barcode in fastq file from HT sequencing
Try to find exact match where expected
Otherwise, use local alignment to identify barcodes
'''

#TODO: maybe write to file when no match and use waterall instead of water to improve efficiency


from subprocess import call
from Bio import SeqIO
import re
#from basic_seq_manip import revcomp


### PARAMETERS


output_folder = "/Users/Toro/Documents/Dropbox/Data/Cas9/Miseq/2014-02-04 OLS in LB 4 and 5/output/"
mapping_file  = "/Users/Toro/Documents/Dropbox/Data/Cas9/Miseq/2014-02-04 OLS in LB 4 and 5/mapping.csv"
barcodes_file = "/Users/Toro/Documents/Dropbox/Data/Cas9/Miseq/2014-02-04 OLS in LB 4 and 5/barcodes.csv"
fastq_R1      = "/Users/Toro/Documents/Dropbox/Data/Cas9/Miseq/2014-02-04 OLS in LB 4 and 5/OLS-in-LB-4-5_S1_L001_R1_001.fastq"
fastq_R2      = ""
bc_offset     = 0
water_path    = "/usr/local/bin/water" 
aln_threshold = 0.8
start_bc      = 0     # which sequence position to start looking for barcode
stop_bc       = 16    # which sequence position to stop looking for barcode (e.g. maximim spacer + length bc)
ROI_start_R1  = 23    # distance from barcode to start of region of interest
ROI_length    = 20    # length of the region of interest
ROI_start_R2  =  0
ROI_offset    =  3    # extra nt to include to account for potential indels
verbose       = True

#### FUNCTIONS


def find_perfect_barcode(seq, orientation="r1", offset=-1):
    """
    find exact match at exact proper position
    position can be adjusted by offset
    return barcode_id, stop position and alignment score (-1s if no match)
    """
    for bc in barcodes[orientation]:
        start = barcodes[orientation][bc]-offset
        stop  = barcodes[orientation][bc]+len(bc)-offset
        if bc == str(seq[start:stop]):
            return barcodes[orientation][bc], stop, 1
    return -1, -1, -1


def find_degenerate_barcode(seq, orientation="r1", offset=-1):
    """
    find degenerate match first where expected
    Then by scanning the region (expect all barcodes to be of the same size)
    """ 

    # where expected
    for bc in degenerate_barcodes[orientation]:
        bc_len = len(bc)
        start = degenerate_barcodes[orientation][bc]-offset
        stop  = degenerate_barcodes[orientation][bc]+len(bc)-offset
        if bc == str(seq[start:stop]):
            return degenerate_barcodes[orientation][bc], stop, (bc_len-1)/float(bc_len)
    
    #scan sequence 
    matches = []
    for i in range(len(seq)-bc_len):
        subseq = seq[i:i+bc_len].tostring()
        if subseq in degenerate_barcodes[orientation]:
            matches.append((degenerate_barcodes[orientation][subseq],i+bc_len))
    
    # return if only one match with arbitrary score
    if len(matches)==1:
        return matches[0][0], matches[0][1]+8, 0.5 
    
    return -1, -1, -1



def water_bc(seq, orientation, threshold=0.8):
    """
    Align barcodes to sequence and return higest high above threshold
    """
    match = {}
    # write sequence to temp file
    h=open("%s/temp" % output_folder, "w")
    h.write(str(seq))
    h.close()
    
    # call alignment and read output file
    call('%s "%s/temp" "%s/bc_%s.fas" -gapopen 0.5 -gapextend 5 "%s/temp.aln" -sformat2 fasta -aformat fasta -awidth3 200 -auto -snucleotide1 -snucleotide2' % (water_path, output_folder, output_folder, orientation, output_folder), shell=True)
    h=open("%s/temp.aln" % output_folder, "r")
    lines = h.readlines()
    h.close()
    
    #parse alignments and get scores 
    for i in range(0,len(lines),4):
        if len(lines[i+1]) >= bc_len:
            score = 0.
            for j in range(len(lines[i+1])):
                if lines[i+1][j] == lines[i+3][j]:
                    score+=1
            id = lines[i+2].strip()[1:]
            score /= len(lines[i+1])
            match[score] = (int(lines[i+2][1]), lines[i+1], lines[i+3], score)
    hi_score = max(match.keys())  
    motif = match[hi_score][1].strip().replace("-","")
    pos   = seq.find(motif)+len(motif)
    if hi_score >= threshold:
        return match[hi_score][0], pos, hi_score
    else:
        return -1, pos, hi_score



if __name__ == "__main__":


    ##### INITIALIZE OUTPUTS AND DATA STRUCTURE
    
    # Determine if paired-end (ie is a R2 file provided)
    paired = True
    if fastq_R2 == "":
        paired = False
    
    # get mapping
    mapping = {}
    ids     = {"r1":[],"r2":[]}
    h = open(mapping_file)
    h.readline()
    for l in h:
        data = l.strip().split(",")
        data[0] = int(data[0])
        data[1] = int(data[1])
        data[2] = str(data[2])
        mapping[(data[0],data[1])] = data[2]
        if not data[0] in ids["r1"]:
            ids["r1"].append(data[0])
        if paired and not data[1] in ids["r2"]:
            ids["r2"].append(data[1])
    h.close()
    
    # get necessary barcodes according to mapping
    barcodes = {"r1":{},"r2":{}}
    h = open(barcodes_file)
    h.readline()
    len_bc = 0
    for l in h:
        data = l.strip().split(",")
        data[0] = str(data[0]).lower()
        data[1] = str(data[1]).upper()
        data[2] = int(data[2])
        if data[2] in ids[data[0]]:
            barcodes[data[0]][data[1]]=data[2]
            bc_len = len(data[1]) 
    h.close()
    
    # write barcodes to fasta files (used for alignment)
    for dir in ("r1","r2"):
        h = open("%s/bc_%s.fas" % (output_folder, dir) ,"w")
        for seq in barcodes[dir]:
            h.write(">%s\n%s\n" % (barcodes[dir][seq],seq))
        h.close()
    
    
    # Produce degenerate barcodes to speed up matching via hash table
    degenerate_barcodes={"r1":{},"r2":{}}
    for dir in barcodes:
        for seq in barcodes[dir]:                   
            for i in range(len(seq)):
                for base in ("A","T","G","C","N"):
                    if base != seq[i]:
                        split_seq = list(seq)
                        split_seq[i]=base
                        degenerate_barcodes[dir]["".join(split_seq)]=barcodes[dir][seq] 
    
    # fastq outputs
    h_dmpx_R1 = {}
    h_dmpx_R2 = {}
    for map in mapping:
        h_dmpx_R1[map] = open(output_folder+"/R1_"+mapping[map]+"_dmpx.fastq","w")
        if paired:
            h_dmpx_R2[map] = open(output_folder+"/R2_"+mapping[map]+"_dmpx.fastq","w")
    
    # log
    h_log = open(output_folder+"/dmpx.log", "w")
    
    #########   MAIN ANALYSIS
    
    records_R1 = SeqIO.parse(fastq_R1, "fastq")
    if paired:
        records_R2 = SeqIO.parse(fastq_R2, "fastq")
     
    nbr   = 0
    bcmm  = 0
    if paired:
        nobc1  = 0
        nobc2  = 0
        nobc12 = 0
    else:
        nobc   = 0 
    
    
    readperbc = {}
    for map in mapping:
        readperbc[mapping[map]] = 0
    
    for record_R1 in records_R1:
        
        nbr+=1
        
        if paired:
            record_R2 = records_R2.next()
        
        if verbose and nbr % 10000 == 0:
            if paired:
                print """
    # %i
    no bc R1   = %.2f%%
    no bc R2   = %.2f%%
    no bc      = %.2f%%
    wrong bc   = %.2f%%
    %% per bc  = %s
    """ % (nbr, nobc1/float(nbr)*100, nobc2/float(nbr)*100, nobc12/float(nbr)*100, bcmm/float(nbr)*100, " | ".join(["%s: %.2f" % (key, readperbc[key]/float(nbr)*100) for key in readperbc if readperbc[key] != 0 ]))
            
            else:
                print """
    # %i
    no bc      = %.2f%%
    wrong bc   = %.2f%%
    %% per bc  = %s
    """ % (nbr, nobc/float(nbr)*100, bcmm/float(nbr)*100, " | ".join(["%s: %.2f" % (key, readperbc[key]/float(nbr)*100) for key in readperbc if readperbc[key] != 0 ]))
            

            
            
        # Assess quality of the reads
    
    #     average_q_R1 = sum(record_R1.letter_annotations["phred_quality"])/float(len(record_R1))
    #     average_q_R2 = sum(record_R1.letter_annotations["phred_quality"])/float(len(record_R2))
    #     if average_q_R1 <= 30:
    #         print average_q_R1
    #     if average_q_R2 <= 30:
    #         print average_q_R2
    

        
        # first try exact barcode match
        bc_R1, pos_R1, score_R1 = find_perfect_barcode(record_R1.seq[start_bc:stop_bc], "r1", bc_offset)
        # if failed, search 1 bp degenerate via hash table (extend the region searched to allow for 1 insertion)
        if bc_R1 == -1:
            bc_R1, pos_R1, score_R1 = find_degenerate_barcode(record_R1.seq[start_bc:stop_bc+1], "r1", bc_offset)
        # if still fail look for pairwise alignment (extend the region searched to allow for 2 insertion)
        #if bc_R1 == -1:
            #bc_R1, pos_R1, score_R1 = water_bc(record_R1.seq[start_bc:stop_bc+2], "r1", aln_threshold)
        
        # Check final results if only single end data 
        if not paired:
            map = (bc_R1, 0)
            if bc_R1 == -1:
                h_log.write("no bc,%s\n" % (record_R1.id))
                nobc += 1
                continue
            elif not map in mapping:
                bcmm+=1
                h_log.write("unexpected bc,%s,%s\n" % (bc_R1, record_R1.id))
                continue
        
        
        # Similarly check for downstream barcodes if paired end
        else:
            # first try exact barcode match
            bc_R2, pos_R2, score_R2 = find_perfect_barcode(record_R2.seq[start_bc:stop_bc], "r2", bc_offset)
            # if failed, search 1 bp degenerate via hash table (extend the region searched to allow for 1 insertion)
            if bc_R2 == -1:
                bc_R2, pos_R2, score_R2 = find_degenerate_barcode(record_R2.seq[start_bc:stop_bc+1], "r2", bc_offset)
            # if still fail look for pairwise alignment (extend the region searched to allow for 2 insertion)
            #if bc_R2 == -1:
                #bc_R2, pos_R2, score_R2 = water_bc(record_R2.seq[start_bc:stop_bc+2], "r2", aln_threshold)

            # check bc consistency
            map = (bc_R1, bc_R2)
            if bc_R1 == bc_R2 == -1:
                h_log.write("no bc,%s\n" % (record_R1.id))
                nobc12 += 1
                continue
            elif bc_R1 == -1:
                #print "no bc up\n%s\n" % (record_R1.format("fastq"))
                h_log.write("no R1 bc,%s\n" % (record_R1.id))
                nobc1 += 1
                continue
            elif bc_R2 == -1:
                #print "no bc up\n%s\n" % (record_R1.format("fastq"))
                h_log.write("no R2 bc,%s\n" % (record_R1.id))
                nobc2 += 1
                continue
            elif not map in mapping:
                bcmm+=1
                h_log.write("unexpected bc,%s,%s,%s,%s\n" % (bc_R1, bc_R2, record_R1.id, record_R2.id))
            continue
        
        
        #### WRITE UPDATED FASTQ
        # only executed if passed check above above
        # trim bc and most of priming sequence
        # leave few nucleotides in case of indel in the priming region (introduced by primer)

        SeqIO.write(record_R1[pos_R1+ROI_start_R1-ROI_offset:pos_R1+ROI_start_R1+ROI_length+ROI_offset], h_dmpx_R1[map], "fastq")
        if paired:
            SeqIO.write(record_R2[pos_R2+ROI_start_R2-ROI_offset:pos_R2+ROI_start_R2+ROI_length+ROI_offset], h_dmpx_R2[map], "fastq")
        
        readperbc[mapping[map]]+=1
    
    if paired:
        toprint = """
    read pair = %i
    no bc R1  = %s
    no bc R2  = %s
    no bc     = %s
    wrong bc  = %s
    %% per bc = %s
    
    """ % (nbr, nobc1, nobc2, nobc12, bcmm, " | ".join(["%s: %.2f" % (key, readperbc[key]/float(nbr)) for key in readperbc if readperbc[key] != 0 ]))
         
    else:
        toprint="""
    read #    = %i
    no bc     = %s
    wrong bc  = %s
    %% per bc = %s
    
    """ % (nbr, nobc, bcmm, " | ".join(["%s: %.2f" % (key, readperbc[key]/float(nbr)) for key in readperbc if readperbc[key] != 0 ]))
    h_log.write(toprint)
    print toprint