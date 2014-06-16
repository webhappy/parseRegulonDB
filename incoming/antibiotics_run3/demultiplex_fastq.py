'''
Created on Dec 15, 2012

@author: cambray

Find barcodes/barcode in fastq file from HT sequencing
Try to find exact match where expected
Otherwise, use local alignment to identify barcodes
'''

from subprocess import call
from Bio import SeqIO
import sys
import csv
from parameters import *
from Bio import pairwise2


#### FUNCTIONS


def find_perfect_barcode(seq, orientation="r1", offset=-1):
    """
    find exact match at exact proper position
    position can be adjusted by offset
    return barcode_id, stop position and alignment score (-1s if no match)
    """
    for bc in barcodes[orientation]:
        start = barcodes[orientation][bc]+offset
        stop  = barcodes[orientation][bc]+len(bc)+offset
        if bc == str(seq[start:stop]):
            return barcodes[orientation][bc], stop, 1
    return -1, -1, -1

def find_degenerate_barcode(seq, orientation="r1", offset=-1):
    """
    find degenerate match where expected
    """

    # where expected
    start = offset
    stop  = MAX_BARCODE_LENGTH + offset
    cand = str(seq[start:stop])
    if cand in degenerate_barcodes[orientation] :
        bc_id = degenerate_barcodes[orientation][cand]
        bc_len = MAX_BARCODE_LENGTH - BARCODE_PAD_LENGTH[bc_id]
        return bc_id, start+bc_len, (bc_len-1)/float(bc_len)

    return -1, -1, -1

def water_bc(seq, orientation, id, threshold=0.8):
    """
    Align barcodes to sequence and return higest above threshold
    """
    match = {}
    # write sequence to temp file
    h=open("%s/temp%i" % (output_folder, id), "w")
    h.write(str(seq))
    h.close()
    
    # call alignment and read output file
    water_params = '-sformat2 fasta -aformat fasta -awidth3 200 -auto -snucleotide1 -snucleotide2'
    call('%s "%s/temp%i" "%s/bc_%s.fas" -gapopen 4 -gapextend 10 "%s/temp%i.aln" %s'
         % (water_path, output_folder, id, output_folder, orientation, output_folder, id, water_params), shell=True)
    h=open("%s/temp%i.aln" % (output_folder, id), "r")
    lines = h.readlines()
    h.close()
    
    #parse alignments and get scores 
    for i in range(0,len(lines)-1,4):
        if bc_len-2 <= len(lines[i+1]) <= bc_len+2:
            score = 0.
            for j in range(len(lines[i+1])):
                if lines[i+1][j] == lines[i+3][j]:
                    score+=1
            id = lines[i+2].strip()[1:]
            score /= max(bc_len, len(lines[i+1]))
            match[score] = (int(lines[i+2][1]), lines[i+1], lines[i+3], score)
    if match:
        hi_score = max(match.keys())  
        motif = match[hi_score][1].strip().replace("-","")
        pos   = seq.find(motif)+len(motif)
        if hi_score >= threshold:
            return match[hi_score][0], pos, hi_score
    return -1, -1, -1

def printUnpairedStatus():
    print """
        # %i
        <Q30       = %.2f%%
        no bc      = %.2f%%
        wrong bc   = %.2f%%
        %% per bc  = %s
        Num attempted by water = %i
        Num found by water = %i
        """ % (total_count, badq / float(total_count) * 100, no_barcode_count / float(total_count) * 100, barcode_mismatch_count / float(total_count) * 100, " | ".join(
        ["%s: %.2f" % (key, readperbc[key] / float(total_count) * 100) for key in readperbc if readperbc[key] != 0]),
        num_attempted_by_water, num_found_by_water)

MAX_BARCODE_LENGTH = 0   # Initialize while making the hash table
BARCODE_PAD_LENGTH = {}

def recurse_all_combinations (orig, cur, storeTable, valToAssign, deviations_allowed = 1):
    if len(cur) == len(orig):  # Base this to end recursion
        if cur in storeTable and valToAssign != storeTable[cur]:
            print cur,'is already in the hash table mapping to', storeTable[cur], mapping[(storeTable[cur])]
            raise
        else:
            storeTable[cur] = valToAssign
            return

    LETTERS = ['A', 'T', 'G', 'C', 'N']

    cur_idx = len(cur)
    if orig[cur_idx] == 'N':
        for k in xrange(0,len(LETTERS)):
            recurse_all_combinations(orig, cur+LETTERS[k], storeTable, valToAssign, deviations_allowed)
    else:
        recurse_all_combinations(orig, cur+orig[cur_idx], storeTable, valToAssign, deviations_allowed)

        if deviations_allowed > 0:
            for k in xrange(0,len(LETTERS)):
                if LETTERS[k] == orig[cur_idx]:
                    continue
                recurse_all_combinations(orig, cur+LETTERS[k], storeTable, valToAssign, deviations_allowed-1)

def generateBarcodeHashTable (candidates):
    """Generate dict mapping strings to barcode ID

    :param candidates: Hashmap of barcode string to barcode ID
    :raise:
    """
    global MAX_BARCODE_LENGTH
    ret = {}

    max_length = 0
    for k in candidates.iterkeys():
        if len(k) > max_length:
            max_length = len(k)

    MAX_BARCODE_LENGTH = max_length
    print "Generating barcode hashtable, max_length =",MAX_BARCODE_LENGTH
    modified_candidates = {}
    for k in candidates.iterkeys():
        BARCODE_PAD_LENGTH[candidates[k]] = max_length-len(k)
        to_append = 'N'*(max_length-len(k))
        modified_candidates[k+to_append] = candidates[k]

    for seq in modified_candidates.iterkeys():
        ret[seq] = modified_candidates[seq]
        recurse_all_combinations(seq,'',ret,modified_candidates[seq])
    print "Done making hashtable, contains",len(ret),'items'
    return ret

if __name__ == "__main__":
    
    ##### INITIALIZE OUTPUTS AND DATA STRUCTURE
    
    # Read STDIN (superseed parameter)
    if len(sys.argv)>1:
        file_ids = [int(id) for id in sys.argv[1].split("-")]
    
    # get mapping
    mapping = {}  # Tuple of (id_R1, id_R2) maps to Description (valid filename)
    ids     = {"r1":[],"r2":[]}
    h = open(mapping_file)  # mapping.csv contains id_R1, id_R2, Description
    h.readline()
    for l in h:
        data = l.strip().split(",")
        data[0] = int(data[0])
        data[1] = int(data[1])
        data[2] = str(data[2])
        mapping[(data[0],data[1])] = data[2]
        if not data[0] in ids["r1"]:
            ids["r1"].append(data[0])
    h.close()

    # get necessary barcodes according to mapping
    global bc_len
    barcodes = {"r1":{},"r2":{}}
    barcodeID_to_barcodeString = {}
    h = open(barcodes_file)  # barcodes.csv contains orientation,sequence,id
    h.readline()
    for l in h:
        data = l.strip().split(",")
        data[0] = str(data[0]).lower()
        data[1] = str(data[1]).upper()
        data[2] = int(data[2])
        if data[2] in ids[data[0]]:
            barcodeID_to_barcodeString[data[2]] = data[1]
            barcodes[data[0]][data[1]]=data[2]
            bc_len = len(data[1]) # all barcodes assumed to have same length!
                                  # Notice this assumption in above line!!!
    h.close()
     
    # write barcodes to fasta files (used for alignment)
    for direction in (["r1"]):
        h = open("%s/bc_%s.fas" % (output_folder, direction) ,"w")
        for seq in barcodes[direction]:
            h.write(">%s\n%s\n" % (barcodes[direction][seq],seq))
        h.close()
    
    
    # Produce degenerate barcodes to speed up matching via hash table
    degenerate_barcodes={"r1":{},"r2":{}}
    degenerate_barcodes['r1'] = generateBarcodeHashTable(barcodes['r1'])

    # I/O
    fastq_R1s   = {}
    h_dmpx_R1s  = {}
    h_logs       = {}

    print "The following file(s) will be analyzed:"
    
    for id in file_ids:
        # inputs
        fastq_R1s[id] = "%s_%03d.fastq" % (fastq_R1, id)
        print fastq_R1s[id]

        # logs
        h_logs[id] = open("%s/%03d_dmpx.log" % (output_folder, id), output_mode)
        
        # fastq outputs
        # This is all GC's original method of dumping to new fastq files
        h_dmpx_R1s[id] = {}

        for map in mapping:
            h_dmpx_R1s[id][map] = open("%s/R1_%s_%03d_dmpx.fastq" % (output_folder, mapping[map], id), output_mode)


    # outFiles = {}  # Map from description to csv file handler
    # for m, v in mapping.iteritems():
    #     outFiles[v] = csv.writer(open("output/%s_reads.csv" % v,'wb'))

    #########   MAIN ANALYSIS
    for id in file_ids:
        records_R1 = SeqIO.parse(fastq_R1s[id], "fastq")

        total_count   = 0
        barcode_mismatch_count  = 0
        badq  = 0
        no_barcode_count = 0
        num_attempted_by_water = 0
        num_found_by_water = 0

        readperbc = {}
        for map in mapping:
            readperbc[mapping[map]] = 0

        for record_R1 in records_R1:
            total_count += 1
            if verbose and total_count % 10000 == 0:
                printUnpairedStatus()

            # Assess quality of the reads
            # Only keep paired reads with one mate average score >= Q30        
            flagQ = False
            average_q_R1 = sum(record_R1.letter_annotations["phred_quality"])/float(len(record_R1))
            if average_q_R1 < 30:
                h_logs[id].write("R1Q%.2f,%s\n" % (average_q_R1, record_R1.id))
                flagQ = True
            if flagQ:
                badq += 1
                continue
            
            # first try exact barcode match
            #bc_R1, pos_R1, score_R1 = find_perfect_barcode(record_R1.seq[start_bc:stop_bc], "r1", bc_offset)
            # if failed, search 1 bp degenerate via hash table (extend the region searched to allow for 1 insertion)
            bc_R1, pos_R1, score_R1 = find_degenerate_barcode(record_R1.seq[start_bc:stop_bc+1], "r1", bc_offset)
            # if still fail look for pairwise alignment (extend the region searched to allow for 2 insertion)
            if bc_R1 == -1:
                num_attempted_by_water += 1
                bc_R1, pos_R1, score_R1 = water_bc(record_R1.seq[start_bc:stop_bc+2], "r1", id, aln_threshold)
                #print record_R1.seq[0:30], bc_R1, barcodeID_to_barcodeString[bc_R1]
                #if total_count > 10:
                #    adfdaf()
                if bc_R1 != -1:
                    num_found_by_water += 1
            
            # Check that we can map this read's multiplex barcode
            map_key = (bc_R1, 0)
            if bc_R1 == -1:
                h_logs[id].write("no bc,%s\n" % (record_R1.id))
                no_barcode_count += 1
                continue
            elif not map_key in mapping:
                barcode_mismatch_count += 1
                h_logs[id].write("unexpected bc,%s,%s\n" % (bc_R1, record_R1.id))
                continue

            pos_R2 = 0  # since there is no barcode for read_2
            #### WRITE UPDATED FASTQ
            # only executed if passed check above
            # trim bc and most of priming sequence
            # leave few nucleotides in case of indel in the priming region (introduced by primer)

            #print record_R1.seq,'has barcode',bc_R1
            # r1_match = find_match(str(record_R1.seq).lower(), allOligos, xrange(pos_R1+ROI_start_R1-2,pos_R1+ROI_start_R1+3))
            # r2_match = find_match(str(record_R2.seq).lower(), allOligos, xrange(pos_R2+ROI_start_R2-2,pos_R2+ROI_start_R2+3))
            SeqIO.write(record_R1[pos_R1+ROI_START_R1-ROI_OFFSET:pos_R1+ROI_START_R1+ROI_LENGTH+ROI_OFFSET], h_dmpx_R1s[id][map_key], "fastq")
            readperbc[mapping[map_key]] += 1
            
        toprint="""
        read #    = %i
        no bc     = %s
        wrong bc  = %s
        %% per bc = %s
        
        """ % (total_count, no_barcode_count, barcode_mismatch_count, " | ".join(["%s: %.2f" % (key, readperbc[key]/float(total_count)) for key in readperbc if readperbc[key] != 0 ]))
        h_logs[id].write(toprint)
        print toprint
