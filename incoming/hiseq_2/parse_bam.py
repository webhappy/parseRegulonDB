'''
Created on Dec 16, 2012

@author: cambray

Parse bam file and use MD tag to classify sequences
Keep track of:
- correct read
- mutated sequence count
- others

Use new aln (muscle if reference not certain)
Produce a different kind of mutant id that would need fix to fit with the one generated from MD tag

'''

import pysam
from Bio import AlignIO, SeqIO
import subprocess
import sys
from parameters import *
import pandas as pd
import numpy as np


def get_oriented_seq(read1, read2):
    # deal with orientation for unmapped read
    seq1 = read1.seq
    seq2 = None
    if read2:
        seq2 = read2.seq
    if read1.is_read2:
        if not read1.is_reverse:
            seq1 = revcomp(seq1)
    elif read2 and not read2.is_reverse:
        seq2 = revcomp(seq2)
    return seq1, seq2
def check_multimere(read1,read2):
    """
    Check for correct scar on either read (hard coded for now)
    """
    seq1, seq2 = get_oriented_seq(read1, read2)

    if "TATGGTACC" in seq1[0:15] or "GGATCCACC" in seq1[-20:]:
        return True
    if read2:
        if "TATGGTACC" in seq2[0:15] or "GGATCCACC" in seq2[-20:]:
            return True
    return False
def refine_aln1(ref, read1):
    seq1, seq2 = get_oriented_seq(read1, None)
    # use muscle to get alignment
    cline = "muscle -gapopen -1000 -weight1 none -weight2 none -maxiters 1 -diags"
    child = subprocess.Popen(cline,
                             stdin=subprocess.PIPE,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE,
                             shell=(sys.platform!="win32"))
    child.stdin.write(">ref\n%s\n>1\n%s" % (ref, seq1))
    child.stdin.close()
    align = dict([(a.name, a.seq) for a in AlignIO.read(child.stdout, "fasta")])
    #calculate score wrt reference as percent identity
    score   = [0., 0.]
    i=0
    while align["ref"][i]=="-":
        i+=1
    start = i
    i=0
    while align["ref"][-i]=="-":
        i+=1
    stop = len(align["ref"])
    if i != 0:
        stop += -i+1
    map_aln2read1 = {}
    i = -1
    for j in range(len(align["1"])):
        if align["1"][j] != "-":
            i+=1
            map_aln2read1[j]=i
        else:
            map_aln2read1[j]="-"
    for i in range(start, stop):
        if align["ref"][i] == align["1"][i]:
            score[0]+=1
    length = len(align["ref"][start:stop])
    score[0]/=length
    # lower acceptable score threshold for low quality
    threshold1 = 0.8
    #if average_quality1 <= 30:
    #    threshold1 = 0.6
    if score[0]<threshold1:
        return("bad aln1")
    else:
        mutations = {}
        j=-1
        for i in range(start,stop):
            j+=1
            if (align["1"][i] != align["ref"][i]):
                if map_aln2read1[i] != "-":
                    mutations[j] = "%s>%s|%i" % (align["ref"][i], align["1"][i], ord(read1.qual[map_aln2read1[i]])-33)
                else:
                    mutations[j] = "%s>%s|-" % (align["ref"][i], align["1"][i])
        if not mutations:
            return ("good",score)
        else:
            return ("mutated",score, align["1"][start:stop].tostring(), len(mutations), mutations.viewitems())
def refine_aln2(ref, read1, read2):
    """
    redo alignment of two reads to ref and check consistant mutations
    only consider if both read aln are above thresholds
    return "good" in no consistent mutations
    return "mutated" with info if consistent mutations
    """
    seq1, seq2 = get_oriented_seq(read1, read2)
    # use muscle to get alignment
    cline = "muscle -gapopen -1000 -weight1 none -weight2 none -maxiters 1 -diags"
    child = subprocess.Popen(cline,
                             stdin=subprocess.PIPE,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE,
                             shell=(sys.platform!="win32"))
    child.stdin.write(">ref\n%s\n>1\n%s\n>2\n%s" % (ref, seq1, seq2))
    child.stdin.close()
    align = dict([(a.name, a.seq) for a in AlignIO.read(child.stdout, "fasta")])
    #calculate score wrt reference as percent identity
    score   = [0., 0.]
    i=0
    while align["ref"][i]=="-":
        i+=1
    start = i
    i=0
    while align["ref"][-i]=="-":
        i+=1
    stop = len(align["ref"])
    if i != 0:
        stop += -i+1
    map_aln2read1 = {}
    map_aln2read2 = {}
    map_aln2ref   = {}
    i1 = i2 = i3 = -1
    for j in range(len(align["1"])):
        if align["1"][j] != "-":
            i1+=1
            map_aln2read1[j]=i1
        else:
            map_aln2read1[j]="%s-"%i1
        if align["2"][j] != "-":
            i2+=1
            map_aln2read2[j]=i2
        else:
            map_aln2read2[j]="%s-"%i2
        if align["ref"][j] != "-":
            i3+=1
            map_aln2ref[j]=i3
        else:
            map_aln2ref[j]="%s-"%i3
    for i in range(start, stop):
        if align["ref"][i] == align["1"][i]:
            score[0]+=1
        if align["ref"][i] == align["2"][i]:
            score[1]+=1
    length = len(align["ref"][start:stop])
    score[0]/=length
    score[1]/=length
    threshold1 = threshold2 = 0.8
    if score[0]<=threshold1 or score[1]<=threshold2:
        # not dealt with
        return("bad aln",score)
    else:
        ismutant = False
        mutant = list(align["ref"][start:stop])
        j=-1
        for i in range(start,stop):
            j+=1
            # check mutations present on both reads
            if align["1"][i] != align["ref"][i] and align["2"][i] != align["ref"][i]:
                if align["1"][i] == align["2"][i]:
                    mutant[j]=align["1"][i].lower()
                    ismutant = True
        mutseq = "".join(mutant)
#         print align["ref"][start:stop]
#         print mutseq
#         print
        if  not ismutant:
            return ("good",score,"", 2)
        else:
            return ("mutated", score, mutseq, 2)
def parse_MD_tag(read):
    """
    MD tag describes alignment of the reference sequence
    return
    """
    tag = read.opt("MD")
    digit_flag  = False
    flag_delet  = 0
    deletions   = []
    mut=""
    if tag.isdigit():
        return False, deletions
    current_pos = 0
    for c in tag:
        if c.isdigit():
            flag_delet = False
            if not digit_flag:
                digit_flag = True
                position = ""
            position += c
        elif c != "^":
            if digit_flag:
                digit_flag = False
                l = len(mut.replace("-",""))
                mut+=read.query[l:l+int(position)]
                current_pos+=int(position)-1
            if flag_delet:
                mut+="-"
                deletions.append(current_pos+flag_delet-1)
                flag_delet+=1
            else:
                current_pos+=1
                mut+=read.query[current_pos].lower()

        else:
            flag_delet = 1
    if digit_flag:
        l = len(mut.replace("-",""))
        mut+=read.query[l:l+int(position)]
    return mut, deletions
def check_mutations(read, mate):
    """
    calculate position of mutation of full ref
    """
    read_mut = ""
    mate_mut = ""
    if "MD" in [t[0]for t in read.tags]:
        read_mut, read_deletions = parse_MD_tag(read)
    if "MD" in [t[0]for t in mate.tags]:
        mate_mut, mate_deletions = parse_MD_tag(mate)
    # no mutations
    if not mate_mut and not read_mut:
        return False, "", samfile.getrname(read.rname), 2
    elif not read_mut:
        return False, "", samfile.getrname(read.rname), 1
    elif not mate_mut:
        return False, "", samfile.getrname(read.rname), 1
    # check mutation if same ref
    if read.rname != mate.rname:
        #print "different hits:", samfile.getrname(read.rname), samfile.getrname(mate.rname)
        return "diff_hit", "", samfile.getrname(read.rname), 0
    # realign
    roffset2 = sum([1 for d in read_deletions])
    moffset1 = sum([1 for d in mate_deletions])
    moffset2 = sum([1 for d in mate_deletions if d>115])
    # check boundaries to extract exact insert
    # that work for most reads, realign the others
    read_bound = "NA"
    mate_bound = "NA"
    for i in (0,-1,1,-2,2):
        up = (12+i,15+i)
        if  ((read_mut[up[0]:up[1]] =="ATG") or (not read_mut[up[0]:up[1]].isupper())) \
        and ((read_mut[up[0]+99:up[1]+99] =="GGT") or (not read_mut[up[0]+99:up[1]+99].isupper())):
            read_bound = (up[0], up[1]+99)
            continue
    for i in (0,-1,1,-2,2):
        up = (15+moffset1+i,18+moffset1+i)
        if  ((mate_mut[up[0]:up[1]] =="ATG") or (not mate_mut[up[0]:up[1]].isupper())) \
        and ((mate_mut[up[0]+99:up[1]+99] =="GGT") or (not mate_mut[up[0]+99:up[1]+99].isupper())):
            mate_bound = (up[0], up[1]+99)
            continue
    if read_bound == "NA" or mate_bound == "NA":
#         print read.query, read.opt("MD"), samfile.getrname(read.rname)
#         print mate.query, mate.opt("MD"), samfile.getrname(mate.rname)
#         print read_mut
#         print mate_mut
#         if read_bound != "NA":
#             print read_mut[read_bound[0]:read_bound[1]]
#         if mate_bound != "NA":
#             print mate_mut[mate_bound[0]:mate_bound[1]]
#         print
        return "same_hit", "", samfile.getrname(read.rname), 0
#     print read_mut[read_bound[0]:read_bound[1]], read.opt("MD"), samfile.getrname(read.rname)
#     print mate_mut[mate_bound[0]:mate_bound[1]], mate.opt("MD"), samfile.getrname(mate.rname)
#     print
    # get consensus mutations
#     print read.query
#     print read_mut, read.opt("MD")
#     print mate_mut
    read_mut = read_mut[read_bound[0]:read_bound[1]]
    mate_mut = mate_mut[mate_bound[0]:mate_bound[1]]
#     print read_mut
#     print mate_mut
#     print
    consensus = ""
    mutated   = False
    if len(read_mut) != len(mate_mut):
        # problem somewhere... hard to figure out, better re-align
        return "mutated_length", "", samfile.getrname(read.rname), 0
    for i in range(len(read_mut)):
        if read_mut[i] == mate_mut[i]:
            consensus += read_mut[i]
            if read_mut[i].islower():
                mutated = True
        elif read_mut[i].islower() and mate_mut[i].islower():
            consensus += "n"
        elif read_mut[i].islower():
            consensus += mate_mut[i]
        else:
            consensus += read_mut[i]
    return mutated, consensus[3:-3], samfile.getrname(read.rname), 2
def check_single_end(read):
    if not "MD" in [t[0] for t in read.tags]:
        return samfile.getrname(read.rname), ("good",)
    else:
        mutseq, deletions = parse_MD_tag(read)
        if not mutseq and not deletions:
            return samfile.getrname(read.rname), ("good",)
        else:
            return samfile.getrname(read.rname), ("mutated", "NA", mutseq, "NA")
def check_mates(read, mate):
    if check_multimere(read, mate):
        return None, ("multi",)

    # get aln results from MD tag
    mutated, mutseq, ref, support = check_mutations(read, mate)
    if mutated == False:
        return ref, ("good",)
    elif mutated == True:
        return ref, ("mutated", "NA", mutseq, support)
    elif mutated == "same_hit" or mutated == "mutated_length":
        return ref, refine_aln2(refs[ref][26:122], read, mate)
    # not multimere but mapped to two different sequences: redo possible alns
    elif mutated == "diff_hit":
        score1 = refine_aln2(refs[ref][26:122], read, mate)
        score2 = refine_aln2(refs[samfile.getrname(mate.rname)][26:122], mate, read)
        hi_score = score1
        if sum(score2[1])>sum(score1[1]):
            # reassign ref
            ref = samfile.getrname(mate.rname)
            hi_score = score2
        return ref, hi_score
    else:
        print "problem", mutated, mutseq, ref, support
        return None, ("?",)

####################
####################

if __name__ == "__main__":
    
    # Read STDIN (superseed parameter)
    if len(sys.argv)>1:
        file_ids = [int(id) for id in sys.argv[1].split("-")]
    
    print "File to analyze : %s" % " ".join([str(id) for id in file_ids])
    print
    
    # Determine if paired-end (ie is a R2 file provided?)
    paired = True
    if fastq_R2 == "":
        paired = False
    
    # get mapping
    mapping = {}  # Tuple of (id_R1, id_R2) maps to Description
    h = open(mapping_file)
    h.readline()
    for l in h:
        data = l.strip().split(",")
        data[0] = int(data[0])
        data[1] = int(data[1])
        data[2] = str(data[2])
        mapping[(data[0],data[1])] = data[2]
    h.close()

    maps = [str(mapping[map]) for map in mapping]
    maps.sort()  # sort descriptions
        
    # parse reference database
    db_fasta = open(db_fasta)  # Assumes structure is > line followed by non-> line
    refs = {}
    for l in db_fasta:
        if ">" in l:
            name = l[1:].strip()
        if not ">" in l:
            refs[name] = l.strip()
            #print name, refs[name]
    db_fasta.close()
    ref_names = refs.keys()
    ref_names.sort()
    
    # initialize output
    h_mapped = open("%s/mapping_summary_%s.csv" % (output_folder, "-".join([str(id) for id in file_ids])), mode="w")
    h_mapped.write("ref,%s\n" % ",".join("count.%s,freqpm.%s" % (map, map) for map in maps))
    for ref in ref_names:
        h_mapped.write("%s\n" % ref)
    h_mapped.close() 
    
    
    # initialize output data structure
    # indexed by reference, then by conditions, then file id, then by tag
    hits = {}
    for ref in ref_names:
        hits[ref]={}
        for map in maps:
            hits[ref][map] = []

    
    # initialize stats
    total     = {}
    good      = {}
    mut       = {}
    mut_total = {}
    mutated   = {}
    unmap     = {}
    if paired:
        diff     = {}
        map1     = {} # to store first in pair of mates
        unmap1   = {}
    
    for map in maps:
        total[map]     = 0.
        good[map]      = 0.
        mut[map]       = 0.
        mut_total[map] = 0.
        mutated[map]   = {}
        unmap[map]     = {}
        if paired:
            diff[map]   = 0.
            map1[map]   = {} # to store first in pair of mates
            unmap1[map] = {}


    # Initialize pandas df
    data = np.zeros((len(refs),len(maps)))
    df = pd.DataFrame(data,ref_names,maps)
    df['seq'] = pd.Series([refs[x][6:26] for x in ref_names],ref_names)
    
    # MAIN analysis
    
    for map in maps:
        
        print "== barcode #%s ==" % map
        mode = output_mode
        files = ["%s/%s_%03d_dmpx_sorted_1" % (output_folder, map, id) ] #, "%s/%s_%03d_dmpx_sorted_2" % (output_folder, map, id)]

        for file_name in files:
            
            print "== file #%i ==" % id
            samfile = pysam.Samfile(file_name, "rb" )
            # iterate over reads
            for read in samfile:
                total[map] += 1
                ref = read.qname
                
                # update stdout stats
                if total[map]%10000 == 0:
                    if not paired:
                        print "%i - good: %.2f%% ; mut: %.2f%%, mut id: %i ; unmapped: %.2f%%" % (total[map], good[map]/float(total[map])*100, mut_total[map]/float(total[map])*100, mut[map], len(unmap[map])/float(total[map])*100)
                    else:
                        print "%i - good: %.2f%% ; mut: %.2f%%, mut id: %i ; diff. hit: %.2f%% ; unmapped: %.2f%%" % (total[map], 2*good[map]/float(total[map])*100, 2*mut_total[map]/float(total[map])*100, mut[map], diff[map]/float(total[map])*200, len(unmap[map])/float(total[map])*100)
    
                # categorize
                # single-end
                if not paired:
                    if read.is_unmapped:
                        ref, decision = 0, (0,)
                        unmap[map][ref] = read
                    else:
                        ref, decision = check_single_end(read)
                
                # paired-end
                else:
                    # TODO
                    pass
                
                # check status if single read or paired  
                if decision[0] == "good":
                    good[map]+=1
                    hits[ref][map].append(ref)
                    df.loc[ref, map] += 1
                elif decision[0] == "multi":
                    diff[map]+=1
                elif decision[0] == "mutated":
                    mut_total[map]+=1
                    if mutated[map].has_key(decision[2]):
                        mutated[map][decision[2]]["cover"]+=1
                    else:
                        mutated[map][decision[2]]={"cover":1,"ref":ref}
                        mut[map]+=1
            
            samfile.close()
            
            # print log, unmapped and mutants outputs separately 
            log = open("%s/mapping_%s.log" % (output_folder, map), mode)

            non_zero = [len(hits[ref][map]) for ref in ref_names if len(hits[ref][map])!=0]
            if len(non_zero) == 0:
                to_print = "Barcode "+str(map)+" has no counts"
            else:
                to_print  = "Barcode %s %s\n#design hit: %i ; average depth: %.2f; median depth: %i\n" % (map,map,len(non_zero),float(sum(non_zero))/len(non_zero), non_zero[len(non_zero)/2])
                if not paired:
                    to_print +="""
                    total read             : %i
                    good                   : %i (%.2f %%)
                    mutant #               : %i
                    mutated                : %i (%.2f %%)
                    unmapped               : %i""" % (total[map], good[map], good[map]/total[map]*100, mut[map], mut_total[map], mut_total[map]/total[map]*100, 2*len(unmap[map]))
                else:
                    to_print +="""
                    total read             : %i
                    good                   : %i (%.2f %%)
                    mutant #               : %i
                    mutated                : %i (%.2f %%)
                    different hit          : %i
                    semi-mapped            : %i
                    semi-unmapped          : %i
                    unmapped               : %i""" % (total[map], good[map]*2, good[map]*2/total[map]*100, mut[map], mut_total[map], mut_total[map]*2/total[map]*100, diff[map]*2, len(map1[map]), len(unmap1[map]), 2*len(unmap[map]))
            
            print "\n----------\n"+to_print+"\n----------\n"
            log.write(to_print+"\n\n")
            log.close()
     
            unmapped_file = open("%s/unmapped_%s.log" % (output_folder, map), mode)
            if not paired:
                unmapped_file.write("read\n")
            else:
                unmapped_file.write("R1,R2\n")
            for ref in unmap[map]:
                if not paired:
                    unmapped_file.write("%s\n" % str(unmap[map][ref]))
                else:
                    unmapped_file.write("%s,%s\n" % (str(unmap[map][ref][0]), str(unmap[map][ref][1])))
            unmapped_file.close()
    
            mutant_file =open("%s/mutant_%s.log" % (output_folder, map), mode)
            mutant_file.write("coverage,ref,seq\n")
            for m in mutated[map]:
                mutant_file.write("%i,%s,%s\n" % (mutated[map][m]["cover"], mutated[map][m]["ref"], m))
            mutant_file.close()
            
            df.to_csv('counts.csv')