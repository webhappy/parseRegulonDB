'''
Created on Dec 16, 2012

@author: cambray

Parse bam file and use MD tag to classify sequences
Keep track of:
- correct read
- mutated sequence count
- others

'''

import pysam
import sys
from parameters import *


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


def check_single_end(read):
    if not "MD" in [t[0] for t in read.tags]:
        return samfile.getrname(read.rname), ("good",)
    else:
        mutseq, deletions = parse_MD_tag(read)
        if not mutseq and not deletions:
            return samfile.getrname(read.rname), ("good",)
        else:
            return samfile.getrname(read.rname), ("mutated", "NA", mutseq, "NA")
            
        
        

          
#     # check quality  
#     average_qual_read = sum([(ord(q) - 33) for q in read.qual])/len(read.qual)
#     average_qual_mate = sum([(ord(q) - 33) for q in mate.qual])/len(mate.qual)

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
    mapping = {}
    h = open(mapping_file)
    h.readline()
    for l in h:
        data = l.strip().split(",")
        data[0] = int(data[0])
        data[1] = int(data[1])
        data[2] = str(data[2])
        mapping[(data[0],data[1])] = data[2]
    h.close()
    
        
    # parse reference database
    db_fasta = open(db_fasta)
    refs = {}
    for l in db_fasta:
        if ">" in l:
            name = l[1:].strip()
        if not ">" in l:
            refs[name] = l.strip()
            #print name, refs[name]
    db_fasta.close()
    
    
    # initialize output data structure
    # indexed by reference, then by conditions, then file id, then by tag
    hits = {}
    for ref in refs:
        hits[ref]={}
        for map in mapping:
            hits[ref][map] = []

    
    # initialize stats
    total     = {}
    good      = {}
    mut       = {}
    mut_total = {}
    mutated   = {}
    unmap     = {}
    
    
    # MAIN analysis
    
    for map in mapping:
        
        total[map]     = 0.
        good[map]      = 0.
        mut[map]       = 0.
        mut_total[map] = 0.
        mutated[map]   = {}
        unmap[map]     = {}
        
        print "== barcode #%s ==" % mapping[map]
        mode = output_mode
        
        for id in file_ids:
            
            print "== file #%i ==" % id
            
            samfile = pysam.Samfile("%s/%s_%03d_dmpx_sorted.bam" % (output_folder, mapping[map], id), "rb" )     
            
            # iterate over reads
            for read in samfile:
                
                total[map] += 1
                ref = read.qname
                
                # update stdout stats
                if total[map]%10000 == 0:
                    print "%i - good: %.2f%% ; mut: %.2f%%, mut id: %i ; unmapped: %.2f%%" % (total[map], good[map]/float(total[map])*100, mut_total[map]/float(total[map])*100, mut[map], len(unmap[map])/float(total[map])*100)
                
                if read.is_unmapped:
                    ref, decision = 0, (0,)
                    unmap[map][ref] = read
                else:
                    ref, decision = check_single_end(read)
                
                
                # check status if single read or paired  
                if decision[0] == "good":
                    good[map]+=1
                    hits[ref][map].append(ref)
                elif decision[0] == "mutated":
                    mut_total[map]+=1
                    if mutated[map].has_key(decision[2]):
                        mutated[map][decision[2]]["cover"]+=1
                    else:
                        mutated[map][decision[2]]={"cover":1,"ref":ref}
                        mut[map]+=1
            
            
            samfile.close()
            
            # print log, unmapped and mutants outputs separately 
            log = open("%s/mapping_%s.log" % (output_folder, mapping[map]), mode)

            non_zero = [len(hits[ref][map]) for ref in refs if len(hits[ref][map])!=0]
            
            to_print  = "Barcode %s %s\n#design hit: %i ; average depth: %.2f; median depth: %i\n" % (mapping[map],map,len(non_zero),float(sum(non_zero))/len(non_zero), non_zero[len(non_zero)/2])
            
            to_print +="""
            total read             : %i
            good                   : %i (%.2f %%)
            mutant #               : %i
            mutated                : %i (%.2f %%)
            unmapped               : %i""" % (total[map], good[map], good[map]/total[map]*100, mut[map], mut_total[map], mut_total[map]/total[map]*100, 2*len(unmap[map]))
           
            print "\n----------\n"+to_print+"\n----------\n"
            log.write(to_print+"\n\n")
            log.close()
     
            unmapped_file = open("%s/unmapped_%s.log" % (output_folder, mapping[map]), mode)
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
    
            mutant_file =open("%s/matant_%s.log" % (output_folder, mapping[map]), mode)
            mutant_file.write("coverage,ref,seq\n")
            for m in mutated[map]:
                mutant_file.write("%i,%s,%s\n" % (mutated[map][m]["cover"], mutated[map][m]["ref"], m))
            mutant_file.close()
            
            mode = "a"
        
    # Compile coverage and freq data per conditions
    h_mapped = open("%s/mapping_summary_%s.csv" % (output_folder, "-".join([str(id) for id in file_ids])), mode)
    h_mapped.write("ref,count_%s,freq_%s\n" % (",count_".join([str(mapping[map]) for map in mapping]),",freq_".join([str(mapping[map]) for map in mapping])))
    total_reads = {}
    for map in mapping:
        total_reads[map] = 0.
        for ref in refs:
            total_reads[map]+=len(hits[ref][map])
    for ref in refs:
        h_mapped.write("%s,%s,%s\n" % (ref,
                                       ",".join([str(len(hits[ref][map])) for map in mapping]),
                                       ",".join([str(len(hits[ref][map])*(10**6)/total_reads[map]) for map in mapping])))
    h_mapped.close()

