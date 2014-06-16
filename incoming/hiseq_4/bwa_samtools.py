'''
Created on Jan 14, 2013

@author: cambray
'''


from subprocess import call
import sys
from parameters import *
from os import system


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
    
    # main analysis
    
    for id in file_ids:
        
        print "== file #%i ==" % id
        
        for map in mapping:
            
            # I/O
            print "= barcode %s =" % mapping[map]
            
            fastq_R1 = "%s/R1_%s_%03d_dmpx.fastq" % (output_folder, mapping[map], id)
            sai_R1   = "%s/R1_%s_%03d_dmpx.sai" % (output_folder, mapping[map], id)
            
            if paired:
                fastq_R2 = "%s/R2_%s_%03d_dmpx.fastq" % (output_folder, mapping[map], id)
                sai_R2   = "%s/R2_%s_%03d_dmpx.sai" % (output_folder, mapping[map], id)
            
            sam_out_1 = "%s/%s_%03d_dmpx_1.sam" % (output_folder, mapping[map], id)
            bam_out_1 = "%s/%s_%03d_dmpx_sorted_1" % (output_folder, mapping[map], id)

            if paired:
                sam_out_2 = "%s/%s_%03d_dmpx_2.sam" % (output_folder, mapping[map], id)
                bam_out_2 = "%s/%s_%03d_dmpx_sorted_2" % (output_folder, mapping[map], id)

            print "Analyzing:\n%s\n%s\n" % (fastq_R1, fastq_R2)

            ## use bwa for aligning and mapping of both read of a pair

            print "Calling bwa"
            
            if make_index:
                system('./bwa index %s' % db_fasta)
            make_index=False
    
            #align with BWA 
            call(['./bwa aln "%s" "%s" > "%s"' % (db_fasta, fastq_R1, sai_R1)], shell=True)
            if paired:
                call(['./bwa aln "%s" "%s" > "%s"' % (db_fasta, fastq_R2, sai_R2)], shell=True)
            
            #finalize mapping
            #if not paired:
            #    call('./bwa samse "%s" "%s" "%s" > "%s"' % (db_fasta, sai_R1, fastq_R1, sam_out), shell=True)
            #else:
                #paired-end mapping: A=does not use insert size estimate for SW alignment; a=size of insert to use instead
            #call('./bwa sampe -A -a 150 "%s" "%s" "%s" "%s" "%s" > "%s"' % (db_fasta, sai_R1, sai_R2, fastq_R1, fastq_R2, sam_out), shell=True)
            call('./bwa samse "%s" "%s" "%s" > "%s"' % (db_fasta, sai_R1, fastq_R1, sam_out_1), shell=True)
            if paired:
                call('./bwa samse "%s" "%s" "%s" > "%s"' % (db_fasta, sai_R2, fastq_R2, sam_out_2), shell=True)

            # get rid of intermediate sai files
            print "Cleaning"
            call('rm "%s"' % sai_R1, shell=True)
            if paired:
                call('rm "%s"' % sai_R2, shell=True)
            
            # convert sam to bam and sort using samtools
            print "Calling SAMtools"
            #convert sam to bam, and pipe to sort (b: output to BAM, t: tab delimited sam input, u: uncompressed BAM for piping)
            #call('./samtools view -ut "%s" "%s" | samtools sort -n - "%s"' % (db_fasta, sam_out, bam_out), shell=True)
            call('./samtools view -ut "%s" "%s" > "%s"' % (db_fasta, sam_out_1, bam_out_1), shell=True)
            if paired:
                call('./samtools view -ut "%s" "%s" > "%s"' % (db_fasta, sam_out_2, bam_out_2), shell=True)
            
            # get rid of intermediate sam
            print "cleaning"
            call('rm "%s"' % sam_out_1, shell=True)
            if paired:
                call('rm "%s"' % sam_out_2, shell=True)
            print