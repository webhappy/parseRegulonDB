import os
import subprocess

for k in range(1,5):
    if k == 1:
        continue
    cur_dir = 'hiseq_'+str(k)
    os.chdir(cur_dir)
    print "**********************"
    print "Now in",cur_dir
    subprocess.call('python demultiplex_fastq.py', shell=True)
    print "Done de-multiplexing"
    subprocess.call('python bwa_samtools.py', shell=True)
    print "Done calling bwa_samtools"
    subprocess.call('python parse_bam.py', shell=True)
    print 'Done parse_bam'
    os.chdir('../')

# print '-----Now importing into MySQL DB-----'
# for k in range(1,5):
#     subprocess.call('python importCounts.py', shell=True)
#     print "Imported hiseq file #",k