'''
Created on Sep 19, 2013

@author: cambray
'''

### PATHS
output_folder = r"/Users/Toro/Documents/Dropbox/Data/Cas9/Miseq/2014-02-04 OLS in LB 4 and 5/output"
mapping_file  = r"/Users/Toro/Documents/Dropbox/Data/Cas9/Miseq/2014-02-04 OLS in LB 4 and 5/mapping.csv"
barcodes_file = r"/Users/Toro/Documents/Dropbox/Data/Cas9/Miseq/2014-02-04 OLS in LB 4 and 5/barcodes.csv"
fastq_R1      = r"/Users/Toro/Documents/Dropbox/Data/Cas9/Miseq/2014-02-04 OLS in LB 4 and 5/OLS-in-LB-4-5_S1_L001_R1_001.fastq"
fastq_R2      = ""
file_ids      = [1]
db_fasta      = r"/Users/Toro/Documents/Dropbox/Data/Cas9/Miseq/OLS Master key - unique oligos 32992 with 5bp context.fas"
water_path    = r"/usr/local/bin/water"
output_mode   = "w"


### POSITIONS
bc_offset     = 0     # account for discrepancy between bc id and number of N
start_bc      = 0     # which sequence position to start looking for bc
stop_bc       = 16    # which sequence position to stop looking for (e.g. maximum spacer + length bc)
ROI_start_R1  = 23    # distance from barcode to start of region of interest
ROI_length    = 20    # length of the region of interest
ROI_start_R2  = 0
ROI_offset    =  3    # extra nt to include to account for potential indels


### OTHERS
aln_threshold = 0.8
make_index    = True
verbose       = True