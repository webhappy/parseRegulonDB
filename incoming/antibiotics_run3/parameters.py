'''
Created on Aug 18, 2013

@author: cambray
'''

### PATHS
output_folder = "output_with_"  # No trailing slash
mapping_file  = "mapping.csv"
barcodes_file = "barcodes.csv"
fastq_R1      = "OLS-antibiotics-3-test_S1_L001_R1_001.fastq/OLS-antibiotics-3-test_S1_L001_R1"
fastq_R2      = ""
file_ids      = [1]
db_fasta      = "../all_oligos.fas"
water_path    = "/usr/local/bin/water"
output_mode   = "w"


### POSITIONS
bc_offset     = 0     # account for discrepancy between bc id and number of N
start_bc      = 0     # which sequence position to start looking for bc
stop_bc       = 16    # which sequence position to stop looking for (e.g. maximum spacer + length bc), only used for water
ROI_START_R1  = 23    # distance from barcode to start of region of interest, it's gcaactctctactgtttctccat for us
ROI_LENGTH    = 20    # length of the region of interest
ROI_OFFSET    =  5    # extra nt to include to account for potential indels (adds this amount on both sides)


### OTHERS
aln_threshold = 0.8
make_index    = True
verbose       = True

def find_match (seq, table, possible_starts):
    for k in possible_starts:
        if seq[k:k+20] in table:
            return seq[k:k+20]
    return None