from pandas.io.parsers import *
import numpy as np
import mysql.connector
import pandas as pd
from pandas.io.sql import *
import shutil

cnx = mysql.connector.connect(user='davidc', password='mysql_password', host='127.0.0.1',  database='CRISPR')
cursor = cnx.cursor(True)

runIDs = {}
for k in range(1,9):
    runIDs[k] = 25+k

def generate_parameter_string (filename,runID):  # filename must not contain the final _001.fastq part
    return '''
### PATHS
output_folder = "output"  # No trailing slash
mapping_file  = "mapping.csv"
barcodes_file = "barcodes.csv"
fastq_R1      = "%s"
fastq_R2      = ""
file_ids      = [1]
db_fasta      = "../all_oligos.fas"
water_path    = "/usr/local/bin/water"
output_mode   = "w"
RUN_ID = %i


### POSITIONS
bc_offset     = 0     # account for discrepancy between bc id and number of N
start_bc      = 0     # which sequence position to start looking for bc
stop_bc       = 16    # which sequence position to stop looking for (e.g. maximum spacer + length bc), only used for water
ROI_START_R1  = 23    # distance from barcode to start of region of interest, it's gcaactctctactgtttctccat for us
ROI_LENGTH    = 20    # length of the region of interest
ROI_OFFSET    =  0    # extra nt to include to account for potential indels (adds this amount on both sides)


### OTHERS
aln_threshold = 0.8
make_index    = True
verbose       = True
''' % (filename,runID)

files = ['IT001_S1_L001_R1_001.fastq','TN7-1-IT002_S2_L001_R1_001.fastq','TN7-2-IT003_S3_L001_R1_001.fastq','JR1-IT004_S4_L001_R1_001.fastq']

# Make mapping.csv and barcodes.csv for each of the directories
for j in range(1,5):
    shutil.copyfile('06132014/'+files[j-1],'hiseq_'+str(j)+'/'+files[j-1])
    mapping_file = open('hiseq_'+str(j)+'/mapping.csv','w')
    barcodes_file = open('hiseq_'+str(j)+'/barcodes.csv','w')
    barcodes_file.write('orientation, sequence, id\n')
    mapping_file.write('id_R1, id_R2, Description\n')
    cursor.execute('Select sampleID, barcode, description from samples where runID=%s'%runIDs[j])
    for (sampleID, barcode, description) in cursor:
        mapping_file.write(str(sampleID)+',0,'+description + '\n')
        barcodes_file.write('r1,'+barcode+','+str(sampleID)+'\n')
    parameters_file =open('hiseq_'+str(j)+'/parameters.py','w')
    parameters_file.write(generate_parameter_string(files[j-1][0:-10],runIDs[j]))

    # All done
    parameters_file.close()
    mapping_file.close()
    barcodes_file.close()
