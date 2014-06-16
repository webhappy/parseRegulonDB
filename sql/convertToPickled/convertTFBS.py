import mysql.connector
from Bio import SeqIO
import pickle

refFasta = SeqIO.read("sequence.fasta", 'fasta')

cnx = mysql.connector.connect(user='root',
                              host='127.0.0.1',
                              database='t')
PROM_BOUNDARY = 35;
cursor = cnx.cursor()
cursor.execute('select SITE.site_id,site_sequence,site_posleft,site_posright,ri_function,CONFORMATION.transcription_factor_id, transcription_factor_name from SITE left join `REGULATORY_INTERACTION` ON SITE.`site_id`=`REGULATORY_INTERACTION`.`site_id` left join CONFORMATION on `REGULATORY_INTERACTION`.`conformation_id`=`CONFORMATION`.`conformation_id` left join TRANSCRIPTION_FACTOR ON CONFORMATION.`transcription_factor_id`=`TRANSCRIPTION_FACTOR`.`transcription_factor_id` where site_posleft>0 ORDER BY site_id asc')
siteContents={};
for (id,seq,left,right,function,tfID,tfName) in cursor: #actually don't care about function for now
    left=int(left)
    right=int(right)
    left_boundary = left
    right_boundary=right

    assert (len(refFasta) >= left_boundary >= 0 and 0 <= right_boundary <= len(refFasta))
    siteSeq=''.join([c for c in seq if c.isupper()])
    assert(str(refFasta.seq[left-1:right])==siteSeq.encode('latin-1')) #  Since refFasta.seq is 0-indexed and left and right are 1-indexed
    if id in siteContents:
        #append
        siteContents[id]['TFs'].append(tfName)
    else:
        siteContents[id]={'TFs':[tfName],'id':id,'left':left_boundary,'right':right_boundary}

ret=[]
for k in siteContents.keys():
    v=siteContents[k]
    name=''
    for tfName in v['TFs']:
        name=name+tfName+'-'
    name=name[0:len(name)-1] #delete trailing -
    ret.append((k,name,v['left'],v['right']))
    assert (len(refFasta) >= left >= 0 and 0 <= right <= len(refFasta))

pickle.dump(ret, open('../convertToPickled/TFBS_pickled.txt', 'wb'))
print len(ret)
cnx.close()
