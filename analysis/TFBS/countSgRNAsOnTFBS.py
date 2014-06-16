import pandas.io.sql as psql
import pandas.io.parsers
import mysql.connector

cnx = mysql.connector.connect(user='root', host='127.0.0.1', database='t')
sites=psql.read_frame('select SITE.site_id,site_sequence,site_posleft,site_posright,ri_function,CONFORMATION.transcription_factor_id, transcription_factor_name '
                      'from SITE left join `REGULATORY_INTERACTION` ON SITE.`site_id`=`REGULATORY_INTERACTION`.`site_id` left join CONFORMATION on `REGULATORY_INTERACTION`.`conformation_id`=`CONFORMATION`.`conformation_id` left join TRANSCRIPTION_FACTOR ON CONFORMATION.`transcription_factor_id`=`TRANSCRIPTION_FACTOR`.`transcription_factor_id` where site_posleft>0 ORDER BY site_id asc',cnx)
# anaerobic_d24_1_LR
expData=pandas.io.parsers.read_csv(open('../../anaerobic.csv','rb'))
expData['mean']=expData[[elem for elem in expData.columns if elem[-3:]=='_LR']].mean(axis=1)
expData['inTFBS']=[False for x in xrange(len(expData))]

def getSgRNAsWithinRange (sgRNAs,left, right):
    selector=( ((sgRNAs['sgRNA_pos']>=left)&(sgRNAs['sgRNA_pos']<=right))
                      | ((sgRNAs['sgRNA_pos']+20<=right) & (sgRNAs['sgRNA_pos']+20>=left)) )
    candidates=sgRNAs[selector]
    sgRNAs.inTFBS[selector]=True
    return candidates[['seq','mean']]

hashtable={}
countGood=0
CUTOFF=-3.0
for (site_id,site_seq,site_left,site_right,function,tf_id,tf_name) in sites.itertuples(index=False):
    sgRNAs=getSgRNAsWithinRange(expData,site_left,site_right)
    #print len(sgRNAs),'sgRNAs within',site_left,'and',site_right
    for (seq,mean) in sgRNAs.itertuples(index=False):
        if seq not in hashtable:
            hashtable[seq]=mean
            if abs(mean)>abs(CUTOFF):
                print tf_name,'at',site_left,'and',site_right,'has sgRNA with value',mean
                countGood+=1

print len(hashtable),'items in sgRNA hashtable'
print countGood,'are below the CUTOFF',CUTOFF
expData.to_csv('sgRNAS_in_TFBS.csv')