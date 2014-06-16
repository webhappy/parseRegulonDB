library('edgeR')
library(ggplot2)
library(plyr)
library(reshape2)
allData <- read.csv('All data - 32992 with LRs.csv',quote="")
counts <- allData[,c(5,6,8,15:17)] 
counts <- allData[,c(44:45,41:42)] # aerobic data
group <- factor(c(1,1,2,2))
row.names(counts) <- allData$seq
y <- DGEList(counts=counts,group=group,genes=allData$seq)
y <- calcNormFactors(y)
y <- estimateCommonDisp(y)
y <- estimateTagwiseDisp(y)
et <- exactTest(y)
topTags(et)
head(counts)
fdr <- p.adjust(et$table[,3],method='fdr')
frame()
par(mfrow=c(2,1),oma=c(1,2,0,0)+.1,mar=c(4,0,1,1)+.5)
hist(fdr,50,xlab="FDR per sgRNA",main="")
hist(et$table[,3],50,xlab="P-val per sgRNA",main="")

par(mfrow=c(1,1))
sum(et$table[,3]<.05)
hist(et$table[et$table[,3]<.05,1],xlab="log fold-change for sgRNAs with P-val < .05",main="Log fold-change of 8,026 sgRNAs")
hist(et$table[fdr<.05,1],xlab="log fold-change for sgRNAs with FDR < .05")
hist(exp(fdr))

# Save data into MySQLDB ----
library(RMySQL)
con <- dbConnect(MySQL(),dbname="t")
finalTable <- data.frame(allData$RegulonDBID,allData$pos,as.character(allData$seq),et$table$PValue,et$table$logFC,apply(counts[,1:length(counts)],1,sum))
colnames(finalTable) <- c("RegulonDBID","pos","seq","pval",'logFC','countsSum')
#dbWriteTable(con,'anaerobic_PVALS',finalTable,overwrite=TRUE,row.names=FALSE)
dbWriteTable(con,'aerobic_PVALS',finalTable,overwrite=TRUE,row.names=FALSE)

#Look at sgRNAs with low counts and weak p-vals
counts$sum <- apply(counts,1,sum)
finalTable$seq <- as.character(finalTable$seq)
finalTable <- merge(finalTable,subset(counts,select='sum'),by.y='row.names',by.x='seq')
weakChanges <- merged[merged$pval > .05,]

#Look at big fold-changes only (find # sgRNAs with strong effect and low p-val)----
big <- finalTable[et$table$logFC < -2,]
hist(big$pval)
lowProbBigFC <- big[big$pval >= .05,]
rownames(lowProbBigFC) <- as.character(lowProbBigFC$seq)
t <- merge(lowProbBigFC,counts,by='row.names')

significantSgRNAs <- finalTable[finalTable$pval < .05,]
