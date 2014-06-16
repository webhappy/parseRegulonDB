library('edgeR')
library(ggplot2)
library(plyr)
library(reshape2)
data_0502 <- read.csv('20140502.csv',quote="")
rownames(data_0502) <- data_0502$seq
allData <- read.csv('All data - 32992 with LRs.csv',quote="")
rownames(allData) <- allData$seq

initCounts <- function (counts) {
  group <- factor(c(1,2))
  row.names(counts) <- data_0502$seq
  y <- DGEList(counts=counts,group=group,genes=data_0502$seq)
  y <- calcNormFactors(y)
  counts$before <-(counts[,1]+1) * y$samples[1,'norm.factors']
  counts$after <- (counts[,2]+1) * y$samples[2,'norm.factors']
  return(counts)
}

calcFC <- function (counts) {
  return(data.frame(counts$after / counts$before, row.names=row.names(counts)))
}

counts30 <- data.frame(data_0502$count_aerobic.t0.1_control, data_0502$count_aerobic.t30.1_control)
counts30 <- initCounts(counts30)
FC30 <- calcFC(counts30)

counts60 <- data.frame(data_0502$count_aerobic.t0.1_control, data_0502$count_aerobic.t60.1_control)
counts60 <- initCounts(counts60)
FC60 <- calcFC(counts60)

countsNor <- data.frame(data_0502$count_aerobic.t0.1_control, data_0502$count_aerobic.t180.6_Nor)
countsNor <- initCounts(countsNor)
FCNor <- calcFC(countsNor)

### prepare output file for web UI
countsNor <- data.frame(data_0502$count_aerobic.t0.1_control, data_0502$count_aerobic.t180.6_Nor)
y <- DGEList(counts=countsNor,group=factor(c(1,2)),genes=data_0502$seq)
y <- calcNormFactors(y)
et <- exactTest(y, dispersion=.1^2)
counts <- countsNor # alias counts to actual variable
message <- paste( counts[,1],' -> ',counts[,2],sep='')
finalTable <- data.frame(as.character(data_0502$seq),et$table$PValue,et$table$logFC,message)
colnames(finalTable) <- c('Seq','Pval','LogFC','message')
write.csv(finalTable,'0502_nor.csv',row.names=FALSE)

# investigate
counts <- merge(allData[,c('seq','count_t4_3')], data_0502[,c('seq','count_aerobic.t0.1_control','count_aerobic.t30.1_control',
                   'count_aerobic.t60.1_control', 'count_aerobic.t180.1_control')],by='seq')
counts <- data.frame(counts[,2:ncol(counts)],row.names=counts[,1])
temp <- DGEList(counts=counts,group=factor(c(1:5)),genes=data_0502$seq)
temp <- calcNormFactors(temp)
countsNormalized <- countsOnly
rownames(countsNormalized) <- data_0502$seq

### prep normalized counts
countsOnly <- data_0502[,grep('^count',colnames(data_0502))]
temp <- DGEList(counts=countsOnly, group=factor(1:7), genes=row.names(countsOnly))
temp <- calcNormFactors(temp)
for (i in 1:ncol(temp)) {  countsNormalized[,i] <- countsOnly[,i] * temp$samples$norm.factors[i] }
countsNormalized <- countsNormalized[,c(6,4,2,7,1,3,5)]
colnames(countsNormalized) <- c('t0_control','t30_control','t60_control','t180_control','t30_Rif','t180_Nor','t180_Gent')

countsAdj <- countsNormalized + 1
countsAdj <- log2(countsAdj)
