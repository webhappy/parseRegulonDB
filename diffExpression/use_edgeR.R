library('edgeR')

allData <- read.csv('../All data - 32992 with LRs.csv',quote="")
counts <- allData[,c('count_t0_2','count_t8_2')]
y <- DGEList(counts=counts,group=c(1,2),genes=allData$seq)
