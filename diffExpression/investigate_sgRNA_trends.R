library('edgeR')
library('reshape')
library(ggplot2)

allData <- read.csv('selected.csv')
counts <- data.frame(allData$X1_aerobic_t30_1_control, allData$X1_aerobic_t60_1_control, allData$X1_aerobic_t180_1_control,
                     allData$X1_aerobic_t30_2_chlor, allData$X1_aerobic_t60_2_chlor, allData$X1_aerobic_t180_2_chlor, row.names=allData$seq)
freqs <- data.frame(apply((counts+1), 2, function(x){return(x/sum(x))}) )
ratios <- data.frame(t30=freqs[,4]/freqs[,1], t60=freqs[,5]/freqs[,2], t180=freqs[,6]/freqs[,3], row.names=row.names(counts))
counts$t30 = freqs[,4]/freqs[,1]
counts$t60 = freqs[,5]/freqs[,2]
counts$t180 = freqs[,6]/freqs[,3]

sums = apply(counts[,1:6], 1, sum)
hist(log2(sums))
sum(sums>500)
subset <- (counts$t30 < 2^-2) & sums>300
sum(subset)
pairs(log2(counts[subset,7:9]))
mean(counts[subset,1])
mean(counts[,1])
