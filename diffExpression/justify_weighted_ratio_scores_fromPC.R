library(ggplot2)
library(reshape2)

allData <- read.csv('selected.csv')
counts <- data.frame(allData$X1_aerobic_t30_1_control, allData$X1_aerobic_t60_1_control, allData$X1_aerobic_t180_1_control,
                      allData$X1_aerobic_t30_2_chlor, allData$X1_aerobic_t60_2_chlor, allData$X1_aerobic_t180_2_chlor, row.names=allData$seq)
counts <- data.frame(allData$X1_aerobic_t30_1_control, allData$X1_aerobic_t60_1_control, allData$X1_aerobic_t180_1_control,
                     allData$X1_aerobic_t30_6_nor, allData$X1_aerobic_t60_6_nor, allData$X1_aerobic_t180_6_No, allData$X2_aerobic_t180_6_nor, row.names=allData$seq)
freqs <- data.frame(apply((counts+1), 2, function(x){return(x/sum(x))}) )
#ratios <- data.frame(t30=freqs[,4]/freqs[,1], t60=freqs[,5]/freqs[,2], t180=freqs[,6]/freqs[,3], row.names=row.names(counts))
counts$t30 = freqs[,4]/freqs[,1]
counts$t60 = freqs[,5]/freqs[,2]
counts$t180 = freqs[,6]/freqs[,3]
counts$t180_2 = freqs[,7]/freqs[,3]

##
countsChlor = counts

pairs(log2(counts[8:11]))
cor(log2(counts[8:11]))

t180_too_low <- ( log2(ratios$t60) - log2(ratios$t180) ) > 2
subset <- counts[t180_too_low,]
sumsOfCounts_subset <- apply(subset[,1:6],1,sum)
sumsofCounts <- apply(counts[,1:6],1,sum)
boxplot(data.frame(sumsofCounts,sumsOfCounts_subset))
# in this subset with high t60, how often is t30 close to t180
hist(log2(subset$t30/subset$t180))

countsLong <- melt(log2(counts[8:11]))
qplot(value, data=countsLong,color=variable, geom='density')
qplot(x=variable,y=value, data=countsLong, geom='boxplot')
