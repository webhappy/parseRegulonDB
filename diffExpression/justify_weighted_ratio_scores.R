library('edgeR')
library('reshape')
library(ggplot2)

allData <- read.csv('selected.csv')
counts <- data.frame(allData$X1_aerobic_t30_1_control, allData$X1_aerobic_t60_1_control, allData$X1_aerobic_t180_1_control,
                      allData$X1_aerobic_t30_2_chlor, allData$X1_aerobic_t60_2_chlor, allData$X1_aerobic_t180_2_chlor, row.names=allData$seq)
counts <- data.frame(allData$X1_aerobic_t30_1_control, allData$X1_aerobic_t60_1_control, allData$X1_aerobic_t180_1_control,
                    allData$X1_aerobic_t30_6_nor, allData$X1_aerobic_t60_6_nor, allData$X1_aerobic_t180_6_Nor, row.names=allData$seq)
freqs <- data.frame(apply((counts+1), 2, function(x){return(x/sum(x))}) )
ratios <- data.frame(t30=freqs[,4]/freqs[,1], t60=freqs[,5]/freqs[,2], t180=freqs[,6]/freqs[,3], row.names=row.names(counts))
counts$t30 = freqs[,4]/freqs[,1]
counts$t60 = freqs[,5]/freqs[,2]
counts$t180 = freqs[,6]/freqs[,3]

# Assign to chlor ####
chlorCounts <- counts
chlorRatios <- ratios

pairs(log2(ratios))
cor(log2(ratios))

t180_too_low <- ( log2(ratios$t60) - log2(ratios$t180) ) > 2
subset <- counts[t180_too_low,]
pairs(log2(subset[7:9]))
sumsOfCounts_subset <- apply(subset[,1:6],1,sum)
sumsofCounts <- apply(counts[,1:6],1,sum)
hist(log2(sumsofCounts))
summary(sumsofCounts)
boxplot(sumsofCounts, sumsOfCounts_subset,names=c('Sums of counts for all sgRNAs','Sums of counts for subset'))
t.test(sumsofCounts,sumsOfCounts_subset)
# in this subset with high t60, how often is t30 close to t180
hist( log2(subset$t30/subset$t180) )

sgRNAs_with_lots_of_counts <- (sumsofCounts > 500)
pairs(log2(ratios[sgRNAs_with_lots_of_counts,]))
cor(log2(ratios[sgRNAs_with_lots_of_counts,]))

sgRNAs_with_low_counts <- (sumsofCounts<50)
sum(sgRNAs_with_low_counts)
pairs(log2(counts[sgRNAs_with_low_counts,7:9]))
cor(log2(counts[sgRNAs_with_low_counts,7:9]))

sgRNAs_with_high_change_at_t30 <- abs(log2(ratios$t30)) > 2 & abs(log2(ratios$t180))<.5
sum(sgRNAs_with_high_change_at_t30)
pairs(log2(counts[sgRNAs_with_high_change_at_t30,7:9]))
temp <- ratios[sgRNAs_with_high_change_at_t30, ]
temp2 <- abs(log2(ratios$t30)) > 2
temp2 <- abs(log2(ratios$t30))
temp3 <- abs(log2(ratios$t180)) < .5
temp <- counts[sgRNAs_with_high_change_at_t30,]
temp$log2 <- log2(temp$t30)
write.csv(counts[sgRNAs_with_high_change_at_t30,], 'over_active_at_t30.csv')

# Reshape ratios so I can plot overlay with ggplot ####
ratiosLong <- melt(log2(ratios))
qplot(value, data=ratiosLong, geom='density', color=variable, xlab='Log2 ratio', main='Comparing distribution of sgRNA effect for different time points')
counts180 = (counts$allData.X1_aerobic_t180_1_control+counts$allData.X1_aerobic_t180_6_Nor)
plot(log2(counts$t180) ~ log2(counts180) )

temp <- hist(log2(counts$t180),200)
which(temp$counts > 2500)
temp$mids[49]
temp$xname
