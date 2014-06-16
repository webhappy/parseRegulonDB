allData <- read.csv('selected.csv')
library('edgeR')
library(ggplot2)

counts <- allData[,-1]
rownames(allData) <- allData$seq

temp <- DGEList(counts=counts, group=factor(1:ncol(counts)), genes=allData[,1])
temp <- calcNormFactors(temp)
countsNormalized <- counts
for (i in 1:ncol(temp)) {  countsNormalized[,i] <- counts[,i] * temp$samples$norm.factors[i] }

### investigate the Norflax data
plot(countsNormalized$X1_aerobic_t180_6_Nor, countsNormalized$X2_aerobic_t180_6_nor,xlab='Rep 1, t180_Norflox',ylab='Rep 2, t180_Norflox')
extremePt <- countsNormalized$X1_aerobic_t180_6_Nor > 40000
plot(countsNormalized$X1_aerobic_t180_6_Nor[!extremePt], countsNormalized$X2_aerobic_t180_6_nor[!extremePt],xlab='Rep 1, t180_Norflox',ylab='Rep 2, t180_Norflox')

### determine if Nor is changing

#countsNormalized <- countsNormalized + 1
subset <- countsNormalized[ , c('X1_aerobic_t30_6_nor','X1_aerobic_t60_6_nor','X1_aerobic_t180_6_Nor') ]
subset <- apply(subset, 2, function(x) x/countsNormalized$X1_aerobic_t0_1_control)
subset[countsNormalized$X1_aerobic_t0_1_control < .001, ] <- NA
summary(subset)
#subset <- subset[subset[,3]<500,]
boxplot(log(subset[,3]))

enriched <- allData[subset[,3]>10 & !is.na(subset[,3]), c('seq', 'X1_aerobic_t0_1_control','X1_aerobic_t30_6_nor','X1_aerobic_t60_6_nor','X1_aerobic_t180_6_Nor')]
depleted <-allData[subset[,3]<.1 & !is.na(subset[,3]) & allData[,'X1_aerobic_t0_1_control']>20, c('seq', 'X1_aerobic_t0_1_control','X1_aerobic_t30_6_nor','X1_aerobic_t60_6_nor','X1_aerobic_t180_6_Nor')]

apply(counts,2,function(x){sum(x==0)})
colSums(counts)

### convert to frequencies ####
par(mar=c(7,5,1,1))
boxplot(counts,las=2)
freqs <- apply(counts, 2, function(x){return(x/sum(x))})
freqs <- data.frame(freqs)

### convert to adjustedCounts ####
countsNormalized <- data.frame(apply(counts, 2, function(x){return(x*2000000/sum(x))}))

# adjusted pseudo counts ####
countsNormalized <- data.frame(apply(counts, 2, function(x){x=x+1;return(x*2000000/sum(x))}))
ratios <- data.frame(ratio30=countsNormalized$X1_aerobic_t30_2_chlor/countsNormalized$X1_aerobic_t30_1_control,
                     ratio60=countsNormalized$X1_aerobic_t60_2_chlor/countsNormalized$X1_aerobic_t60_1_control,
                     ratio180=countsNormalized$X1_aerobic_t180_2_chlor/countsNormalized$X1_aerobic_t180_1_control,
                     chlor30=countsNormalized$X1_aerobic_t30_2_chlor,chlor60=countsNormalized$X1_aerobic_t60_2_chlor,chlor180=countsNormalized$X1_aerobic_t180_2_chlor,
                     cont0=countsNormalized$X1_aerobic_t0_1_control,cont30=countsNormalized$X1_aerobic_t30_1_control,cont60=countsNormalized$X1_aerobic_t60_1_control,
                     cont180=countsNormalized$X1_aerobic_t180_1_control,row.names=allData$seq)
ratios <- data.frame(ratio30=countsNormalized$X1_aerobic_t30_2_chlor/countsNormalized$X1_aerobic_t30_1_control,
                     ratio60=countsNormalized$X1_aerobic_t60_2_chlor/countsNormalized$X1_aerobic_t60_1_control,
                     ratio180=countsNormalized$X1_aerobic_t180_2_chlor/countsNormalized$X1_aerobic_t180_1_control,
                     wt30=sqrt(countsNormalized$X1_aerobic_t30_2_chlor+countsNormalized$X1_aerobic_t30_1_control),
                     wt60=sqrt(countsNormalized$X1_aerobic_t60_2_chlor+countsNormalized$X1_aerobic_t60_1_control),
                     wt180=sqrt(countsNormalized$X1_aerobic_t180_2_chlor+countsNormalized$X1_aerobic_t180_1_control)                     )
scores <- apply(ratios, 1, function(x){(x[1]*x[4]+x[2]*x[5]+x[3]*x[6])/sum(x[4:6])})
ratios <- data.frame(ratio30=countsNormalized$X1_aerobic_t30_2_chlor/countsNormalized$X1_aerobic_t30_1_control,
                     ratio60=countsNormalized$X1_aerobic_t60_2_chlor/countsNormalized$X1_aerobic_t60_1_control,
                     ratio180=countsNormalized$X1_aerobic_t180_2_chlor/countsNormalized$X1_aerobic_t180_1_control,
                     chlor30=countsNormalized$X1_aerobic_t30_2_chlor,chlor60=countsNormalized$X1_aerobic_t60_2_chlor,chlor180=countsNormalized$X1_aerobic_t180_2_chlor,
                     cont0=countsNormalized$X1_aerobic_t0_1_control,cont30=countsNormalized$X1_aerobic_t30_1_control,cont60=countsNormalized$X1_aerobic_t60_1_control,
                     cont180=countsNormalized$X1_aerobic_t180_1_control,row.names=allData$seq)
ratios$scores <- scores
pairs(log2(ratios[,1:3]))


enriched <- freqs[,3]/freqs$X1_aerobic_t0_1_control

