library(ggplot2)
data <- read.csv('offTargetResults.csv')
summary(data)

qplot(data$mean,color=data$hasOffTargets,geom='density',xlab='Mean LR of each sgRNA')+theme_bw()
qplot(data$mean,color=data$onlyPerfectOffTargets,geom='density')+theme_bw()

portion <- data[data$numOff>0 & data$mean < -2,]
lm(portion$mean ~ portion$numOff)
plot(portion$numOff,portion$mean,main='Does having more off-targets increase fitness defect?',sub='Limiting to 51 sgRNAs that hit more than once w/ fitness defect stronger than -2',ylab='Mean LR of each sgRNA',xlab='# of places sgRNA hits (can be non-exact)')
