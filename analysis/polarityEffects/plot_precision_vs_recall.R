#data <- read.csv('performance_replicates2_3_aerobic.csv')
#data <- read.csv('replicates4and5_aerobic.csv')
data <- read.csv('essentialsPerformance.csv')
#data <- data[100:nrow(data),]
shading = colorRampPalette(c("#3300ff",'black'))(len)
plot(data$Specificity, data$Recall, xlab='Specificity=Found non-essentials/Known non-essential', ylab='Recall=Found essential/Known essential',col=shading)
plot(data$Specificity, data$Precision, xlab='Specificity=Found non-essentials/Known non-essential', ylab='Precision=Found essential/Called essential',col=shading)

library(ggplot2)
ggplot(data,aes(x=Recall,y=Specificity))+xlab("Recall = TP / # known essentials")+ylab("Specificty = True negative rate")+geom_point()+geom_text(data=data[seq(1,nrow(data),20),],aes(x=Recall,y=Specificity,label=Cutoff),hjust=0,vjust=1)
ggplot(data,aes(x=Recall,y=Precision))+xlab("Recall = TP / # known essentials")+ylab("Precision = TP / # positive")+geom_point()+geom_text(data=data[seq(1,nrow(data),20),],aes(x=Recall,y=Precision,label=Cutoff),hjust=0,vjust=1)
