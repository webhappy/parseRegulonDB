library(ggplot2)
data <- read.csv('sgRNAS_in_TFBS.csv')
summary(data)
hist(data$mean)
temp <- cbind(rep(TRUE,nrow(data)),data$inTFBS=='False')
qplot(data$mean,geom='density',labels=c('Not in TFBS','In TFBS'),color=data$inTFBS,xlab='Mean LR score over 3 replicates for each sgRNA')

ggplot(data)+theme(legend.position='left')+geom_density(aes(x=mean,color='blue'),data=data[data$inTFBS=='True',],color='blue',alpha=.2)+geom_density(aes(x=mean,color='green'),color='green',data=data,alpha=.2)+scale_color_manual(values=c('In TFBS','All sgRNAs'))
