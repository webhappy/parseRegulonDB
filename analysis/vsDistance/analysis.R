library(ggplot2)
data=read.csv('fitnessVsDist.csv',header=FALSE)
colnames(data) <- c('gene','position','effect')
plot(data$position,data$effect)

ggplot(data,aes(x=position,y=effect,color=gene))+geom_line()+xlim(-20,100)

plot(factor(table(data$gene)),xlab="# of sgRNAs in this gene",ylab="# of genes")

fits=read.csv('fittingPerGene.csv',header=FALSE)
colnames(fits) <- c('gene','intercept','slope')
hist(fits$slope,breaks=20,xlim=c(-.03,.03),main="Slopes of linear fit for each gene of effect vs distance from ATG")
