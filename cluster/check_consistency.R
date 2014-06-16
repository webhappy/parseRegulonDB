genes <- read.csv('consistency_of_genes.csv')
row.names(genes) <- genes$Gene
genes <- data.frame(genes[,-1])
boxplot(genes[,seq(2,12,2)],las=2)
boxplot(genes[,seq(1,11,2)],las=2)
boxplot(genes[,seq(2,12,2)]/genes[,seq(1,11,2)], las=2)
seq(1,11,2)

pairs(genes[,c(1,3,5)])
pairs(genes[,c(7,9,11)])