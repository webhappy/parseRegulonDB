library(RMySQL)
con <- dbConnect(MySQL(),user="root",dbname="t")
genes <- dbReadTable(con, 'GENE')
# to hack in logFC
genes$logFC <- genes$median
hist(genes[genes$operon_pos==1,]$logFC)

genes <- genes[!is.na(genes$logFC) & genes$operon_pos>0,]
genes$operon_pos <- as.factor(genes$operon_pos)
table(genes$operon_pos)
hist(genes$avg_lr,xlab="Average fitness effect at each gene",main="Averaging over all replicates and sgRNA's for each gene")
summary(genes)
library(ggplot2)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# To use for fills, add
scale_fill_manual(values=cbPalette)
# To use for line and point colors, add
scale_color_manual(values=cbPalette)
data=genes[as.numeric(genes$operon_pos) < 7,]
qplot(data=data, avg_lr,geom='density',color=operon_pos)

sgRNAs <- read.csv('sgRNAs_in_essential_genes.csv')
sgRNAs$operon_pos <- as.factor(sgRNAs$operon_pos)
sgRNAs$geneName <- as.factor(sgRNAs$geneName)
qplot(sgRNAs$avg_lr,xlab="Average effect of each sgRNA")
qplot(sgRNAs$normalized_lr,xlab="Average effect of each sgRNA normalized to strongest sgRNA for each gene")+ggtitle(expression(atop("Normalized distribution of sgRNAs",atop("Each gene has one sgRNA with effect size=1"))))
summary(sgRNAs)
nlevels(sgRNAs$geneName)

badSgRNAs_normalized <- sgRNAs[sgRNAs$normalized_lr < .3, ]
badSgRNAs <- sgRNAs[ sgRNAs$avg_lr > -5,]
-16*.3
sgRNAs$isBad <- sgRNAs$normalized_lr < .3
sgRNAs$isBad <- factor(sgRNAs$isBad, levels=c(FALSE,TRUE),labels=c("Good sgRNA","Bad sgRNA"))

qplot(sgRNAs$normalized_lr,xlab="Average effect of each sgRNA normalized to strongest sgRNA for each gene",geom='density',color=sgRNAs$isBad)+scale_colour_manual(values=cbPalette)
qplot(sgRNAs$distance_ATG,xlab="Distance for each sgRNA from ATG",geom='density',color=sgRNAs$isBad)+scale_colour_manual(values=cbPalette)
boxplot(sgRNAs$distance_ATG ~ sgRNAs$isBad,ylab="Distance from ATG",xlab="sgRNA has normalized effect > .3")
plot(table(sgRNAs$operon_pos) ~ sgRNAs$isBad,ylab="Operon position",xlab="sgRNA has normalized effect > .3")
qplot(operon_pos,data=sgRNAs,facets=isBad ~ .,geom="histogram",xlab="Operon position for gene",main="Bad status is not a function of operon position")+scale_x_discrete()
