library(ggplot2)

gross <- read.csv('gross.csv')
grossSubset <- data.frame(gross$ECKID,gross$Gene,gross$GLUCOSE...UNSPECIFIED,gross$GLYCEROL...UNSPECIFIED)
keio <- read.csv('keio.csv')
keio$LB_22hr[keio$LB_22hr=='N.A.']=NA
keio$LB_22hr <- as.numeric(as.character(keio$LB_22hr))
keio$MOPS_24hr[keio$MOPS_24hr=='N.A.']=NA
keio$MOPS_24hr <- as.numeric(as.character(keio$MOPS_24hr))
keioSubset <- data.frame(keio$eckID,keio$LB_22hr,keio$MOPS_24hr)

palsson <- read.csv('palsson_microarray_data_w_info_080826.csv')
palsson$X <- as.character(palsson$X)
palsson$X <- as.character(sapply(strsplit(palsson$X,split="_"),'[[',1))
palssonAerobicMean <- apply(cbind(palsson$ec_aer_wild_O_a_SciA.CEL,palsson$ec_aer_wild_O_b_SciA.CEL,palsson$ec_aer_wild_O_c_SciA.CEL),1,mean)
palssonAnaerobicMean <- apply(cbind(palsson$ec_aer_wild_nO_a_SciA.CEL,palsson$ec_aer_wild_nO_b_SciA.CEL,palsson$ec_aer_wild_nO_c_SciA.CEL,palsson$ec_aer_wild_nO_d_SciA.CEL),1,mean)
palssonSubset <- data.frame(tolower(palsson$X),palssonAerobicMean,palssonAnaerobicMean)

merged <- merge(grossSubset,keioSubset,by.x='gross.ECKID',by.y='keio.eckID',all=TRUE)
merged$gross.Gene <- tolower(merged$gross.Gene)
merged <- merged[!is.na(merged$gross.Gene),]
merged <- merge(merged, palssonSubset,by.x='gross.Gene',by.y='tolower.palsson.X.',all.x=TRUE)

#plots ----
hist(keio$LB_22hr)
hist(keio$MOPS_24hr)
hist(gross$GLUCOSE...UNSPECIFIED)
plot(merged$gross.GLUCOSE...UNSPECIFIED,merged$keio.LB_22hr)

library(RMySQL)
#pull in genes data from mySQL
con <- dbConnect(MySQL(),dbname='t')
genes <- dbReadTable(con,'GENE')
PEC <- dbReadTable(con,'PEC')
PEC$gene_name <- tolower(PEC$gene_name)
genes$gene_name <- tolower(genes$gene_name)
genes <- genes[!is.na(genes$gene_name),]
mergedWithGenes <- merge(merged,genes,by.x='gross.Gene',by.y='gene_name',all.y=T)
mergedWithGenes$isEssential <- mergedWithGenes$gross.Gene %in% PEC$gene_name
qplot(mergedWithGenes$gross.GLUCOSE...UNSPECIFIED,mergedWithGenes$keio.LB_22hr,color=mergedWithGenes$isEssential)
qplot(mergedWithGenes$keio.LB_22hr,geom='density',color=mergedWithGenes$isEssential)

#pull in list of cases we're investigating----
oddCases <- read.csv('oddCases.csv')
oddCases$Gene <- tolower(as.character(oddCases$Gene))
oddCases$Seq <- as.character(oddCases$Seq)
finalTable <- merge(oddCases,mergedWithGenes,by.x='Gene',by.y='gross.Gene',all.x=TRUE)
write.csv(finalTable,'oddCases_with_data.csv')

#look at the data----
strongChanges <- mergedWithGenes$mean < -1
sum(strongChanges)
strongChanges <- mergedWithGenes[strongChanges,]
plot(strongChanges$keio.LB_22hr,strongChanges$mean,xlab='Keio OD in LB @ 22 hrs',ylab='Mean of significant sgRNAs')
plot(strongChanges$gross.GLUCOSE...UNSPECIFIED,strongChanges$mean,xlab='Gross Glucose score',ylab='Mean of significant sgRNAs')
plot(strongChanges$palssonAnaerobicMean,strongChanges$mean,xlab='Palsson anaerobic mean microarray expression',ylab='Mean of significant sgRNAs')
