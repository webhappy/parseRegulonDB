library('RMySQL')
library('reshape')
library('plyr')
con <- dbConnect(MySQL(), dbname='t')
genes <- dbReadTable(con,'GENE')
first=genes[genes$operon_pos==1,]
hist(genes$avg_lr,main="Finding cutoff to define essential genes",xlab="Average fitness effect of all sgRNA/replicates for each gene")
#hist(first)
essentialFirst<-first[first$avg_lr < -3 & !is.na(first$avg_lr),]

exptResults <- dbReadTable(con,'EXPT_RESULTS')
exptCast <- cast(exptResults,pos+strand ~ replicate)
exptCast$means <- apply(exptCast[, 3:5],1,mean)

# gene contains gene_posleft, gene_posright, and gene_strand
getDistanceFromStartOfGene <- function (pos, strand, gene) {
  if ( strand == 'forward' ) {
    # subtract pos from gene's left
    return(pos - gene$gene_posleft)
  } else{
    return(gene$gene_posright - pos)
  }
}

points=data.frame()
addPoints <- function (gene) {
  cur <- exptCast[(exptCast$pos < gene$gene_posright) 
                  & exptCast$pos > gene$gene_posleft & 
                    exptCast$strand==gene$gene_strand,]
  # contains pos, strand, and means
  apply(cur,2, function (x){print(x)})#cat(getDistanceFromStartOfGene(x$pos,x$strand,gene)))
}

s=addPoints(essentialFirst[2,])
gene <- essentialFirst[2,]
cur <- exptCast[(exptCast$pos < gene$gene_posright) 
                & exptCast$pos > gene$gene_posleft & 
                  exptCast$strand==gene$gene_strand,]
