# Make resulting file
library('edgeR')
library(RMySQL)

con <- dbConnect(MySQL(),dbname="t")
sgRNAs <- dbReadTable(con,'sgRNA')

#Read all data----
allData <- read.csv('All data - 32992 with LRs.csv',quote="")

exportExptResults <- function (colNamesBefore,colNamesAfter,repNames,fileName) {
  # Generates a file that is ready to be used in the web-app for browsing in the Javascript interface
  # Notice that variable allData is used within here!
  # I assume all new experiments will be recorded as columns within the All data spreadsheet
  # This version generates a message assuming that Before and After are paired.
  #
  # Args:
  #  colNamesBefore: List with Names of columns to use in allData that contain BEFORE states
  #  colNamesAfter: List with Names of columns to use in allData that contains AFTER states
  #  repNames: List with String descriptions, only used to generate the message that is shown on mouse-over of the sgRNA
  #  fileName: String containing output file
  counts <- allData[,c(colNamesBefore,colNamesAfter)]
  group <- factor(c(rep(1,length(colNamesBefore)),rep(2,length(colNamesAfter))))
  row.names(counts) <- allData$seq
  y <- DGEList(counts=counts,group=group,genes=allData$seq)
  y <- calcNormFactors(y)
  y <- estimateCommonDisp(y)
  y <- estimateTagwiseDisp(y)
  et <- exactTest(y)
  allData <- merge(allData,sgRNAs,by='seq')
  et$table$Seq=rownames(et$table)
  allData <- merge(allData,et$table,by.x='seq',by.y='Seq')
  for (i in 1:length(repNames)) {
    if (i>1){
      message <- paste(message,'<br>')
    } else{
      message <- rep('',length(et$table$logFC))
    }
    message <- paste(message, repNames[i],': ',allData[,colNamesBefore[i]],' -> ',allData[,colNamesAfter[i]],sep='')
  }
  finalTable <- data.frame(as.character(allData$seq),allData$pos.y,allData$strand.y,allData$PValue,allData$logFC,message)
  colnames(finalTable) <- c('Seq','Pos','Strand','Pval','LogFC','message')
  write.csv(finalTable,fileName,row.names=FALSE)
}

exportExptResults(c('count_O2_t0_1','count_O2_t0_2','count_O2_t0_3'),c('count_O2_d24_1','count_O2_d24_2','count_O2_d24_3'),c('Replicate 1','Replicate 2','Replicate 3'),'anaerobic_0314.csv')
exportExptResults(c('count_t0_2','count_t0_3'),c('count_t8_2','count_t8_3'),c('Replicate 2','Replicate 3'),'aerobic_2and3.csv')
exportExptResults(c('count_t0_4','count_t0_5'),c('count_t8_4','count_t8_5'),c('Replicate 4','Replicate 5'),'aerobic_4and5.csv')
