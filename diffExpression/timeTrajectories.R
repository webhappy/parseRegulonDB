library('edgeR')
library(ggplot2)
library(plyr)
library(reshape2)
allData <- read.csv('All data - 32992 with LRs.csv')

# time trajectories ----
times <- cbind(allData$count_t0_3, allData$count_t4_3,allData$count_t8_3, allData$count_t13_3)
labels <- factor(c(1,2,3,4),levels=c(1,2,3,4),labels=c('t0','t4','t8','t13'))
y <- DGEList(counts=times,group=labels)
y <- calcNormFactors(y)
pseudoCounts <- equalizeLibSizes(y)
effSize <- pseudoCounts$common.lib.size
effFraction <- (pseudoCounts$pseudo.counts) / effSize
# preserved <- effFraction[effFraction[,1] > 1e-10 & effFraction[,3] > 1e-10 & effFraction[,2] > 1e-10 & effFraction[,4] > 1e-10,]
ratios <- (pseudoCounts$pseudo.counts[,2:4]+1) /(1+pseudoCounts$pseudo.counts[,1])
colnames(ratios) <- c("T4","T8",'T13')

#melt
ratiosLong <- melt(ratios)
ratiosLong$Var2 <- factor(ratiosLong$Var2,c("T4","T8",'T13'),ordered=TRUE)

# investigate only extreme trajectories ----
ggplot(data=ratiosLong[log2(ratios[,3]) < -3,],aes(x=Var2,y=log2(value)))+geom_violin()+xlab('Time point')
ggplot(data=ratiosLong[log2(ratios[,3]) > -3 & log2(ratios[,3]) < -2 ,],aes(x=Var2,y=log2(value)))+geom_violin()+xlab('Time point')


# plot time series, % change versus time ----
subset <- effFraction[log2(ratios[,3]) < -2,]
colnames(subset) <- c('T0','T4','T8','T13')
#Divide by final value
progress <- (subset[,1:4]-subset[,1])/(subset[,4]-subset[,1])
temp <- melt(progress)
temp$Var2 <- factor(temp$Var2,c("T0",'T4','T8','T13'),c("T=0",'T=4 hrs','T=8 hrs','T=13 hrs'),ordered=TRUE)
dat <- ddply(temp,"Var2",summarise,mean=mean(value),sd=sd(value))
ggplot(temp,aes(x=temp$Var2,y=temp$value))+ylab("Average progress to final fraction")+xlab("Time point")+scale_y_continuous(breaks=seq(0,1,.2))+stat_summary(fun.data='mean_sdl',mult=1)+geom_point(data=dat,aes(x=Var2,y=mean,size=c(3)))+theme(legend.position="none")
ggplot(dat,aes(x=Var2))+geom_point(aes(x=Var2,y=mean))+geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd))



ggplot(data=ratiosLong[log2(ratioT13) < -1 & log2(ratioT13) > -2 ,],aes(x=Var2,y=log2(value),group=Var1))+geom_line()+stat_summary(fun.y=mean, colour="red", geom="line", aes(group = 1))


par(mfrow=c(1,1))
hist(log2((ratioT8/ratioT13)[log2(ratioT13) < -3 ]),main="Counts at T8 vs counts at T13 for subset of sgRNAs")
hist(log2((ratioT4/ratioT13)[log2(ratioT13) < -3 ]),main="Counts at T4 vs counts at T13 for subset of sgRNAs")

#colnames(ratios) <- 1:3
adjCounts <- cbind(times[,1] * y$samples$norm.factors[1], times[,2] * y$samples$norm.factors[2], times[,3] * y$samples$norm.factors[3])



length(ratioT13)
sum(sub)

# copy-pasted helper ####
# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}