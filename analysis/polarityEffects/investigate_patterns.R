df = read.csv('genes_vs_operon_pos.csv', row.names=1)
hist(df[,2])

subset = !is.na(df$pos3) & df$pos3 < -3
sum(subset)
t = df[subset,]
