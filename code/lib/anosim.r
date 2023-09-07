data <- read.table("04.otu_table_with_taxonomy.txt",sep="\t",header=T,row.names=1)
groups <- read.table("group.list",sep="\t",header=F)
colnames(groups) <- c("samples","group")

library(vegan)
ar = anosim(t(data),groups$group,permutations = 999, distance = "bray")
summary(ar)
dc = data.frame(dis=ar$dis.rank, class=ar$class.vec)
