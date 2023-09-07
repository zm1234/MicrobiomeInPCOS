library(psych)

otu <- read.table("tax_table.txt",sep="\t",header=T,row.names=1)
otu <- otu[-dim(otu)[1],]
factor <- read.table("factor.input.txt",sep="\t",header=T,row.names=1)
occor = corr.test(t(otu),factor,alpha=.05,method="spearman",adjust="none")
#occor = corr.test(otu,use="pairwise",method="pearson",adjust="fdr",alpha=.05)
occor.r = occor$r # 取相关性矩阵R值
write.table(occor.r,file="factor.tax.corr.xls",sep="\t",col.names = NA, row.names = TRUE,quote = FALSE,fileEncoding = "UTF-8")
occor.p = occor$p # 取相关性矩阵p值
write.table(occor.p,file="factor.tax.corr.p.xls",sep="\t",col.names = NA, row.names = TRUE,quote = FALSE,fileEncoding = "UTF-8")
occor.r.ori <- occor.r
occor.p.ori <- occor.p
occor.r[occor.p>0.05] = 0
occor.r[is.na(occor.r)] = 0
write.table(occor.r,file="factor.tax.corr.screen.ori.txt",sep="\t",col.names = NA, row.names = TRUE,quote = FALSE,fileEncoding = "UTF-8")
mark <-which(as.numeric(rowSums(occor.r)) == 0)
if(length(mark) > 0 ){occor.r <- occor.r[-mark,]} # 产生无缺失值的行
mark <-which(as.numeric(colSums(occor.r)) == 0)
if(length(mark) > 0 ){occor.r <- occor.r[,-mark]} # 产生无缺失值的列
occor.r <- occor.r.ori[row.names(occor.r),colnames(occor.r)]
occor.p <- occor.p.ori[row.names(occor.r),colnames(occor.r)]
occor.p[is.na(occor.p)] = 1
write.table(occor.r,file="factor.tax.corr.screen.txt",sep="\t",col.names = NA, row.names = TRUE,quote = FALSE,fileEncoding = "UTF-8")
write.table(occor.p,file="factor.tax.corr.p.screen.txt",sep="\t",col.names = NA, row.names = TRUE,quote = FALSE,fileEncoding = "UTF-8")
