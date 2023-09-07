suppressPackageStartupMessages(library("optparse"))

#get options 
option_list <- list(
        make_option("--infile", action="store",default=NULL, help="The input file for analysis"),
        make_option("--out", action="store",default=NULL, help="The output file")
)

opt<-parse_args(OptionParser(usage="%prog [options] --infile otu_table.txt --out out.pdf ", option_list=option_list))
if(is.null(opt$infile)){
    cat ("Use  %prog -h for more help info\nThe author:Junru Chen, chenjr@geneworks.cn\n")
    quit("no")
}

library(RColorBrewer)
otu.table <- read.table(opt$infile,header=T,sep="\t")
name<-colnames(otu.table)[2:(length(colnames(otu.table))-1)]
colors<-palette(rainbow(length(name)))
otu.tab.len<-length(colnames(otu.table))
otu.matrix<-as.matrix(otu.table[,2:(otu.tab.len-1)])
colnames(otu.matrix)<-name
otu.mat.col.len<-length(colnames(otu.matrix))
sum.data<-apply(otu.matrix,2,sum)
otu<-matrix(nrow=length(otu.matrix[,1]),ncol=otu.mat.col.len)
times=otu.mat.col.len%/%35
otu.order<-matrix(nrow=length(otu.matrix[,1]),ncol=otu.mat.col.len)
otu[,1:otu.mat.col.len]<-otu.matrix[,1:otu.mat.col.len]/sum.data[1:otu.mat.col.len]
for(i in 1:otu.mat.col.len){
  otu.order[,i]<- rev(sort(otu[,i]))
}
x.ranges<-vector(length=length(otu[1,]))
otu=otu.order
for(i in 1:otu.mat.col.len){
  x.ranges[i]<-length(otu[which(otu[,i]>0)])
}
otu.order<-matrix(nrow=length(otu.matrix[,1]),ncol=otu.mat.col.len)
x.range<-max(x.ranges)*1.1

cairo_pdf(file=opt$out)
x<-seq(1:length(otu[which(otu[,1]>0)]))
par(mai=c(1.2,1.2,0.8,0.8))
plot(x,log(otu[which(otu[,1]>0),1],base=10),axes=F,ylab="Relative Abundance",xlab="Species Rank",type='l',lwd=2,xlim=c(0,x.range),cex.lab=1,cex.axis=1.5)
axis(1,)
axis(2,at=0:-5,labels=10^(0:-5))
for(i in 1:otu.mat.col.len){
  x<-seq(1:length(otu[which(otu[,i]>0)]))
  lines(x,log(otu[which(otu[,i]>0),i],base=10),type='l',lty=1,lwd=2,col=colors[i])
}
legend("topright",lty=1,lwd=1.5,legend=name,cex=0.8,ncol=times+1,col=colors)
box()
dev.off()
