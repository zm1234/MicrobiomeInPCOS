suppressPackageStartupMessages(library("optparse"))

#get options 
option_list <- list(
        make_option("--input1", action="store",default=NULL, help="The input2"),
	make_option("--input2", action="store",default=NULL, help="The input2"),
        make_option("--prefix", action="store",default="./venn", help="The output file prefix,  [default %default]"),
	make_option("--title", action="store",default="F vs S", help="title")
)

opt<-parse_args(OptionParser(usage="%prog [options] --input tax.txt", option_list=option_list))
if(is.null(opt$input1) || is.null(opt$input2)){
    cat ("Use  %prog -h for more help info\n")
    quit("no")
}
outdir<-paste(dirname(opt$prefix),"/",sep="")
if(!file.exists(outdir)){   dir.create(outdir) }

prefix=opt$prefix

library(VennDiagram)
venn.plot <- function(dataList, groupList, outfile, title){
# Data for Draw
for(i in 1:length(groupList)){
        assign(paste("venn_",groupList[i],sep=""),dataList[[i]])
}
zuhe <- groupList
input <- list()
for(k in 1:length(zuhe)){
        if(zuhe[k] %in% groupList){
            input <- append(input,list(get(paste("venn_",zuhe[k],sep=""))))
        }
}
names(input) <- zuhe[1:length(input)]
cols <- rainbow(length(input))
# Draw plot
venn <- venn.diagram(input,fill=cols,cat.fontface='bold',cat.col=cols,margin=0.1,filename=NULL)
# venn <- venn.diagram(input,fill=cols,cat.fontface='bold',cat.col=cols,margin=0.1,filename=NULL, main=title)
# Output
pdf(file=outfile)
grid.draw(venn)
dev.off()
}

print(opt$input1)
print(opt$input2)
data1<-read.table(opt$input1,header=T,row.names=1,check.names=F,sep="\t")
data2<-read.table(opt$input2,header=T,row.names=1,check.names=F,sep="\t")

# Group
groupList=c("F", "S")

topn=100
if (dim(data1)[1]>=topn && dim(data2)[1]>=topn){
dataList <- list()
dataList[[1]] <- rownames(data1)[1:topn]
dataList[[2]] <- rownames(data2)[1:topn]
outfile=paste0(prefix, '.top100.pdf')
venn.plot(dataList, groupList, outfile, title)
}

dataList <- list()
dataList[[1]] <- rownames(data1)
dataList[[2]] <- rownames(data2)
outfile=paste0(prefix, '.pdf')
venn.plot(dataList, groupList, outfile, paste0(title, ' (top', topn, ')'))

