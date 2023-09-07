suppressPackageStartupMessages(library("optparse"))

#get options 
option_list <- list(
        make_option("--infile", action="store",default=NULL, help="The input file for analysis"),
        make_option("--outdir", action="store",default="./", help="The output dirctory,  [default %default]"),
        make_option("--out", action="store",default="Chao1_index", help="The output filename,  [default %default]")
)

opt<-parse_args(OptionParser(usage="%prog [options] --infile chao1.avg.txt ", option_list=option_list))

if(is.null(opt$infile)){
    cat ("Use  %prog -h for more help info\nThe author:Junru Chen, chenjr@geneworks.cn\n")
    quit("no")
}

#set options

outdir<-paste(opt$outdir,"/",sep="")
if(!file.exists(outdir)){   dir.create(outdir) }
setwd(outdir)

library(ggplot2)

pdf(file=paste(opt$out,".pdf",sep=""))
chao1 <- read.table(opt$infile,header=T,sep="\t")
colnames(chao1) <- c("Sequences_Per_Sample","Sample","chao1")
p <- ggplot(chao1,aes(x= Sequences_Per_Sample,y=chao1,color=Sample))
p + geom_line(aes(color=factor(Sample)))+ labs(x='Sequences Per Sample',y=opt$out)+theme(legend.background=element_blank(),panel.background = element_rect(fill="white",colour="black",linetype=1,size=1))
dev.off()

#png(file=paste(opt$out,".png",sep=""),type="cairo")
#p + geom_line()+labs(x='Sequences Per Sample',y=opt$out)+theme(legend.background=element_blank(),panel.background = element_rect(fill="white",colour="black",linetype=1,size=1))
#dev.off()

outpdf <- paste(opt$out,".pdf",sep="")
outpng <- paste(opt$out,".png",sep="")
system(paste("cd",outdir,";convert -density 300",outpdf, outpng))
