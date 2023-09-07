suppressPackageStartupMessages(library("optparse"))

#get options 
option_list <- list(
        make_option("--input", action="store",default=NULL, help="The input otu table for analysis"),
        make_option("--outdir", action="store",default="./", help="The output dirctory,  [default %default]"),
        make_option("--venn", action="store",default=NULL,  help="The venn list for analysis")
)

opt<-parse_args(OptionParser(usage="%prog [options] --input otu_table.txt --venn venn.list", option_list=option_list))
if(is.null(opt$input)){
    cat ("Use  %prog -h for more help info\nThe author: Junru Chen, chenjr@geneworks.cn\n")
    quit("no")
}
outdir<-paste(opt$outdir,"/",sep="")
if(!file.exists(outdir)){   dir.create(outdir) }

library(VennDiagram)
data<-read.table(opt$input,header=T,row.names=1,check.names=F,sep="\t")
data = as.matrix(data)
sample_list=colnames(data)
venn_list <- read.table(opt$venn,header=F,sep="\t",fill=T)
venn_list = as.matrix(venn_list)

for(i in 1:length(sample_list)){
	assign(paste("venn_",sample_list[i],sep=""),attr(which(data[,i]!=0),"names"))
}

for(i in 1:dim(venn_list)[1]){
	zuhe <- venn_list[i,]
	input <- list()
	for(k in 1:length(zuhe)){
		if(zuhe[k] %in% sample_list){
			input <- append(input,list(get(paste("venn_",zuhe[k],sep=""))))
		}
	}
	names(input) <- zuhe[1:length(input)]
	cols <- rainbow(length(input))
	venn <- venn.diagram(input,fill=cols,cat.fontface='bold',cat.col=cols,margin=0.1,filename=NULL)

	#for pdf
	name <- gsub('(-{2,5}|-$)','',paste(zuhe,collapse="-"),perl=T)
	out <- paste(i,"_",name,".pdf",sep="")
	pdf(file=out)
	grid.draw(venn)
	dev.off()

	#for png
	#out <- paste(i,"_",name,".png",sep="")
	#png(file=out,type="cairo")
	#grid.draw(venn)
	#dev.off()

	outpng <- paste(i,"_",name,".png",sep="")
	system(paste("cd",outdir,";convert -density 300",out, outpng))
}
