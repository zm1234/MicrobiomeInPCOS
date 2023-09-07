suppressPackageStartupMessages(library("optparse"))

#get options 
option_list <- list(
        make_option("--indir", action="store",default=NULL, help="The dirctory of Relative files"),
        make_option("--outdir", action="store",default="./", help="The output dirctory,  [default %default]"),
        make_option("--prefix", action="store",default="tax_table",  help="The prefix of the input directory's files,  [default %default]"),
        make_option("--top", action="store", type="integer", default=20,  help="The head taxonomy to show in the picture of bar,  [default %default]"),
        make_option("--heatmap", action="store", type="integer", default=20,  help="The head taxonomy to show in the picture of heatmap,  [default %default]"),
        make_option("--group", action="store",default=NULL,  help="The group info for heatmap,format=samplename group"),
        make_option("--mark", action="store",default="y",  help="analysis for heatmap, y|n default=y")
)

opt<-parse_args(OptionParser(usage="%prog [options] --indir", option_list=option_list))
if(is.null(opt$indir)){
    cat ("Use  %prog -h for more help info\nThe author: Junru Chen, chenjr@geneworks.cn\n")
    quit("no")
}

#functions
top.bar <- function(infile, outdir, outfile){
    if(!file.exists(outdir)){dir.create(outdir)}
	setwd(outdir)
	data <- read.table(infile,header=T,row.names=1,sep="\t")
	len=max(nchar(names(data)))
	colname <- colnames(data)[-length(data)]
	rowname <- head(rownames(data),n=opt$top)
	if(opt$top < dim(data)[1]){
		data <- head(as.matrix(data[,-dim(data)[2]]),n=opt$top)
		colnames(data) <- colname
		rownames(data) <- rowname
		Others <- 1- colSums(data)
		data <- rbind(data,Others)
		rownames(data)[dim(data)[1]] <- "Other"
	}else{
		data <- as.matrix(data[,-dim(data)[2]])
		colnames(data) <- colname
		rownames(data) <- rowname
	}
	data<-t(data)
	pdf(file=outfile,family="Helvetica")
        if(len>9){
        par(mar=c(5, 4, 3, 6))
        }else{
        par(mar=c(3, 4, 3, 6))
        }

#	par(mar=c(3, 4, 3, 6))
	x <- barplot(as.matrix(t(data)),legend.text =colnames(data),args.legend = list(xjust=0,cex=0.6),xaxt="n",col = col,las=1,ylab="Relative Abundance",cex.axis=0.8)
	axis(at=x,labels=row.names(data),xpd = T,cex.axis=0.6,font=2,side=1,las=2,)
	dev.off()
	#outfile <- gsub('.pdf$','.png',outfile,perl=T)
	#png(file=outfile,type="cairo",family="Helvetica")
	#par(mar=c(3, 4, 3, 6))
	#x <- barplot(as.matrix(t(data)),legend.text =colnames(data),args.legend = list(xjust=0,cex=0.6),xaxt="n",col = col,las=1,ylab="Relative Abundance",cex.axis=0.8)
	#axis(at=x,labels=row.names(data),xpd = T,cex.axis=0.6,font=2,side=1,las=2,)
	#dev.off()
	outpng <- gsub('.pdf$','.png',outfile,perl=T)
	system(paste("cd",outdir,";convert -density 300",outfile, outpng))
	outfile <- gsub('.pdf$','.xls',outfile,perl=T)
	write.table(t(data),file=outfile,sep="\t",row.names=T)
}

top.heatmap <- function(infile, outdir, outfile){
	if(!file.exists(outdir)){dir.create(outdir)}
	setwd(outdir)
	library(pheatmap)
    	data <-read.table(infile,sep="\t",header=T,row.names=1)
	if(opt$heatmap < dim(data)[1]){
    		data <- head(data[,-dim(data)[2]],n=opt$heatmap)
	}else{
		data <- data[,-dim(data)[2]]
		data <- data[-dim(data)[1],]
	}
	pheatmap(data,scale="row",filename=outfile) 
	
	#outfile <- gsub('.pdf$','.png',outfile,perl=T)
	#png(file=outfile,type="cairo")
	#pheatmap(data,scale="row",filename=NA) 
	#dev.off()
	outpng <- gsub('.pdf$','.png',outfile,perl=T)
	system(paste("cd",outdir,";convert -density 300",outfile, outpng))

	if (file.exists("Rplots.pdf")){file.remove("Rplots.pdf") }
	outfile <- gsub('.pdf$','.xls',outfile,perl=T)
	write.table(data,file=outfile,sep="\t",row.names=T)
}

top.heatmap.group <- function(infile, outdir, outfile,groupfile){
	if(!file.exists(outdir)){dir.create(outdir)}
	setwd(outdir)
	library(pheatmap)
    	data <-read.table(infile,sep="\t",header=T,row.names=1)
	if(opt$heatmap < dim(data)[1]){
    		data <- head(data[,-dim(data)[2]],n=opt$heatmap)
	}else{
		data <- data[,-dim(data)[2]]
		data <- data[-dim(data)[1],]
	}
    	group<-read.table(groupfile,sep="\t",header=F)
	colnames(group) <- c("sample","group")
	annotation_col = data.frame(Group=factor(group$group))
	rownames(annotation_col) = group$sample
	pheatmap(data,scale="row",annotation_col = annotation_col,filename=outfile,cluster_cols=F) 
	
	#outfile <- gsub('.pdf$','.png',outfile,perl=T)
	#png(file=outfile,type="cairo")
	#pheatmap(data,scale="row",annotation_col = annotation_col,filename=NA) 
	#dev.off()
	outpng <- gsub('.pdf$','.png',outfile,perl=T)
	system(paste("cd",outdir,";convert -density 300",outfile, outpng))

	if (file.exists("Rplots.pdf")){file.remove("Rplots.pdf") }
	outfile <- gsub('.pdf$','.xls',outfile,perl=T)
	write.table(data,file=outfile,sep="\t",row.names=T)
}

#set options
indir <- opt$indir
top <- opt$top
heatmap <- opt$heatmap
prefix <- opt$prefix
bar_cols_num <- opt$top + 1
outdir<-paste(opt$outdir,"/",sep="")
if(!file.exists(outdir)){   dir.create(outdir) }
bar_outdir <- paste(outdir,"top","/",sep="")
if(!file.exists(bar_outdir)){   dir.create(bar_outdir) }
heatmap_outdir <- paste(outdir,"heatmap","/",sep="")
if(!file.exists(heatmap_outdir)){   dir.create(heatmap_outdir) }

#set cols
col <- rainbow(bar_cols_num)
col[1] <- "#91b1d7"
col[2] <- "#ff9970"
col[5] <- 'lightgreen'
col[6] <- '#89cff0'
col[8] <- '#f4c2c2'
col[9] <- '#ffe135'
col[11] <- '#fe6f5e'
col[15] <- '#007aa5'
col[19] <- '#9f8170'
col[bar_cols_num] <- "grey90"

#bar picture for ranks
ranks <- c("phylum","class","order","family","genus")
ranks_brief <- c("p","c","o","f","g")

for (i in 1:5) {
	rank <- ranks[i]
	rank_brief <- ranks_brief[i]
	outdir_g <- paste(bar_outdir,rank,sep="")
	infile_g <- paste(indir,"/",prefix,".",rank_brief,".txt",sep="")
	outfile_g<- paste(bar_outdir,rank,"/",rank_brief,opt$top,".relative.dis.pdf",sep="")
	#bar plot
	top.bar(infile=infile_g, outdir=outdir_g, outfile=outfile_g)
	#heatmap plot
	outdir_heatmap_g <- paste(heatmap_outdir,rank,sep="")
	outfile_heatmap_g<- paste(heatmap_outdir,rank,"/",rank_brief,opt$heatmap,".relative.heatmap.pdf",sep="")
	if(is.null(opt$group)){
		if(opt$mark == "y"){
			top.heatmap(infile=infile_g, outdir=outdir_heatmap_g, outfile=outfile_heatmap_g)	
		}
	}else{
		if(opt$mark == "y"){
			top.heatmap.group(infile=infile_g, outdir=outdir_heatmap_g, outfile=outfile_heatmap_g,group=opt$group)
		}
	}
}




