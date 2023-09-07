suppressPackageStartupMessages(library("optparse"))

#get options 
option_list <- list(
        make_option("--indm", action="store",default=NULL, help="The input dm file for pcoa analysis"),
        make_option("--outdir", action="store",default="./", help="The output dirctory,  [default %default]"),
        make_option("--group", action="store",default=NULL,  help="The group info for drawing,format=samplename group")
)

opt<-parse_args(OptionParser(usage="%prog [options] --indm dm.txt --group group.list ", option_list=option_list))

if(is.null(opt$indm) || is.null(opt$group) ){
    cat ("Use  %prog -h for more help info\nThe author: Junru Chen, chenjr@geneworks.cn\n")
    quit("no")
}

#set options

outdir<-paste(opt$outdir,"/",sep="")
if(!file.exists(outdir)){   dir.create(outdir) }

#main scripts
    setwd(outdir)
	library(ggplot2)
    library(vegan)
	data = read.table(opt$indm)
    d=as.dist(data)
	groups = read.table(opt$group, head=F,colClasses=c("character","character"))
	colnames(groups) <- c("sample","group")
    length=length(unique(as.character(groups$group)))
    times1=length%/%8
    res1=length%%8
    times2=length%/%5
    res2=length%%5
    col1=rep(1:8,times1)
    col=c(col1,1:res1)
    pich1=rep(c(15:18,20,7:14,0:6),times2)
    pich=c(pich1,15:(15+res2))
    
# Dendorgram
    h = hclust(d, "average");
    hcd = as.dendrogram(h)
    labelColors = c("#CDB380", "#036564", "#EB6841", "#EDC951","red","green","blue","yellow","black","purple")
    clusMember  <- groups$group
    unique_g <- unique(clusMember)
    n <- 1
    for (i in unique_g) {
        pos <- grep(paste('^',i,'$',sep=""),clusMember)
        clusMember[pos] <- n
        n <- n+1
    }
    clusMember <- as.numeric(clusMember)
    attr(clusMember,"names") <- groups$sample
    colLab <- function(n) {
        if (is.leaf(n)) {
            a <- attributes(n)
            labCol <- labelColors[clusMember[which(names(clusMember) == a$label)]]
            attr(n, "nodePar") <- c(a$nodePar, lab.col = labCol)
        }
        n
    }
    clusDendro = dendrapply(hcd, colLab)
    pdf("Dendrogram.pdf",width=12)
    if(length(unique(groups$group)) == length(unique(groups$sample))){
        plot(hcd,sub="",xlab="",ylab="",horiz=T)
    }else{ plot(clusDendro, horiz=T)}
    dev.off()

# PCoA   		
	pca = cmdscale(d,k=3,eig=T) 
    PC1=pca$points[,1]
    PC2=pca$points[,2]
    write.csv(pca$points,file="PCoA.csv")

    ncol=ncol(groups)
    group1=c()
    group2=c()
    for(i in 1:length(groups$sample)){
        Order=grep(paste0('^',rownames(pca$points)[i],'$'),groups$sample,perl=T)
        group1[i]=groups$group[Order]
        if(ncol==3){
            group2[i]=groups$group2[Order]
        }
    }
    group1=factor(group1,levels=unique(group1))
    group2=factor(group2,levels=unique(group2))
    if(ncol==2){
        plotdata = data.frame(rownames(pca$points),PC1,PC2,group1)
        colnames(plotdata)=c("sample","PC1","PC2","group")
    }else if(ncol==3){
        plotdata = data.frame(rownames(pca$points),PC1,PC2,group1,group2)
        colnames(plotdata)=c("sample","PC1","PC2","group1","group2")
    }
	plotdata$sample = factor(plotdata$sample)
	plotdata$PC1=as.numeric(as.vector(plotdata$PC1))
	plotdata$PC2=as.numeric(as.vector(plotdata$PC2))
	pc1 =floor(pca$eig[1]/sum(pca$eig)*10000)/100
	pc2 = floor(pca$eig[2]/sum(pca$eig)*10000)/100
	pc3 = floor(pca$eig[3]/sum(pca$eig)*10000)/100
	pc1
	pc2
	pc3

	p2<-ggplot(plotdata, aes(PC1, PC2)) +
        geom_point(aes(colour=group,shape=group),size=6)+ 
        scale_shape_manual(values=pich)+
        scale_colour_manual(values=col)+
        geom_hline(yintercept =0,colour="grey80", linetype="dashed") + geom_vline(xintercept =0,colour="grey80", linetype="dashed") +
        labs(title="PCoA - PC1 vs PC2") + xlab(paste("PC1 ( ",pc1,"%"," )",sep="")) + ylab(paste("PC2 ( ",pc2,"%"," )",sep=""))+
        theme(text=element_text(family="Arial",size=18))+
       	theme(panel.background = element_rect(fill='white', colour='black'), panel.grid=element_blank())
	
	cairo_pdf("PCoA12.pdf",height=12,width=15)
	p2
    dev.off()
	

    p5<-ggplot(plotdata, aes(PC1, PC2)) +
        geom_text(aes(label=sample),size=5,family="Arial",hjust=0.5,vjust=-1)+ 
        geom_point(aes(colour=group,shape=group),size=6)+ 
        scale_shape_manual(values=pich)+
        scale_colour_manual(values=col)+
        labs(title="PCoA - PC1 vs PC2") + xlab(paste("PC1 ( ",pc1,"%"," )",sep="")) + ylab(paste("PC2 ( ",pc2,"%"," )",sep=""))+
        theme(text=element_text(family="Arial",size=18))+
        geom_hline(yintercept =0,colour="grey80", linetype="dashed") + geom_vline(xintercept =0,colour="grey80", linetype="dashed") +
       	theme(panel.background = element_rect(fill='white', colour='black'), panel.grid=element_blank())

        cairo_pdf("PCoA12-2.pdf",height=12,width=15)
    p5
    dev.off()
