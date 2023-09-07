library(pheatmap)
data <- read.table("factor.tax.corr.screen.txt",sep="\t",header=T,row.names=1)
p <- read.table("factor.tax.corr.p.mark.xls",sep="\t",header=T,row.names=1)
pheatmap(data, display_numbers=p,notecol="white",fontsize=5,cellwidth=10,cellheight=10, number_color = "white", fontsize_number = 6, cluster_rows =T,file="corr.pdf",color=colorRampPalette(c("green","black","red"))(100))
