############################################################
#imports and parameters
############################################################


library(inline)
source("../_functions.R")
library(pheatmap)
#create a text file containing the gene symbols of DDR genes
#taken from: https://www.mdanderson.org/documents/Labs/Wood-Laboratory/human-dna-repair-genes.html
ddrgenes=as.character(read.table("signatures/DDR_genes_mdanderson.txt",header=T,sep="\t")[,1])
data=read.table("RNASeq_U2OS_data.xls",header=TRUE)

if(!dir.exists("Extended_figure_11a")){
	dir.create("Extended_figure_11a")
}








DGE=data[,grepl("DGE_",colnames(data))|grepl("ID_",colnames(data))]
pos_NA=is.na(DGE$ID_SYMBOL)
DGE=DGE[which(!pos_NA),]
rownames(DGE)=DGE$ID_SYMBOL

#extract only SYMBOLS of DDR genes list.
#MISSING: SHLD3 MRE11A SLX1B TFIIH
pos_DDR=match(ddrgenes,rownames(DGE))
pos_DDR=pos_DDR[!is.na(pos_DDR)]
DGE_DDR=DGE[pos_DDR,]

#write table of DDR genes with complete information
write.table(DGE_DDR,file="Extended_figure_11a/RNA_seq_DDR.xls",col.names=T,row.names=F,sep="\t",quote=F)


DEG_siCDK12=DGE_DDR[,grepl("_siCdk12$",colnames(DGE_DDR))][,c(2,6)]
colnames(DEG_siCDK12)=c("log2FoldChange","padj")
DEG_siCDK12_OHT=DGE_DDR[,grepl("_siCdk12_OHT$",colnames(DGE_DDR))][,c(2,6)]
colnames(DEG_siCDK12_OHT)=c("log2FoldChange","padj")
DEG_OHT=DGE_DDR[,grepl("_OHT$",colnames(DGE_DDR))& !grepl("_siCdk12",colnames(DGE_DDR))][,c(2,6)]
colnames(DEG_OHT)=c("log2FoldChange","padj")
listdegs=list(DEG_siCDK12,DEG_OHT,DEG_siCDK12_OHT)
names(listdegs)=c("DEG_siCDK12","DEG_OHT","DEG_siCDK12_OHT")


#excluded also: MSH5 SLX1A DNTT beause 0 expression or NA log2FC/padj
RNA_result_obj=matrixFromRNAresults(resultsList=listdegs,padj_thresh=1,log2FCthresh=0,centers=4)
plotRNAresults(matrixResults=RNA_result_obj,quantile_saturation=0.98,display_significance="none",
				palette=c("darkgreen","grey80","darkred"),signif_thresh=1,breaks=10,fileName="Extended_figure_11a/Extended_figure_11a.pdf",rownames_thresh=50)








