#import data

library(TxDb.Hsapiens.UCSC.hg19.knownGene)

data=read.table("RNASeq_U2OS_data.xls",header=TRUE)

#create a text file containing the gene symbols of DDR genes
#taken from: https://www.mdanderson.org/documents/Labs/Wood-Laboratory/human-dna-repair-genes.html
signature_DDR_genes_manderson="signatures/DDR_genes_mdanderson.txt"

if(!dir.exists("Extended_figure_11c")){
	dir.create("Extended_figure_11c")
}









#scatter correlation between gene length and log2FC upon siCDK12 (hihlight UP and DOWN genes 
#with different colors)
gene_list=genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
#match genes in RNA-Seq and gene_list
pos=match(gene_list$gene_id,data$ID_ENTREZ)
elementMetadata(gene_list)=data[pos,]
gene_list$gene_length=width(gene_list)
#filter only significant and only genes longer than 10000 bp
gene_list_signif=gene_list[!is.na(gene_list$DGE_padj_siCdk12)& gene_list$DGE_padj_siCdk12<0.01]
gene_list_signif=gene_list_signif[gene_list_signif$gene_length>10000]

#only DDR genes. 
ddrgenes=as.character(read.table(signature_DDR_genes_manderson,header=T,sep="\t")[,1])
pos_DDR=match(ddrgenes,gene_list$ID_SYMBOL)
pos_DDR=pos_DDR[!is.na(pos_DDR)]
DGE_DDR=gene_list[pos_DDR]

write.table(as.data.frame(DGE_DDR),file="Extended_figure_11c/DDR_genes_complete.xls",
				sep="\t",quote=F,row.names=F,col.names=T)
pos_signif_UP= (!is.na(DGE_DDR$DGE_padj_siCdk12)& DGE_DDR$DGE_padj_siCdk12<0.01)&DGE_DDR$DGE_log2FoldChange_siCdk12>0.5
pos_signif_DOWN= (!is.na(DGE_DDR$DGE_padj_siCdk12)& DGE_DDR$DGE_padj_siCdk12<0.01)&DGE_DDR$DGE_log2FoldChange_siCdk12<(-0.5)
pdf("Extended_figure_11c/Extended_figure_11c_left.pdf")
plot(log10(DGE_DDR$gene_length),DGE_DDR$DGE_log2FoldChange_siCdk12,
		col=ifelse(pos_signif_DOWN,"blue",
					ifelse(pos_signif_UP,"red","black"))
		,pch=19,ylab="RNA-Seq log2FC siCDK12",xlab="log10 gene length (bp)",main=paste0("DDR genes: ",length(DGE_DDR)))
legend("topright",legend=c("DEG up 0.5","DEG down -0.5"),pch=19,col=c("red","blue"))
legend("bottomright",legend=c(
								paste("Pearson cor:",round(cor(log10(DGE_DDR$gene_length),DGE_DDR$DGE_log2FoldChange_siCdk12,use="complete.obs"),2)),
								paste("Spearman cor:",round(cor(log10(DGE_DDR$gene_length),DGE_DDR$DGE_log2FoldChange_siCdk12,use="complete.obs",method="spearman"),2))),
					lty=1,col="black")
dev.off()




DGE_DDR_UP=DGE_DDR[pos_signif_UP,]
DGE_DDR_DOWN=DGE_DDR[pos_signif_DOWN,]

#function for the dotplot (stripchart)
dotplot<-function(Objectlist,labs,widthlines=0.2,center="median",...){
	stripchart(Objectlist,vertical=TRUE,method="jitter",xaxt="n",...)
	axis(1,at=1:length(Objectlist),labels=labs)
	for(i in 1:length(Objectlist)){
		segments(x0=i-widthlines,x1=i+widthlines,y0=median(Objectlist[[i]]),y1=median(Objectlist[[i]]),lwd=2)
	}
}
pdf("Extended_figure_11c/Extended_figure_11c_right.pdf")
#wilcox.test(DGE_DDR_UP$gene_length,DGE_DDR_DOWN$gene_length,paired=F)
dotplot(Objectlist=list(DGE_DDR_UP$gene_length,DGE_DDR_DOWN$gene_length),col=c("red","blue"),
		labs=c("DEG up siCDK12","DEG down siCDK12"),pch=19,ylab="gene length(bp)",
		main="p=0.1322, wilcox test")
dev.off()

















