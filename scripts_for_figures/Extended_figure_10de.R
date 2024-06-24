##################################################################
#import data and libraries
##################################################################
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
data=read.table("RNASeq_U2OS_data.xls",header=TRUE)


if(!dir.exists("Extended_figure_10de")){
	dir.create("Extended_figure_10de")
}







#boxplot of log2FC (DEG UP and DOWN) showing gene length
gene_list=genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
#match genes in RNA-Seq and gene_list
pos=match(gene_list$gene_id,data$ID_ENTREZ)
elementMetadata(gene_list)=data[pos,]
gene_list$gene_length=width(gene_list)
#filter only significant and only genes longer than 10000 bp
gene_list_signif=gene_list[!is.na(gene_list$DGE_padj_siCdk12)& gene_list$DGE_padj_siCdk12<0.01]
gene_list_signif_UP=gene_list_signif[gene_list_signif$DGE_log2FoldChange_siCdk12>0.5]
gene_list_signif_DOWN=gene_list_signif[gene_list_signif$DGE_log2FoldChange_siCdk12<(-0.5)]
pdf("Extended_figure_10de/boxplot_DEG_gene_length.pdf")
boxplot(gene_list_signif_UP$gene_length,gene_list_signif_DOWN$gene_length,outline=F,varwidth=T,notch=T,
		col=c("red","blue"),ylab="gene length (bp)",
		names=c(paste("DEG up siCDK12 (",length(gene_list_signif_UP),")"),paste("DEG down siCDK12 (",length(gene_list_signif_DOWN),")")),
		main="p<2.2e-16, wilcox test")
dev.off()
#wilcox.test(gene_list_signif_UP$gene_length,gene_list_signif_DOWN$gene_length,paired=F)













#boxplot gene expression level between DEG UP/DOWN/NOdeg in ctrl condition
gene_list=genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
#match genes in RNA-Seq and gene_list
pos=match(gene_list$gene_id,data$ID_ENTREZ)
elementMetadata(gene_list)=data[pos,]
gene_list$gene_length=width(gene_list)
avg_ctrl_cpm=apply(as.matrix(elementMetadata(gene_list)[ ,grep("cpm_ctrl",colnames(elementMetadata(gene_list))) ]),1,mean)
gene_list$avg_ctrl_cpm=avg_ctrl_cpm
gene_list$RPKM=gene_list$avg_ctrl_cpm/gene_list$gene_length

gene_list_signif=gene_list[!is.na(gene_list$DGE_padj_siCdk12)& gene_list$DGE_padj_siCdk12<0.01]
gene_list_signif_UP=gene_list_signif[gene_list_signif$DGE_log2FoldChange_siCdk12>0.5,]
gene_list_signif_DOWN=gene_list_signif[gene_list_signif$DGE_log2FoldChange_siCdk12<(-0.5),]


pdf("Extended_figure_10de/boxplot_DEG_expressionLevel.pdf")
boxplot(gene_list_signif_UP$RPKM,gene_list_signif_DOWN$RPKM,outline=F,varwidth=T,notch=T,
		col=c("red","blue"),ylab="avg ctrl RPKM expression",
		names=c(paste("DEG up siCDK12 (",length(gene_list_signif_UP),")"),paste("DEG down siCDK12 (",length(gene_list_signif_DOWN),")")),
		main="p=1.234e-12, wilcox test")
dev.off()

#wilcox.test(gene_list_signif_UP$avg_ctrl_cpm,gene_list_signif_DOWN$avg_ctrl_cpm,paired=F)













