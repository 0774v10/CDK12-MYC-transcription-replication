###take DEGs after Myc OE (DEG up) and analyse how they behave in other conditions.

########################################################################
#imports and parameters
########################################################################

library(inline)
source("../_functions.R")
padjthresh=0.01
log2fcthresh=0.5
data=read.table("RNASeq_U2OS_data.xls",header=TRUE)

if(!dir.exists("Extended_figure_13")){
	dir.create("Extended_figure_13")
}


DGE=data[,grepl("DGE_",colnames(data))|grepl("ID_",colnames(data))]

#identification of genes regulated by Myc: DEGs in OHT condition
pos_tokeep=!is.na(DGE$DGE_padj_OHT)
DGE=DGE[pos_tokeep,]
pos_signif=DGE$DGE_padj_OHT<padjthresh
DGE=DGE[pos_signif,]
pos_UP=DGE$DGE_log2FoldChange_OHT>log2fcthresh
pos_DOWN=DGE$DGE_log2FoldChange_OHT<(-log2fcthresh)
DEG_UP=DGE[pos_UP,]
DEG_DOWN=DGE[pos_DOWN,]
DEG_UP=DEG_UP[order(-DEG_UP$DGE_log2FoldChange_OHT),]
DEG_DOWN=DEG_DOWN[order(-DEG_DOWN$DGE_log2FoldChange_OHT),]

print(paste("Number of DEG UP after Myc overexpression:",nrow(DEG_UP)))
print(paste("Number of DEG DOWN after Myc overexpression:",nrow(DEG_DOWN)))
#write DEG tables
write.table(DEG_UP,file=paste("Extended_figure_13/DEG_UP_OHT_",log2fcthresh,".xls",sep=""),sep="\t",quote=F,row.names=T,col.names=NA)
write.table(DEG_DOWN,file=paste("Extended_figure_13/DEG_DOWN_OHT_",log2fcthresh,".xls",sep=""),sep="\t",quote=F,row.names=T,col.names=NA)


#produce random sample of DEGs for Pol2 profile and ChroKit:
set.seed(123)
sample_UP_chrokit=as.character(sample(DEG_UP$ID_SYMBOL,500,replace=FALSE))
write.table(sample_UP_chrokit,file="Extended_figure_13/sample_UP_chrokit.txt",sep="\t",quote=F,row.names=F,col.names=F)



pos_NA=is.na(DEG_UP$ID_SYMBOL)
DEG_UP=DEG_UP[which(!pos_NA),]
rownames(DEG_UP)=DEG_UP$ID_SYMBOL

#heatmap log2FC DEGs UP OHT:
part_siCDK12=DEG_UP[,grepl("_siCdk12$",colnames(DEG_UP))][,c(2,6)]
colnames(part_siCDK12)=c("log2FoldChange","padj")
part_siCDK12_OHT=DEG_UP[,grepl("_siCdk12_OHT$",colnames(DEG_UP))][,c(2,6)]
colnames(part_siCDK12_OHT)=c("log2FoldChange","padj")
part_OHT=DEG_UP[,grepl("_OHT$",colnames(DEG_UP))& !grepl("_siCdk12",colnames(DEG_UP))][,c(2,6)]
colnames(part_OHT)=c("log2FoldChange","padj")
listdegs=list(part_siCDK12,part_OHT,part_siCDK12_OHT)
names(listdegs)=c("siCDK12","OHT","siCDK12_OHT")


#heatmap log2FC of Myc deg UP (and other conditions)
RNA_result_obj=matrixFromRNAresults(resultsList=listdegs,padj_thresh=padjthresh,log2FCthresh=log2fcthresh,orderingDriver=2,orderingtype="ranking")
plotRNAresults(matrixResults=RNA_result_obj,quantile_saturation=0.98,display_significance="lateral",
				palette=c("darkgreen","grey80","darkred"),signif_thresh=padjthresh,breaks=10,fileName="Extended_figure_13/Extended_figure_13a.pdf",rownames_thresh=50)


#boxplot of the ranked values (as if it was "clsuter 1")
pdf(paste("Extended_figure_13/Extended_figure_13b.pdf",sep=""))
par(mar=c(11,5,4,4))
boxplot(RNA_result_obj$mat,ylab="log2 fold change",las=2,notch=TRUE,main=paste("MYC overexpressed"))
abline(h=0,lty=2,col="red")
dev.off()

