
########################################################################
#import data and set parameters
########################################################################
padjthresh=0.01
log2fcthresh=0.5
log2fcthreshNODEG=0.2
data=read.table("RNASeq_U2OS_data.xls",header=TRUE)

if(!dir.exists("Extended_figure_10a")){
	dir.create("Extended_figure_10a")
}






#find DEGs
DGE=data[,grepl("DGE_",colnames(data))|grepl("ID_",colnames(data))]

#identification of genes regulated by CDK12: DEGs in siCDK12 condition
pos_tokeep=!is.na(DGE$DGE_padj_siCdk12)
DGE=DGE[pos_tokeep,]
pos_signif=DGE$DGE_padj_siCdk12<padjthresh
pos_NO= DGE$DGE_log2FoldChange_siCdk12<log2fcthreshNODEG &DGE$DGE_log2FoldChange_siCdk12>(-log2fcthreshNODEG) &
		DGE$DGE_log2FoldChange_OHT<log2fcthreshNODEG &DGE$DGE_log2FoldChange_OHT>(-log2fcthreshNODEG) &
		DGE$DGE_log2FoldChange_siCdk12_OHT<log2fcthreshNODEG &DGE$DGE_log2FoldChange_siCdk12_OHT>(-log2fcthreshNODEG) 
DEG_NO=DGE[pos_NO,]
DGE=DGE[pos_signif,]
pos_UP=DGE$DGE_log2FoldChange_siCdk12>log2fcthresh
pos_DOWN=DGE$DGE_log2FoldChange_siCdk12<(-log2fcthresh)
DEG_UP=DGE[pos_UP,]
DEG_DOWN=DGE[pos_DOWN,]
DEG_UP=DEG_UP[order(-DEG_UP$DGE_log2FoldChange_siCdk12),]
DEG_DOWN=DEG_DOWN[order(-DEG_DOWN$DGE_log2FoldChange_siCdk12),]


print(paste("Number of siCDK12 DEG UP:",nrow(DEG_UP)))
print(paste("Number of siCDK12 DEG DOWN:",nrow(DEG_DOWN)))

#write DEG tables
write.table(DEG_UP,file=paste("Extended_figure_10a/DEG_UP_siCDK12_",log2fcthresh,".xls",sep=""),sep="\t",quote=F,row.names=T,col.names=NA)
write.table(DEG_DOWN,file=paste("Extended_figure_10a/DEG_DOWN_siCDK12_",log2fcthresh,".xls",sep=""),sep="\t",quote=F,row.names=T,col.names=NA)

#produce random sample of DEGs for Pol2 profile and ChroKit:
set.seed(123)
sample_UP_chrokit=as.character(sample(DEG_UP$ID_SYMBOL,500,replace=FALSE))
sample_DOWN_chrokit=as.character(sample(DEG_DOWN$ID_SYMBOL,500,replace=FALSE))
sample_NO_chrokit=as.character(sample(DEG_NO$ID_SYMBOL,500,replace=FALSE))
write.table(sample_UP_chrokit,file="Extended_figure_10a/sample_UP_chrokit.txt",sep="\t",quote=F,row.names=F,col.names=F)
write.table(sample_DOWN_chrokit,file="Extended_figure_10a/sample_DOWN_chrokit.txt",sep="\t",quote=F,row.names=F,col.names=F)
write.table(sample_NO_chrokit,file="Extended_figure_10a/sample_NO_chrokit.txt",sep="\t",quote=F,row.names=F,col.names=F)




#plot DEG log2FC
pdf("Extended_figure_10a/Extended_figure_10a.pdf",width=12,height=8)
par(mar=c(5,5,4,4))
plot(c(DEG_UP$DGE_log2FoldChange_siCdk12,DEG_DOWN$DGE_log2FoldChange_siCdk12),pch=20, ylab="log2 fold change",xlab="Ranked DEGs",
	col=c(rep("red",nrow(DEG_UP)),rep("blue",nrow(DEG_DOWN))),cex.axis=2, cex.lab=2.2,main=paste("log2FC thresh:",log2fcthresh,"padj thresh:",padjthresh))
legend("topright",legend=c(paste("DEG UP:",nrow(DEG_UP)),paste("DEG DOWN:",nrow(DEG_DOWN))),col=c("red","blue"),pch=20,cex=2)
dev.off()



