library(org.Hs.eg.db)
library(DESeq2)



#get raw counts for each gene:
U2OS_siCtrl_1_S13671_count=read.table("counts_rnaseq/3052.count",header=TRUE)[,7] 	
U2OS_siCtrl_2_S13672_count=read.table("counts_rnaseq/3051.count",header=TRUE)[,7] 	
U2OS_siCtrl_3_S13673_count=read.table("counts_rnaseq/3050.count",header=TRUE)[,7] 
U2OS_siCdk12_1_1_S13676_count=read.table("counts_rnaseq/3047.count",header=TRUE)[,7] 	
U2OS_siCdk12_1_2_S13677_count=read.table("counts_rnaseq/3046.count",header=TRUE)[,7] 	
U2OS_siCdk12_2_1_S13678_count=read.table("counts_rnaseq/3045.count",header=TRUE)[,7] 	
U2OS_siCdk12_2_2_S13679_count=read.table("counts_rnaseq/3044.count",header=TRUE)[,7] 
U2OS_siCtrl_1_OHT_S13686_count=read.table("counts_rnaseq/3037.count",header=TRUE)[,7] 	
U2OS_siCtrl_2_OHT_S13687_count=read.table("counts_rnaseq/3036.count",header=TRUE)[,7] 	
U2OS_siCtrl_3_OHT_S13688_count=read.table("counts_rnaseq/3035.count",header=TRUE)[,7] 
U2OS_siCdk12_1_1_OHT_S13691_count=read.table("counts_rnaseq/3032.count",header=TRUE)[,7] 	
U2OS_siCdk12_1_2_OHT_S13692_count=read.table("counts_rnaseq/3031.count",header=TRUE)[,7] 	
U2OS_siCdk12_2_1_OHT_S13693_count=read.table("counts_rnaseq/3030.count",header=TRUE)[,7] 	
U2OS_siCdk12_2_2_OHT_S13694_count=read.table("counts_rnaseq/3029.count",header=TRUE)[,7] 





df_counts=cbind(
U2OS_siCtrl_1_S13671_count,U2OS_siCtrl_2_S13672_count,U2OS_siCtrl_3_S13673_count,U2OS_siCdk12_1_1_S13676_count,U2OS_siCdk12_1_2_S13677_count,U2OS_siCdk12_2_1_S13678_count,U2OS_siCdk12_2_2_S13679_count,U2OS_siCtrl_1_OHT_S13686_count,U2OS_siCtrl_2_OHT_S13687_count,U2OS_siCtrl_3_OHT_S13688_count,U2OS_siCdk12_1_1_OHT_S13691_count,U2OS_siCdk12_1_2_OHT_S13692_count,U2OS_siCdk12_2_1_OHT_S13693_count,U2OS_siCdk12_2_2_OHT_S13694_count
)
colnames(df_counts)=c(
"counts_ctrl_1","counts_ctrl_2","counts_ctrl_3","counts_siCdk12_1_1","counts_siCdk12_1_2","counts_siCdk12_2_1","counts_siCdk12_2_2","counts_OHT_1","counts_OHT_2","counts_OHT_3","counts_siCdk12_1_1_OHT","counts_siCdk12_1_2_OHT","counts_siCdk12_2_1_OHT","counts_siCdk12_2_2_OHT"
)
gID=read.table("counts_rnaseq/3052.count",header=TRUE)$Geneid
row.names(df_counts)=gID

cpm_matrix=apply(df_counts,2,function(i){ (i/sum(i))*1000000 })


##############
#DESeq2 analyses
##############

#OHT
df_rawcounts=df_counts[,c(1:3,8:10)]
#build design matrix
df_design=data.frame(condition=as.factor(c(rep("a_ctrl",3),rep("b_OHT",3))), sample=as.factor(colnames(df_rawcounts)))
rownames(df_design)=colnames(df_rawcounts)
#call DEGs for paired samples
dds <- DESeqDataSetFromMatrix(countData = df_rawcounts,
                              colData = df_design,
                              design = ~ condition)
dds <- DESeq(dds)
DEG_OHT <- results(dds)


#siCDK12
df_rawcounts=df_counts[,c(1:3,4:7)]
#build design matrix
df_design=data.frame(condition=as.factor(c(rep("a_ctrl",3),rep("b_siCDk12",4))), sample=as.factor(colnames(df_rawcounts)))
rownames(df_design)=colnames(df_rawcounts)
#call DEGs for paired samples
dds <- DESeqDataSetFromMatrix(countData = df_rawcounts,
                              colData = df_design,
                              design = ~ condition)
dds <- DESeq(dds)
DEG_siCdk12 <- results(dds)


#siCDK12/OHT
df_rawcounts=df_counts[,c(1:3,11:14)]
#build design matrix
df_design=data.frame(condition=as.factor(c(rep("a_ctrl",3),rep("b_siCDk12OHT",4))), sample=as.factor(colnames(df_rawcounts)))
rownames(df_design)=colnames(df_rawcounts)
#call DEGs for paired samples
dds <- DESeqDataSetFromMatrix(countData = df_rawcounts,
                              colData = df_design,
                              design = ~ condition)
dds <- DESeq(dds)
DEG_siCdk12_OHT <- results(dds)





##############
#IDs conversion
##############
entrez=strsplit(rownames(df_counts),split="\\.")
entrez=sapply(entrez,"[[",1)
symbol <- mapIds(org.Hs.eg.db,
                     keys=entrez,
                     column="SYMBOL",
                     keytype="ENTREZID",
                     multiVals="first")

refseq <- mapIds(org.Hs.eg.db,
                     keys=entrez,
                     column="REFSEQ",
                     keytype="ENTREZID",
                     multiVals="first")

ensembl<-mapIds(org.Hs.eg.db,
                     keys=entrez,
                     column="ENSEMBL",
                     keytype="ENTREZID",
                     multiVals="first")

#join the results in a single table
finaldf=data.frame(ID_ENTREZ=rownames(df_counts),ID_ENSEMBL=ensembl,ID_SYMBOL=symbol,ID_REFSEQ=refseq)






 	


##############
#join tables together
##############
colnames(DEG_OHT)=paste("DGE_",colnames(DEG_OHT),"_OHT",sep="")
colnames(DEG_siCdk12)=paste("DGE_",colnames(DEG_siCdk12),"_siCdk12",sep="")
colnames(DEG_siCdk12_OHT)=paste("DGE_",colnames(DEG_siCdk12_OHT),"_siCdk12_OHT",sep="")
DEG_Cdk12KDOHT=cbind(DEG_OHT,DEG_siCdk12)
DEG_Cdk12KDOHT=cbind(DEG_Cdk12KDOHT,DEG_siCdk12_OHT)
data_RNASeq=cbind(finaldf,DEG_Cdk12KDOHT)
data_RNASeq=cbind(data_RNASeq,df_counts)

#append CPM to the final table
colnames(cpm_matrix)=gsub("counts_","cpm_",colnames(cpm_matrix))
data_RNASeq=cbind(data_RNASeq,cpm_matrix)

#save pre-processed data
write.table(data_RNASeq,"RNASeq_U2OS_data.xls",quote=FALSE,row.names=F,col.names=TRUE,sep="\t")
