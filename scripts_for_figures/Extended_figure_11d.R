####################################################################
#import data and source functions.
####################################################################

library(Rcpp)
library(inline)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#import other libraries
library(Rsamtools)
library(qs)
library(bamsignals)

source("../_functions.R")
data=read.table("RNASeq_U2OS_data.xls",header=TRUE)


#here put the path to the bam files of the RNA-Seq for each sample.
#for example, ctrlRNAseq_rep1_bam="path/to/file/ctrlRNAseq_rep1.bam" 
ctrlRNAseq_rep1_bam="" 
ctrlRNAseq_rep2_bam="" 
ctrlRNAseq_rep3_bam="" 
siCdk12RNAseq_rep1_bam="" 
siCdk12RNAseq_rep2_bam="" 
siCdk12RNAseq_rep3_bam="" 
siCdk12RNAseq_rep4_bam=""

#long DDR genes: BRCA1, BARD1, BLM, PALB2, FAN1, RAD51,ATR, ATM
longDDRgenes=c("BRCA1", "BARD1", "BLM", "PALB2", "FAN1", "RAD51","ATR", "ATM")


if(!dir.exists("Extended_figure_11d")){
	dir.create("Extended_figure_11d")
}







####################################################################
#compare RNA-Seq enrichment between NT and shCDK12 for each exon of the
#long DDR genes
####################################################################

#extract exons form hg19 TxDB, listed by their ENTREZ ID
exon_list=exonsBy(TxDb.Hsapiens.UCSC.hg19.knownGene,by=c("gene"))
#map gene symbols desired with RNA-Seq table and, with ENTREZ ID, take exons per gene
pos=match(longDDRgenes,data$ID_SYMBOL)
ENTREZ_longDDR=data$ID_ENTREZ[pos]
exons_longDDR=exon_list[as.character(ENTREZ_longDDR)]
exons_longDDR=lapply(1:length(exons_longDDR),function(i){
	    current=exons_longDDR[[i]]
		current$gene=longDDRgenes[i]
		current=as.data.frame(current)
		return(current)
})
exons_longDDR=do.call(rbind,exons_longDDR)
exons_longDDR=makeGRangesFromDataFrame(exons_longDDR,keep.extra.columns=T)



#now calculate enrichment in exons from RNA-Seq bam files.
cov_ctrl_1=GRbaseCoverage2(Object=exons_longDDR,signalfile=ctrlRNAseq_rep1_bam,signalfileNorm=ctrlRNAseq_rep1_bam)
cov_ctrl_2=GRbaseCoverage2(Object=exons_longDDR,signalfile=ctrlRNAseq_rep2_bam,signalfileNorm=ctrlRNAseq_rep2_bam)
cov_ctrl_3=GRbaseCoverage2(Object=exons_longDDR,signalfile=ctrlRNAseq_rep3_bam,signalfileNorm=ctrlRNAseq_rep3_bam)
cov_siCdk12_1=GRbaseCoverage2(Object=exons_longDDR,signalfile=siCdk12RNAseq_rep1_bam,signalfileNorm=siCdk12RNAseq_rep1_bam)
cov_siCdk12_2=GRbaseCoverage2(Object=exons_longDDR,signalfile=siCdk12RNAseq_rep2_bam,signalfileNorm=siCdk12RNAseq_rep2_bam)
cov_siCdk12_3=GRbaseCoverage2(Object=exons_longDDR,signalfile=siCdk12RNAseq_rep3_bam,signalfileNorm=siCdk12RNAseq_rep3_bam)
cov_siCdk12_4=GRbaseCoverage2(Object=exons_longDDR,signalfile=siCdk12RNAseq_rep4_bam,signalfileNorm=siCdk12RNAseq_rep4_bam)

cov_ctrl_1_tot=makeMatrixFrombaseCoverage(GRbaseCoverageOutput=cov_ctrl_1[[1]],Nbins=1,Snorm=TRUE,key=cov_ctrl_1[[2]],norm_factor=cov_ctrl_1[[3]])
cov_ctrl_2_tot=makeMatrixFrombaseCoverage(GRbaseCoverageOutput=cov_ctrl_2[[1]],Nbins=1,Snorm=TRUE,key=cov_ctrl_2[[2]],norm_factor=cov_ctrl_2[[3]])
cov_ctrl_3_tot=makeMatrixFrombaseCoverage(GRbaseCoverageOutput=cov_ctrl_3[[1]],Nbins=1,Snorm=TRUE,key=cov_ctrl_3[[2]],norm_factor=cov_ctrl_3[[3]])
cov_siCdk12_1_tot=makeMatrixFrombaseCoverage(GRbaseCoverageOutput=cov_siCdk12_1[[1]],Nbins=1,Snorm=TRUE,key=cov_siCdk12_1[[2]],norm_factor=cov_siCdk12_1[[3]])
cov_siCdk12_2_tot=makeMatrixFrombaseCoverage(GRbaseCoverageOutput=cov_siCdk12_2[[1]],Nbins=1,Snorm=TRUE,key=cov_siCdk12_2[[2]],norm_factor=cov_siCdk12_2[[3]])
cov_siCdk12_3_tot=makeMatrixFrombaseCoverage(GRbaseCoverageOutput=cov_siCdk12_3[[1]],Nbins=1,Snorm=TRUE,key=cov_siCdk12_3[[2]],norm_factor=cov_siCdk12_3[[3]])
cov_siCdk12_4_tot=makeMatrixFrombaseCoverage(GRbaseCoverageOutput=cov_siCdk12_4[[1]],Nbins=1,Snorm=TRUE,key=cov_siCdk12_4[[2]],norm_factor=cov_siCdk12_4[[3]])

exons_longDDR$cov_ctrl_1=cov_ctrl_1_tot
exons_longDDR$cov_ctrl_2=cov_ctrl_2_tot
exons_longDDR$cov_ctrl_3=cov_ctrl_3_tot
exons_longDDR$cov_siCdk12_1=cov_siCdk12_1_tot
exons_longDDR$cov_siCdk12_2=cov_siCdk12_2_tot
exons_longDDR$cov_siCdk12_3=cov_siCdk12_3_tot
exons_longDDR$cov_siCdk12_4=cov_siCdk12_4_tot

#convert back to data frame
exons_longDDR_df=as.data.frame(exons_longDDR)
#split
exons_longDDR_df_split=split(exons_longDDR_df,exons_longDDR_df$gene)

for(i in 1:length(exons_longDDR_df_split)){
	current=exons_longDDR_df_split[[i]]
	mat=current[,9:ncol(current)]

	matmeans=sapply(1:nrow(as.matrix(mat)),function(k){
		currentline=as.matrix(mat)[k,]
		type=factor(as.factor(c(1,1,1,2,2,2,2)),labels=c("ctrl","sicdk12"))
		tapply(currentline,type,mean)
	})
	finalmat=cbind(mat,t(matmeans))
	if(unique(current$strand=="-")){
		finalmat=finalmat[nrow(finalmat):1,]
	}
	pdf(paste0("Extended_figure_11d/Plot_enrichment_exons_",names(exons_longDDR_df_split)[i],".pdf"))
	matplot(as.data.frame(finalmat),type="l",col=c(rep("red",3),rep("blue",4),"red","blue"),lty=c(rep(2,7),rep(1,2)),lwd=c(rep(1,7),rep(4,2)),ylab="RNA-Seq signal exons (rpm/bp)",xlab="exons",main=names(exons_longDDR_df_split)[i])
	legend("bottom",legend=c("ctrl","siCdk12","mean ctrl","mean siCDK12"),col=c("red","blue"),lty=c(2,2,1,1),lwd=c(1,1,4,4))
	dev.off()

	stat_mat=lapply(1:nrow(mat),function(k){
		pval=wilcox.test( as.numeric(mat[k,1:3]),as.numeric(mat[k,4:7]),paired=F,exact=T)$p.value
		diffpos=(mean(as.numeric((mat[k,4:7])))-mean(as.numeric(mat[k,1:3]) ))>0
		lll=list(pval,diffpos)
		names(lll)=c("pval","diffpos")
		return(lll)
	})
	pvals=sapply(stat_mat,"[[",1)
	diffpos=sapply(stat_mat,"[[",2)
	if(unique(current$strand=="-")){
		pvals=pvals[length(pvals):1]
		diffpos=diffpos[length(diffpos):1]
	}
	pdf(paste0("Extended_figure_11d/Plot_pval_exons_",names(exons_longDDR_df_split)[i],".pdf"))
	plot(pvals,col=ifelse(diffpos,"blue","red"),pch=19)
	abline(h=0.05,lty=2,col="green")
	legend("topright",legend=c("more transcription in ctrl","more transcription in siCDK12"),
				col=c("red","blue"),pch=19)
	dev.off()
}














