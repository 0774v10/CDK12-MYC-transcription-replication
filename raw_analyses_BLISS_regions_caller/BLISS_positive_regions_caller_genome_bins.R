
#load data and functions
library(inline)
library(GenomicRanges)
library(rtracklayer)
source("../_functions.R")
hg19_2k=read.table("hg19_2k.bed")
divideinto=10



#Here put the bigWig files for each BLISS sample:
listwigs=list(
	"file1.bw",
	"file2.bw",
	"file3.bw",
	"file4.bw",
	"file5.bw",
	"file6.bw",
	"file7.bw",
	"file8.bw"
)
names(listwigs)=c(
	"namefile1",
	"namefile2",
	"namefile3",
	"namefile4",
	"namefile5",
	"namefile6",
	"namefile7",
	"namefile8"
)



#get genomic ranges of binned hg19 genome assembly. Choose 2kb window binning, sice each peak seems to be 50k wide
#(from genome browser)
hg19_2k_range=GRanges(Rle(hg19_2k$V1),IRanges(hg19_2k$V2,hg19_2k$V3))
steps=ceiling(length(hg19_2k_range)/divideinto)
starts=seq(1,steps*divideinto,steps)
ends=starts+steps-1
ends[length(ends)]=length(hg19_2k_range)





df_2k_genome_bins=data.frame(chr=hg19_2k$V1,start=hg19_2k$V2,end=hg19_2k$V3)
mat_toadd=as.data.frame(matrix(rep(NA,length(listwigs)*length(hg19_2k_range)),ncol=length(listwigs)))
colnames(mat_toadd)=names(listwigs)
df_2k_genome_bins=cbind(df_2k_genome_bins,mat_toadd)
for (i in 1:length(listwigs)){
	#for this wig, calculate base enrichment splitted for memory
	currentenrich=listwigs[i]
	currentwig=listwigs[[i]]
	print (paste("Processing wig",names(currentenrich)))
	listfrags=list()
	for (k in 1:divideinto){
		print(paste("......Processing bin",k))
		currentbin=hg19_2k_range[starts[k]:ends[k]]
		basecov=GRbaseCoverage2(Object=currentbin, signalfile=currentwig,signalfileNorm=currentwig,signalControl=NULL,signalControlSpike=NULL,multiplFactor=1e+06)
		result=makeMatrixFrombaseCoverage(GRbaseCoverageOutput=basecov[[1]],Nbins=1,Snorm=FALSE,norm_factor=basecov[[2]])
		listfrags=c(listfrags,result)
	}
	currentcov=do.call("rbind",listfrags)
	df_2k_genome_bins[,i+3]=currentcov
}
colnames(df_2k_genome_bins)[4:ncol(df_2k_genome_bins)]=paste("cov_",colnames(df_2k_genome_bins)[4:ncol(df_2k_genome_bins)],sep="")
write.table(df_2k_genome_bins,file="df_2k_genome_BLISS.xls",sep="\t",quote=F,row.names=F,col.names=T)


#evaluation of densities of enrichments and threshold 
if (!dir.exists("distribution_enrichment_bins")){
	dir.create("distribution_enrichment_bins")
}




#plot frequency of enrichment (try to find a threshold above which call "positive" or "negative" regions)
mat_2k=df_2k_genome_bins[,grepl("cov_",colnames(df_2k_genome_bins))]
inf=quantile(log2(as.matrix(mat_2k[mat_2k!=0])),.001)
sup=quantile(log2(as.matrix(mat_2k)),.999)


#here plot the distribution of the coverage of the various 2k bins.
#adjust this based on the number of the samples in the matrix
pdf("distribution_enrichment_bins/Frequency_enrichments_2k_bin_hg19_BLISS.pdf")
plot(density(log2(mat_2k[,1]),bw=0.05),col="black",xlim=c(inf,sup),ylim=c(0,1.2),lwd=3,xlab="log2 bin coverage")
lines(density(log2(mat_2k[,2]),bw=0.05),col="grey50",lwd=3)
lines(density(log2(mat_2k[,3]),bw=0.05),col="blue4",lwd=3)
lines(density(log2(mat_2k[,4]),bw=0.05),col="blue",lwd=3)
lines(density(log2(mat_2k[,5]),bw=0.05),col="green4",lwd=3)
lines(density(log2(mat_2k[,6]),bw=0.05),col="green",lwd=3)
lines(density(log2(mat_2k[,7]),bw=0.05),col="red4",lwd=3)
lines(density(log2(mat_2k[,8]),bw=0.05),col="red",lwd=3)
abline(v=4.6,col="red",lty=2)
legend("topleft",legend=colnames(mat_2k),col=c("black","grey50","blue4","blue","green4","green","red4","red"),lwd=3)
dev.off()







#Use the set tresholds to define regions with "significant number of UMI BLISS events".
#This threshold is decided based on the previous plot.
#use 2k window bins hg19. A rudimental peak caller.
thresh=4.6
maxgap=1
thresh_val=2^thresh
df_2k_genome_bins=read.table("df_2k_genome_BLISS.xls",header=T)
toprocess=grep("cov_",colnames(df_2k_genome_bins),value=T)

#split chromosome
splitted_binned=split(df_2k_genome_bins,df_2k_genome_bins$chr)


peaksforbams=rep(list(data.frame(chr=character(),start=character(),end=character())),length(toprocess))
names(peaksforbams)=toprocess
for ( j in 1:length(splitted_binned)){
	print(j)
	current=splitted_binned[[j]]
	currentchr=names(splitted_binned)[j]

	for(i in 1:length(toprocess)){
		currentcol=current[,toprocess[i]]
		positives=currentcol>thresh_val
		positives[positives]=1

		#join together positive regions with gap<=maxgap
		#https://masterr.org/r/how-to-find-consecutive-repeats-in-r/
		sequence=rle(positives)
		vals=sequence$values
		reps=sequence$lengths
		tocollapse=reps<=maxgap & vals==0
		expandedcumsum=cumsum(reps)
		endpos=expandedcumsum[tocollapse]
		startpos=(expandedcumsum-reps+1)[tocollapse]

		if(length(endpos)>0 & length(startpos)>0){
			positives_corrected=positives
			#make 0s from startpos to endpos as 1s
			for(k in 1:length(startpos)){
				fromstart=startpos[k]
				toend=endpos[k]
				rang=fromstart:toend
				positives_corrected[rang]=1
			}

			#contiguous 1s for start/end range
			sequence_peaks=rle(positives_corrected)
			vals_peaks=sequence_peaks$values
			reps_peaks=sequence_peaks$lengths
			#find all windows with 1
			expandedcumsum_peaks=cumsum(reps_peaks)
			tocollapse_peaks=vals_peaks==1
			endpos_peaks=expandedcumsum_peaks[tocollapse_peaks]
			startpos_peaks=(expandedcumsum_peaks-reps_peaks+1)[tocollapse_peaks]


			start_genomic_pos=current[,2][startpos_peaks]
			end_genomic_pos=current[,3][endpos_peaks]
			df=data.frame(chr=rep(currentchr,length(start_genomic_pos)),start=start_genomic_pos,end=end_genomic_pos)
			peaksforbams[[i]]=rbind(peaksforbams[[i]],df)			
		}
		

	}

}

#print peaks to files:
for(i in 1:length(peaksforbams)){
	currentname=names(peaksforbams)[i]
	filename=paste("ROI_BLISS_positive_",currentname,".bed",sep="")
	genomiclocations=peaksforbams[[i]]
	write.table(genomiclocations,file=paste("distribution_enrichment_bins/",filename,sep=""),row.names=F,col.names=F,sep="\t",quote=F )

}

