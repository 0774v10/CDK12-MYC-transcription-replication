#useful functions for NGS data analysis.
#if you use them, please aknowledge Ottavio Croci


#plot a dotplot
dotplot<-function(Objectlist,col,labs,widthlines=0.2,center="median",...){
	stripchart(Objectlist,vertical=TRUE,method="jitter",xaxt="n",...)
	axis(1,at=1:length(Objectlist),labels=labs)
	for(i in 1:length(Objectlist)){
		segments(x0=i-widthlines,x1=i+widthlines,y0=median(Objectlist[[i]]),y1=median(Objectlist[[i]]),col="black")
	}
}



#cluster a list of matrixes
clusterMatrixKmeans<-function(matlist,clustinds,numberclusters,startingpoints,iter){
  #matlist: list of matrixes for which the rows must be clustered
  #clustinds: indexes rpresenting the position of the matrixes inside matlist 
  #         that will drive the clustering. All the other matrixes will be reordered accordingly
  #numberclusters: number of clusters for K-means
  #startingpoints: number of starting points for the kmeans clustering
  #iter: number of iteration for k-mean clustering

  #k-means is not deterministic! make it reproducible across different runs
  set.seed(123)
  if (!is.null(clustinds)){
    if(class(clustinds)!="integer"){
      stop("'clustinds' must be an integer or NULL...")
    }
  }
  if(class(matlist)!="list"){
    stop("'matlist' must be a list...")
  }
  for(i in 1:length(matlist)){
    if(class(matlist[[i]])!="matrix"){
      stop("each element of 'matlist' must be a matrix...")
    }
  }
  if (class(numberclusters)!="numeric" & class(numberclusters)!="integer"){
    stop("'numberclusters' must be a number...")
  }
  if(class(startingpoints)!="numeric" & class(startingpoints)!="integer"){
    stop("'startingpoints' must be a number...")
  }else{
    if(startingpoints<=0){
      stop("'startingpoints' must be > 0...")
    }
  }
  if(class(iter)!="numeric" & class(iter)!="integer"){
    stop("'iter' must be a number...")
  }else{
    if(iter<=0){
      stop("'iter' must be > 0...")
    }
  }

  #initialize the big matrx from all the matrixes
  mat <- NULL
  for (i in 1:length(matlist)) {
    mat <- cbind(mat, matlist[[i]])
  }

  #if indexes of clustering are not null,
  if(!is.null(clustinds) & numberclusters>0 & numberclusters<434){
    #prepare the matrix composed only by the matrixes that will drive the clustering
    clmat = NULL      
    for (i in clustinds){
      clmat <- cbind(clmat, matlist[[i]])
    } 
    #treat NAs (if all NAs, return NULL)
    NAcounts <- apply(mat, 1, function(x) length(which(is.na(x))))
    NAinds <- which(NAcounts > ncol(mat) * 0.3)
    if (length(NAinds) == nrow(mat)) 
        return(NULL)
    if (length(NAinds) > 0) {
        mat <- mat[-NAinds, ]
        clmat <- clmat[-NAinds, ]
    }

    #check clmat: no more clusters than 2^n, where n is all the possible combinations 
    #of the rows of the matrix
    maxclusterallowed=table(!duplicated(clmat))["TRUE"]
    if(numberclusters > maxclusterallowed ){
      numberclusters=maxclusterallowed
      print(paste("warning: max clustering allowed is",maxclusterallowed))
    }
    #HERE modify or add kmeans parameters, or alternatively set them by the user from UI
    #if very fast and you want more precise, increase iter.max (10 is default)
    #nstart is the initial random configuration (default of original function=1)
    clustobject <- kmeans(clmat, centers = numberclusters,nstart=startingpoints,iter.max=iter)
    groups=clustobject$centers
    clusters= clustobject$cluster
    ord=order(clusters)
    #you can improve ordering/separation of clusters between them according to groups matrix
    #keep more sparate those more different
    #order the matrix
    mat =mat[ord,]
  
  #otherwise cluster object is null and the order is not given (no clustering, no ranking)
  }else{
    ord=1:nrow(mat)
    clustobject=NULL
  }
  resul=list(mat,ord,clustobject)
  names(resul)=c("mat","ord","clustobject")
  return(resul)
}






########################################################################################
########################################################################################
#GeneSet enrichment analysis
########################################################################################
########################################################################################

#from simple GSEA analysis (investigate gene sets) returns the
#array of q values of the first "top" hits
extractqvalfromGSEAsimple<-function(fileToRead,top=5) {

	if(!file.exists(fileToRead)){
		stop("File doesn't exist")
	}
	if(top>20){
		top=20
	}
	x=readLines(fileToRead)
	#be careful: start with the 11 th line
	linesok=x[10:(10+top-1)]
	linesok=strsplit(linesok,split="\t")
	qvals=as.numeric(unlist(lapply(linesok,"[[",7)))
	names(qvals)=unlist(lapply(linesok,"[[",1))
	return(qvals)
}



processGSEA<-function(conditions,set,dire,pthr=0.05,NESthr=2,topthr=10) {
	#from result of GSEA java app, do a heatmap of the results
	#BE CAREFUL: nomenclature must be: xxxx_condition_setxxxx and set/conditions must be 
	#exactly the same as those given in input. Do it in the dedicated directory.
	#DO NOT put files with similar names of the GSEA results in the same directory
	#it consider both negative and positive gene sets

	#conditions=name of the conditions to put together. must be exactly in the file names
	#set= gene set used, that must correspond to the file name
	#dire=directory of GSEA results
	#pthr=threshold of the FWER p value for the results
	#NESthr=threshold for Normalized enrichment score
	#topthr=only top xxx (according to NES in the results) will be considered
	dircontent=dir(dire)
	if (length(dircontent)==0){
		stop("Nothing in directory...")
	}	
	allROIS=list()
	alltables=list()
	for (k in 1:length(conditions)){
		print (conditions[k])
		dirname=dircontent[grepl(paste(".+_",conditions[k],"_",set,".+",sep=""),dircontent)]
		dir2content=dir(paste(dire,"/",dirname,sep=""))
		if (length(dir2content)==0){
			stop("Nothing in subdirectory...")
		}
		filenameneg=paste(dire,"/",dirname,"/",dir2content[grepl("gsea_report_for_.+neg.+xls",dir2content)],sep="")
		filenamepos=paste(dire,"/",dirname,"/",dir2content[grepl("gsea_report_for_.+pos.+xls",dir2content)],sep="")
		negtable=read.table(filenameneg,sep="\t",header=TRUE)
		postable=read.table(filenamepos,sep="\t",header=TRUE)

		#save entire tables:
		intermtable=rbind(negtable,postable)
		intermtable=intermtable[,c(1,6,9)]
		alltables[[k]]=intermtable

		#filter tables according to FWER.p.val and rbind them. Then take only names and NES (or pval/padj)
		negtable=negtable[negtable$FWER.p.val<pthr,]
		postable=postable[postable$FWER.p.val<pthr,]
		#put rank threshold (ex. only top 5)
		negtable=negtable[negtable$NES< (-NESthr),]
		postable=postable[postable$NES> NESthr,]
		#max top n from GSEA results (NES ranked)
		if (nrow(negtable)>topthr){
			negtable=negtable[1:topthr,]
		}
		if (nrow(postable)>topthr){
			postable=postable[1:topthr,]
		}
		
		finaltable=rbind(negtable,postable)
		finaltable=finaltable[,c(1,6,9)]
		allROIS[[k]]=finaltable
	}
	#now union of all geneset for that category. For each condition, populate the final matrix
	allgenesets=c()
	for (k in allROIS){
		allgenesets=c(allgenesets,as.character(k$NAME))
	}
	allgenesets=unique(allgenesets)
	#For each condition, populate the final matrix
	finalmatrix=matrix(rep(NA,length(conditions)*length(allgenesets)),nrow=length(allgenesets))
	finalmatrix_pvals=finalmatrix
	rownames(finalmatrix)=rownames(finalmatrix_pvals)=allgenesets
	colnames(finalmatrix)=colnames(finalmatrix_pvals)=conditions


	
	#extract Normalized enrichment scores for all the conditions (columns)
	#for each column
	for (k in 1:ncol(finalmatrix)){
		#for each geneset (rows)
		values=sapply(1:nrow(finalmatrix), function(j) {
			x=alltables[[k]]
			n=rownames(finalmatrix)[j]
			if (n %in% as.character(x$NAME)){
				return(  x[as.character(x$NAME)==n,]$NES )
			}else{
				return(0)
			}
			#return array of the NES, padj when found
		})
		finalmatrix[,k]=values
	}

	#extract FWER- pvalues for all the conditions (columns)
	for (k in 1:ncol(finalmatrix_pvals)){
		#for each geneset (rows)
		values=sapply(1:nrow(finalmatrix_pvals), function(j) {
			x=alltables[[k]]
			n=rownames(finalmatrix_pvals)[j]
			if (n %in% as.character(x$NAME)){
				return(  x[as.character(x$NAME)==n,]$FWER.p.val )
			}else{
				#if not found, pvalue is not significant, therefore ==1
				return(1)
			}
			#return array of the NES, padj when found
		})
		finalmatrix_pvals[,k]=values
	}

	return(list(NES=finalmatrix,pvals=finalmatrix_pvals) )

	# mabs=max(abs(max(finalmatrix)),abs(min(finalmatrix)))
	# brk=c( seq( -mabs , mabs,mabs/100))
	# my_palette2 <- colorRampPalette(c("green","black","red"))(n = length(brk)-1 )
	# par(mar=c(6,6,6,10),oma=c(8,8,8,18))
	# heatmap.2(finalmatrix,breaks=brk,col=my_palette2,Colv=colcluster,Rowv=rowcluster)
}





#plot the results (NES/pval)
gseadotplot<-function(gseaObj,top=20,padj_thresh=0.05,NES_thresh=0) {
	#gseaObj: gsea object produced by clusterProfiler package, with GSEA or other functions
	#top: how many top genes (abs value of NES) to plot
	#padj_thresh: threshold of padjusted
	#NES_thresh: threshold of NES (absolute value)
	library(clusterProfiler)
	library(forcats)
	library(ggplot2)
	library(DOSE)
	if(class(top)!="numeric"){
		stop("'top' must be a number, suggested ~10/20")
	}
	if(top<1){
		stop("'top' must be at least 1")
	}
	if(class(padj_thresh)!="numeric" | padj_thresh<0 | padj_thresh>1){
		stop("'padj_thresh' must be a number betwen 0 and 1 (suggested 0.05)")
	}
	if(class(NES_thresh)!="numeric" | NES_thresh<0 ){
		stop("'NES_thresh' must be a positive number")
	}
	if(class(gseaObj)!="gseaResult"){
		stop("'gseaObj' must be a gsea object from clusterProfiler package")
	}
	dot_df=as.data.frame(gseaObj)
	if(nrow(dot_df)==0){
		print("empty dataset...")
		return(0)
	}
	dot_df$type = "upregulated"
	dot_df$type[dot_df$NES < 0] = "downregulated"
	#here, put code to threshold the df according to some parameters
	#select TOP NES results 
	maxvals=abs(dot_df$NES)
	maxvals=maxvals[order(-maxvals)]
	thresh=maxvals[top]
	#if thresh is NA it means that we have < thresh elements and cannot subset based on top
	if(!is.na(thresh)){
		dot_df=dot_df[abs(dot_df$NES)>=thresh,]
	}
	#select padj threshold
	dot_df=dot_df[dot_df$p.adjust<=padj_thresh,]
	#select NES threshold
	dot_df=dot_df[abs(dot_df$NES)>=NES_thresh,]

	if(nrow(dot_df)==0){
		print("maybe threshold parameters too stringent...")
		return(0)
	}

	#calculate the TRUE GeneRatio; find length of each geneset and the numeber of intersecting gnes
	pos=match(as.character(dot_df$ID),names(gseaObj@geneSets))
	lengths=sapply(gseaObj@geneSets,length)
	dot_df$length_geneSet=lengths[pos]
	#find GeneRatio (setSize is the number of intersecting genes of query and geneSet):
	dot_df$GeneRatio=round(dot_df$setSize/dot_df$length_geneSet,2)
	dot_df$neglog10padj=-log10(dot_df$p.adjust)
	minpadj=min(dot_df$neglog10padj)
	maxpadj=max(dot_df$neglog10padj)
	if(!is.na(thresh)){
		limline=max(thresh,NES_thresh)
	}else{
		limline=NES_thresh
	}
	p <- ggplot(dot_df, aes(x = NES, y = fct_reorder(Description, NES))) + 
	               geom_point(aes(size = GeneRatio, color = neglog10padj)) +
	               theme_bw(base_size = 14) +
	        scale_colour_gradient(limits=c(minpadj, maxpadj), low="red") +
	        ylab(NULL) #+
	        #ggtitle("title of the plot")+
	        #geom_vline(xintercept=c(-limline,limline), linetype="dashed",size=1)
	return(p) 
}









########################################################################################
########################################################################################
#bulkRNAseq, compare multiple comparisons and heatmaps
########################################################################################
########################################################################################




#function for high-level bulk RNA-Seq analysis from DESeq2/edgeR results
matrixFromRNAresults<-function(resultsList,padj_thresh=0.01,log2FCthresh=0.5,orderingDriver=1:length(resultsList),orderingtype="kmeans",centers=4,...) {
	#resultsList: list of dataframes of results of differential gene analyses. For example, results from
	#				edgeR or DESeq2. Must have at least a column named "log2FoldChange" and "padj", for
	#				log 2 fold change in gene expression and statistical significance, respectively. Names
	#				of that list will be used as labels in the heatmap and in the same order provided.
	#				rownames of these data frames must be the genes
	#padj_thresh, log2FCthresh: threshold to define DEGs. For the heatmap, the function consider all the genes
	#				that are considered DEGS (with that parameters) in at least one condition
	#orderingDriver: who drives the cluster/ranking? Default: log2 FC of all the comparisons in the list provided. If
	#				only the second comparison should drive the clustering, put 2. If "ranking" in orderingtype: put only one number
	#orderingtype: kind of clustering or ranking. Can be either "kmeans" or "hierarchical" # "ranking" 
	#centers: the number of clusters. Useless if "ranking"
	#...: other parameters for the clustering (for example, centers,nstart and iter.max for kmeans)

	#check input dataframes:
	if(class(resultsList)!="list"){
		stop("'resultsList' must be a list")
	}

	if(length(resultsList)<1){
		stop("'resultsList' must be a list with at least one element")
	}
	for (i in 1:length(resultsList)){
		if(class(resultsList[[i]])!="data.frame" & class(resultsList[[i]])!="matrix"){
			stop(paste("The element number",i,"of 'resultsList' is not a data.frame nor a matrix..."))
		}else{
			#check required columns are present and are numeric
			if(!"padj"%in% colnames(resultsList[[i]]) | ! "log2FoldChange"%in%colnames(resultsList[[i]])){
				stop(paste("The element number",i,"of 'resultsList' does not have the 'padj' and 'log2FoldChange' columns..."))
			}else{
				#check if numeric
				pos_padj=match("padj",colnames(resultsList[[i]]))
				pos_log2FC=match("log2FoldChange",colnames(resultsList[[i]]))
				if ( class(resultsList[[i]][,pos_padj])!="numeric" |  class(resultsList[[i]][,pos_log2FC])!="numeric"){
					stop(paste("The element number",i,"of 'resultsList' have padj or log2FoldChange columns that are not numeric..."))
				}else{
					notNA=resultsList[[i]][,pos_padj] [!is.na(resultsList[[i]][,pos_padj])]
					if(min(notNA)<0 | max(notNA)>1){
						stop(paste("The element number",i,"of 'resultsList' have a padj > 1 or <0..."))
					}
				}
			}
		}
	}

	if(length(orderingDriver)>length(resultsList)){
		stop("'orderingDriver' indexes must not be longer than 'resultsList' length")
	}
	if(class(orderingDriver)!= "numeric" & class(orderingDriver)!= "integer"){
		stop("'orderingDriver' must be number(s)")
	}
	if(any(orderingDriver>length(resultsList) | orderingDriver<1)){
		stop("Some of the 'orderingDriver' are greater than the length of 'resultsList' or <1...")
	}

	if(orderingtype!="kmeans"&orderingtype!="hierarchical" & orderingtype!="ranking"){
		stop("'orderingtype' must be either 'kmeans' or 'hierarchical' or 'ranking'")
	}



	#take DEGs for each result
	degs=as.list(rep(NA,length(resultsList)))
	for(i in 1:length(resultsList)){
		current=resultsList[[i]]
		pos_padj=match("padj",colnames(current))
		pos_log2FC=match("log2FoldChange",colnames(current))		
		degs[[i]]=rownames(current[!is.na(current[,pos_padj]<padj_thresh) & current[,pos_padj ] <padj_thresh &abs(current[,pos_log2FC])>log2FCthresh ,])
	}
	allDEGs=Reduce("union",degs)

	#create the matrix:
	df_log2FC=df_padj=matrix(rep(NA,length(resultsList)*length(allDEGs)),ncol=length(resultsList))
	for(i in 1:length(resultsList)){
		current=resultsList[[i]]
		pos=match(allDEGs,rownames(current))
		df_log2FC[,i]=current[pos,]$log2FoldChange
		df_padj[,i]=current[pos,]$padj
	}
	colnames(df_log2FC)=colnames(df_padj)=names(resultsList)
	rownames(df_log2FC)=rownames(df_padj)=allDEGs

	#remove NAs: genes with log2FC==NA or gene symbols/ID present only in one of the elements provided
	toremove=apply(df_log2FC,1,function(i){any(is.na(i))})
	df_log2FC=df_log2FC[!toremove,]
	df_padj=df_padj[!toremove,]
	print(paste("Removed",table(toremove)["TRUE"],"genes with at least one NA"))

	#randomness. Randomly shuffle the table before clustering and randomly pick colors
	set.seed(123)
	shuffle=sample(nrow(df_log2FC),replace=FALSE)
	df_log2FC=df_log2FC[shuffle,]
	df_padj=df_padj[shuffle,]

	# cluster matrix, reorder and set color scales
	if(orderingtype=="kmeans"){
		clustobject <- kmeans(df_log2FC[,orderingDriver],centers=centers, ...)
		ord=order(clustobject$cluster)
		clustnum_ordered=clustobject$cluster[ord]
	}else if (orderingtype=="hierarchical"){
		clustobject <- hclust(dist(df_log2FC[,orderingDriver],method="euclidean"), ...)
		clst=cutree(clustobject,k=centers)
		ord=clustobject$order
		clustnum_ordered=clst[ord]
	}else{
		#orderingtype=="ranking"
		#if orderingDriver is a single number, simple ranking for that column.
		#if is multiple number, rank by the sum of the rows for the columns selected
		sums=apply(df_log2FC[,orderingDriver,drop=FALSE],1,sum)
		ord=order(-sums)
		clustnum_ordered=rep(1,nrow(df_log2FC))
	}
	
	df_log2FC_ord=df_log2FC[ord,]
	df_padj_ord=df_padj[ord,]

	#transform pvalue of NA in 1
	df_padj_ord[is.na(df_padj_ord)]=1

	#return matrix with values, ordered by clustering, and the clustering
	res=list(df_log2FC_ord,df_padj_ord,clustnum_ordered)
	names(res)=c("mat","padj","clusters")
	return(res)




}

#extract the part of the matrix from matrixFromRNAresults function corresponding to the desired cluster
clustfromRNA<-function(matrixResults,clusters) {
	#matrixResults: result from matrixFromRNAresults function
	#cluster: a number indicating which cluster to extract
	if(class(matrixResults)!="list" | length(matrixResults)!=3){
		stop("maybe 'matrixResults' is not the result of matrixFromRNAresults")
	}
	nclust=	length(table(matrixResults$clusters))
	if(any(clusters>nclust | clusters<1)){
		stop("'clusters' must be in the clusters of matrixFromRNAresults")
	}
	df_log2FC_ord=matrixResults$mat
	df_padj=matrixResults$padj
	clustnum_ordered=matrixResults$clusters	
	pos=clustnum_ordered%in%clusters

	df_log2FC_ord=df_log2FC_ord[pos,]
	df_padj=df_padj[pos,]
	clustnum_ordered=clustnum_ordered[pos]
	res=list(df_log2FC_ord,df_padj,clustnum_ordered)
	names(res)=c("mat","padj","clusters")
	return(res)
}


#function to plot the RNA-Seq bulk multiple comparisons (input is the result of matrixFromRNAresults function)
plotRNAresults<-function(matrixResults,quantile_saturation=0.98,display_significance="lateral",color_clusters=NULL,
						palette=c("darkgreen","grey80","darkred"),signif_thresh=0.01,breaks=10,fileName=NA,rownames_thresh=50,...) {
	#matrixResults: result from matrixFromRNAresults function. This is a list of 3 elements: mat (matrix
	#				of log2FC), padj (matrix of stat. significance), clusters (ehich gene belongs to which cluster)
	#quantile_saturation: quantile value (of log2FC) at which the color intensity of the heatmap is saturated.
	#display_significance: the way to display the padjusted: "lateral" (lateral bars), "inner" (asterisk,
	#						useful when very few genes are displayed), "none" (not displayed)
	#color_clusters: colors to assign to each cluster in the lateral bar. Default NULL: random colors
	#palette: color palette for the log2 fold change representation in the heatmap. Must be array of valid colors
	#signif_thresh: threshold for defining a comparison as "significant"
	#breaks: how many parts the color scale is divided
	#rownames_thresh: threshold of row numbers of the matrix below which show names in the heatmap (default 50)
	#...: parameters for pheatmap
	library(pheatmap)
	if(class(matrixResults)!="list" | length(matrixResults)!=3){
		stop("maybe 'matrixResults' is not the result of matrixFromRNAresults")
	}
	if(quantile_saturation<0 | quantile_saturation>1){
		stop("'quantile_saturation' must be between 0 and 1")
	}
	if(display_significance!="lateral" & display_significance!="inner"&display_significance!="none"){
		stop("'display_significance' must be 'lateral' or 'inner' or 'none'")
	}
	if(!is.null(color_clusters)){
		if(class(color_clusters)!="character"){
			stop("'color_clusters' must be a character vector of valid colors")
		}
		if(length(color_clusters)<length(table(matrixResults$clusters))){
			stop(paste("'color_clusters' must be either NULL or at least of the same length of the number of clusters (",length(table(matrixResults$clusters)),")"))
		}		
	}

	if(signif_thresh>1 | signif_thresh<0){
		stop("'signif_thresh' must be between 0 and 1")
	}

	#extract info from input
	df_log2FC_ord=matrixResults$mat
	df_padj=matrixResults$padj
	clustnum_ordered=matrixResults$clusters

	nclusters=length(table(clustnum_ordered))

	#digitalize padj, to define "significant" and "not significant"
	df_padj[df_padj<=signif_thresh]=0
	df_padj[df_padj>signif_thresh]=1
	colnames(df_padj)=paste("padj",colnames(df_padj))
	df_padj=df_padj[,ncol(df_padj):1]

	#define colors of clustering
	if(is.null(color_clusters)){
		set.seed(123)
		color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
		color=sample(color,nclusters)
	}else{
		color=color_clusters[1:nclusters]
	}

	#determine how padj will be displayed
	if(display_significance=="lateral"){
		x=apply(df_padj,2,as.character)
		rownames(x)=rownames(df_padj)
		df_ann=cbind(as.data.frame(x),data.frame(cluster=as.character(clustnum_ordered)))
		ann_colors = list()
		for(i in 1:ncol(df_padj)){
			ann_colors[[i]]=c("blue4","white")
			names(ann_colors[[i]])=c("0","1")
			names(ann_colors)[i]=colnames(df_padj)[i]
		}
		ann_colors[[length(ann_colors)+1]]=color
		names(ann_colors[[length(ann_colors)]])=as.character(unique(clustnum_ordered))
		names(ann_colors)[length(ann_colors)]="cluster"
	}else {
		df_padj[df_padj==0]="*"
		df_padj[df_padj==1]=""	
		df_ann=data.frame(cluster=as.character(clustnum_ordered))
		rownames(df_ann)=rownames(df_log2FC_ord)
		ann_colors = list(
			cluster=color
		)
		names(ann_colors[[1]])=	as.character(unique(clustnum_ordered))
	}
	
	ambsnummax=quantile(abs(df_log2FC_ord),quantile_saturation)
	brk=c( seq( -ambsnummax , ambsnummax,ambsnummax/breaks))
	my_palette <- colorRampPalette(palette)(n = length(brk)-1 )

	#threshold for showing the names of the genes
	if(nrow(df_log2FC_ord)>rownames_thresh){
		show=FALSE
	}else{
		show=TRUE
	}

	if(display_significance=="lateral"){
		pheatmap(df_log2FC_ord,cluster_rows=FALSE,cluster_cols=FALSE,scale="none",breaks=brk,
					color=my_palette,legend=TRUE,show_rownames = show,annotation_row=df_ann,
					annotation_colors=ann_colors,filename=fileName)		
	}else if(display_significance=="inner"){
		pheatmap(df_log2FC_ord,cluster_rows=FALSE,cluster_cols=FALSE,scale="none",breaks=brk,
					color=my_palette,legend=TRUE,show_rownames = show,annotation_row=df_ann,
					annotation_colors=ann_colors,display_numbers=df_padj_ord,fontsize_number=1,filename=fileName)			
	}else{
		pheatmap(df_log2FC_ord,cluster_rows=FALSE,cluster_cols=FALSE,scale="none",breaks=brk,
					color=my_palette,legend=TRUE,show_rownames = show,annotation_row=df_ann,
					annotation_colors=ann_colors,filename=fileName)		
	}

}





###############################
# gene ontology from ChroKit framework (https://github.com/ocroci/ChroKit)
###############################


#the code for this function was inspired from enricher_internal function in DOSE package 
#(author Guangchuang Yu, http://guangchuangyu.github.io)
GOcalc<-function(gene,terms,minsize=10,maxsize=500,padj_method="BH") {
  #gene: list of genes
  #term: a data frame with term->gene associations, resulting from readGMT function
  #       or a combination of different catgories
  #minsize and maxsize: filter of size of genesets to consider
  #padj_method: method to calculate the padjusted from hypergeometric test
  if(class(gene)!="character"){
    stop("'gene' must be a character vector with gene symbols")
  }

  if(minsize<0 | maxsize<0){
    stop("'minsize' and 'maxsize' must be > 0")
  }
  if(!(padj_method%in%c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none"))){
    stop("wrong padj method given")
  }

  gene=unique(gene)
  gene=gene[!is.na(gene)]

  if(length(gene)==0){
    stop("gens in query is empty...")
  }

  #take only those queried:
  pos=which(as.character(terms$gene)%in%gene)

  if(sum(pos)==0){
    #no genes provided matches those in genesets
    return(NULL)
  }

  
  queried_terms=terms[pos,]
  keyval=split(as.character(queried_terms[,2]),as.character(queried_terms[,1]))
  #retrieve universe (all genes in the term GMT provided)
  #if universe is custom, put here the universe as INTERSECTION.
  #in case, intersect also with keyval
  allgenes=unique(as.character(terms$gene))
  # if(!is.null(background)){
  #   allgenes=intersect(allgenes,background)
  # }
  
  #take all genes of terms specified found with at last one of our genes:
  allkeyval=split(as.character(terms[,2]),as.character(terms[,1]))
  pos=names(allkeyval)%in%names(keyval)
  allkeyval=allkeyval[pos]
  
  #filter keyval and allkeyval for min and max
  sizes=sapply(allkeyval,length)
  idx= sizes>=minsize & sizes<=maxsize
  #if no idx, no geneset satisfies criteria
  if(sum(idx)==0){
    return(NULL)
  }
  
  keyval=keyval[idx]
  allkeyval=allkeyval[idx]

  IDnames=names(keyval)
  
  #do hypergeometric test
  query=sapply(keyval, length)
  totalBG=sapply(allkeyval, length)
  genes_selected=sum(gene%in%as.character(terms$gene))
  
  #calculate pvalues with hypergeometric model
  HyperGeomTests=sapply(1:length(query),function(i){
                                  phyper(q=query[i]-1, m=totalBG[i], n=length(allgenes)-totalBG[i], k=genes_selected, lower.tail = FALSE, log.p = FALSE)
                                }
                        )
  
  tocalc=lapply(1:length(query),  function(i){
                lst=list(unname(query[i][1]),unname(genes_selected),unname(totalBG[i]),unname(length(allgenes)))
                names(lst)=c("overlap","query","geneset","universe")
                return(lst)
  })
  tocalc=do.call("rbind",tocalc)
  

  ## gene ratio: selected genes / genes for each geneset
  ## background ratio: ratio between all genes in a geneset (regardless the query) and the universe (union of all genes of all genesets used)

  gene_ratio = unlist(tocalc[,1])/unlist(tocalc[,2])
  background_ratio= unlist(tocalc[,3])/unlist(tocalc[,4])
  
  HyperGeomTests_padj=p.adjust(HyperGeomTests, method=padj_method)

  geneID = sapply(keyval, function(i) paste(i, collapse="/"))
  
  finaldf=data.frame(Term=IDnames,
                    Genes=geneID,
                    overlap=unlist(tocalc[,1]),
                    query=unlist(tocalc[,2]),
                    geneset=unlist(tocalc[,3]),
                    universe=unlist(tocalc[,4]),
                    gene_ratio=gene_ratio,
                    background_ratio=background_ratio,
                    pval=HyperGeomTests,
                    padj=HyperGeomTests_padj)
  #rank according to padj:
  ord=order(finaldf$padj)
  finaldf=finaldf[ord,]
  
  if(nrow(finaldf)==0){
    return(NULL)
  }
  return(finaldf)

}



#function for filter matrix of results from GO (padj, generatio....)
filterGOres<-function(GOres,padjthresh,generatiothresh,topN) {
  #inputs: the result matrix from GO analysis (serverROI) with thresholds from input GUI
  #filter padj (at least one ROI must be <thresh)

  mat=GOres
  pos_tokeep=rep(TRUE,nrow(mat))

  partdfpadj=mat[,grepl("_padj$",colnames(mat)),drop=FALSE]
  
  if(min(dim(mat))>0){
    tokeeppadj=apply(partdfpadj,1,function(i){any(i<padjthresh)})
    pos_tokeep=pos_tokeep&tokeeppadj
  }

  #filter gene ratio
  if(min(dim(mat))>0){
    partdfgeneratio=mat[,grepl("_gene_ratio$",colnames(mat)),drop=FALSE]
    tokeepgeneratio=apply(partdfgeneratio,1,function(i){any(i>=generatiothresh)})
    pos_tokeep=pos_tokeep&tokeepgeneratio     
  }
  


  #filter topN elements
  if(min(dim(mat))>0){
    #filter with previous filters:
    partdfpadjfilt=partdfpadj[pos_tokeep,]
    minvals=apply(partdfpadj,1,min)
    minvalsfilt=minvals[pos_tokeep]
    #keep the first N positions of the already filtered matrix (top 10 for example of the matrix
    #already filtered for the other parameters)
    minvalsfilt=sort(minvalsfilt)
    if(length(minvalsfilt)>topN){
      #in this case we have more elements than topN => filter
      val=minvalsfilt[topN]
      tokeeptopN=minvals<=val
      pos_tokeep=pos_tokeep&tokeeptopN

    }else{

    }
  }

  return(pos_tokeep)
}














###############################################################################
# functions for coverage (taken from ChroKit: https://github.com/ocroci/ChroKit)
################################################################################






#this function calculate the pileup for each base pair for each range of the ROI
#is a little improvement of the GRbaseCoverage function from compEpiTools (lapply intead of for)
#a better strategy to save RAM could be to keep the norm. factor associated to each BAM, and use it when needed.
#in that case it is possible to store everything as integer, more than 2 billions values 
#https://stackoverflow.com/questions/23660094/whats-the-difference-between-integer-class-and-numeric-class-in-r
GRbaseCoverage2<-function(Object, signalfile,signalfileNorm=NULL,signalControl=NULL,signalControlSpike=NULL,multiplFactor=1e+06)
{
    #Object: genomic range in which calculate the base coverage
    # signalfile= path to the signal file file (the associated signal file index must be present in the same directory if bam file)
    #      or the bw/bigwig file 
    #signalfileNorm= path to the signal file (BAM) for which normalize the signalfile (can be the same file for library normalization,
    #             or the spike in bam, for spike in normalization)
    #signalControl=path to the signal file (BAM) of the input for spike-in normalization (ChIP-Seq). If this option
    #         is provided, signalControlSpike must be also provided
    #signalControlSpike=path to the signal file (BAM) of the spike-in in the input (ChIP-Seq). If this option is
    #         provided, also signalControl must be provided
    #multiplFactor= constant factor to multiply the final normalization coefficient, usually 1 million
    if (!is.character(signalfile)) {
      stop("signalfile has to be a file path of class character...")
    }
    if(!file.exists(signalfile)){
      stop("'signalfile' doesn't exist...")
    }

    if(!is.null(signalfileNorm)){
      if(!is.character(signalfileNorm)){
        stop("'signalfileNorm' has to be a file path of class character...")
      }
      if(!file.exists(signalfileNorm)){
        stop("'signalfileNorm' doesn't exist...")
      }
    }

    #check the existence of control and spikein
    if(!is.null(signalControl)){
      if(!is.character(signalControl)){
        stop("'signalControl' has to be a file path of class character...")
      }
      if(!file.exists(signalControl)){
        stop("'signalControl' doesn't exist...")
      }
    }

    if(!is.null(signalControlSpike)){
      if(!is.character(signalControlSpike)){
        stop("'signalControlSpike' has to be a file path of class character...")
      }
      if(!file.exists(signalControlSpike)){
        stop("'signalControlSpike' doesn't exist...")
      }
    }   
    #if provide signalControl you must provide also signalControlSpike
    if (xor(is.null(signalControl),is.null(signalControlSpike))){
      stop("If 'signalControl' is provided, also 'signalControlSpike' must be provided, and vice versa...")
    }

    #if you provide signalControlSpike and signalControl you must also have provided signalfileNorm
    if(!is.null(signalControl) & is.null(signalfileNorm)){
      stop("If 'signalControl' and 'signalControlSpike' are provided, also 'signalfileNorm' must be provided")
    }


    #select if bam:
    if(substring(signalfile,nchar(signalfile)-3,nchar(signalfile))==".bam"){
      BAMseqs <- names(scanBamHeader(signalfile)[[1]]$targets)
      matchingSeqs <- which(as.character(seqnames(Object)) %in% BAMseqs)
      if (length(matchingSeqs) == 0) 
          return(sapply(width(Object), function(x) rep(0, x)))
      param <- ApplyPileupsParam(which = Object[matchingSeqs],what = "seq")
      pileupFiles=PileupFiles(signalfile)
      coverage <- applyPileups(pileupFiles, FUN = function(x) x,param = param)
      rm(pileupFiles)
      widths <- width(Object[matchingSeqs])
      covList <- list()
      starts <- start(Object[matchingSeqs])
      covList=lapply(1:length(Object[matchingSeqs]), function(i) {
          covx <- coverage[[i]]
          cvec <- rep(0, widths[i])
          inds <- covx$pos - starts[i] + 1
          #sometims max(inds) is > than length of cvec. => NAs in some positions. Why?
          cvec[inds] <- colSums(covx$seq)
          return(cvec)
      })
      covList<-lapply(covList,function(i)as.integer(i))
      rm(coverage)


      ####################################################################
      # normalization
      ####################################################################

      if(!is.null(signalfileNorm)){
          param <- ScanBamParam(flag = scanBamFlag(isUnmappedQuery = FALSE))
          nreads <- countBam(signalfileNorm, param = param)$records
          norm_factor=multiplFactor/nreads

          if(!is.null(signalControl)) {
            nreadCtr=countBam(signalControl, param = param)$records
            nreadCtrSpike=countBam(signalControlSpike, param = param)$records
            norm_factor=norm_factor* nreadCtrSpike/  nreadCtr         
          }
          print (paste("Norm. factor=",norm_factor))  
          #covList <- lapply(covList, function(x) (x*norm_factor) )
      }else{
      	norm_factor=1
      }

      

      if (length(matchingSeqs) < length(Object)) {
          coverageTot <- sapply(width(Object), function(x) list(rep(0, 
              x)))
          coverageTot[matchingSeqs] <- covList
      }
      else coverageTot <- covList

      rm(covList)
      return(list(coverageTot,norm_factor))


    #select if wig:
    }else if (substring(signalfile,nchar(signalfile)-2,nchar(signalfile))==".bw" | tolower(substring(signalfile,nchar(signalfile)-6,nchar(signalfile)))==".bigwig"){
    
      #open wig file only for names and lengths of chromosomes, to find
      #those matching with "Object"
      print ("    Coverage of WIG file...")
      chromosomes=seqlevelsInUse(Object)
      #just open the wig file, without opening all chromosome content. This serves just to know 
      #how long is each chromosome
      grfake=GRanges(Rle(chromosomes[1]),IRanges(1,1))
      wig=import(signalfile,which=grfake,as = 'Rle')
      chrswig=lengths(wig)
      rm(wig)
      common_chromosomes=intersect(chromosomes,names(chrswig))
      #initialize 0s for all ranges (if wig do not have a chr, the remaining are 0s)
      widths=width(Object)
      coverageTot=lapply(1:length(Object),function(i){integer(widths[i])})
      for (i in 1:length(common_chromosomes)){
        #length(as.character(seqnames(Object))) should be equal to length(Object) =>
        #find positions of Object that have that chromosome
        #opening chromosome by chromosome
        pos=as.character(seqnames(Object))==common_chromosomes[i]
        chrsize=chrswig[common_chromosomes[i]]
        grcurrentchr=GRanges(Rle(common_chromosomes[i]),IRanges(1,chrsize))
        wigtempchr=import(signalfile,which=grcurrentchr,as = 'Rle')
        #count coverage for that chromosome in wig in the Object ranges in the same chromosome
        counts=Views(unlist(wigtempchr[[common_chromosomes[i]]]),ranges(Object[pos]))
        rm(wigtempchr)
        #extract the base coverage for all those ranges for ith chromosome.
        #used "as.integer" to save space: all numbers are pure integers
        tmp=viewApply(counts, as.integer)
        
        if(length(counts)==1){
          tmp=list(as.vector(tmp))
        }
        if(class(tmp)=="matrix"){
          ##check whether the order is correct
          tmp=lapply(seq_len(ncol(tmp)), function(i) tmp[,i])
        }
        rm(counts)
        ##IMPORTANT##
        #############
        #change here if numbers can not be integer, but real
        #############
        tmp<-lapply(tmp,function(i)as.integer(i))
        ########################################################
        coverageTot[pos] <- tmp
      }
      #return coverage and norm factor (in case of WIG, simply 1) to keep the same result structure 
      #for the function
      return(list(coverageTot,1))      
    }else{
      stop("Error in retrieving file/extension...")
    }

}




#functions that takes GRbaseCoverage and output GRcoverageInbins output
#taken from GRcoverageInbins function from compEpiTools
# and re-implemented in Rcpp
makeMatrixFrombaseCoverageCPP <- cxxfunction(signature(GRbaseCoverageOutput='List',Nbins='integer',Snorm="integer"), plugin='Rcpp', body = '  
     Rcpp::List xlist(GRbaseCoverageOutput); 
     int nbins = Rcpp::as<int>(Nbins);
     int snormbool=Rcpp::as<int>(Snorm);
     int n = xlist.size(); 
     //std::vector<double> res(n);  
     double res;
     //define final matrix output: ncol=nbin, nrow=n
     //double mat[n][nbins];
     Rcpp::NumericMatrix mat( n , nbins );

     std::vector<int> lengths(n);

     for(int i=0; i<n; i++) {     
         SEXP ll = xlist[i]; 
         Rcpp::NumericVector y(ll);  
         int m=y.size();   
         //divide m (size) by the number of bins: how many elements to be summed for each bin?
         int goodpart=m/nbins;
         //printf("size of a bin: %d\\n",goodpart);

         // for each bin, sum elements in each bin
         for(int k=0; k<nbins; k++){
          res=0;
          for(int j=k*goodpart; j<(k+1)*goodpart; j++){     
              res=res+y[j]; 
          }   
          mat(i,k)=res;
         }

         //populate the array of lengths if Snorm=TRUE
         lengths[i]=m;
     }
     //if Snorm=TRUE (Snorm>0), divide each matrix value for m
     if(snormbool>0){
      for(int i=0; i<n; i++){
        for(int k=0; k<nbins; k++){
          mat(i,k) = mat(i,k)/lengths[i];
        }
      }
     }
       
   return mat;  
') 

#wrapper for makeMatrixFrombaseCoverageCPP function
makeMatrixFrombaseCoverage <-function(GRbaseCoverageOutput,Nbins,Snorm=FALSE,norm_factor=1) {
    # GRbaseCoverageOutput: output from GRbaseCoverage2 function
    # Nbins: the number of bins to ivide the coverage for each range into
    # Snorm: whether to normalize the coverage for each bin for the length of the range (TRUE/FALSE)
    if(Snorm==TRUE){
      norm=1
    }else{
      norm=0
    }
    #here introduce normalization step. We could have done easyly in the CPP function,
    #but it's fast anyway
    result=makeMatrixFrombaseCoverageCPP(GRbaseCoverageOutput,Nbins,Snorm=norm)
    if(norm_factor!=1){
    	result <- result*norm_factor
    }
    return(result)
}
