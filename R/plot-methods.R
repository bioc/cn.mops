# Copyright (C) 2011 Klambauer Guenter 
# <klambauer@bioinf.jku.at>

.convertToFastSegRes <- function(mopsres){
	ol <- IRanges::findOverlaps(segmentation(mopsres),
			normalizedData(mopsres))
	ID <- as.character(values(segmentation(mopsres))$sampleName)
	seg.mean <- values(segmentation(mopsres))$mean
	T <- table(queryHits(ol))
	nn <- length(normalizedData(mopsres))
	cs <- cumsum(T)%%nn
	num.mark <- as.integer(T)
	
	startRow <- c(1,(cs+1)[-length(cs)])
	endRow <- cs
	endRow[which(endRow==0)] <- nn
	grRet <- segmentation(mopsres)
	values(grRet) <- data.frame(ID=ID,num.mark=num.mark,seg.mean=seg.mean,
			startRow=startRow,endRow=endRow,stringsAsFactors=FALSE)
	
	return(grRet)
}

.makeLogRatios <- function(mopsres,mainCN="CN2"){
	X <- IRanges::as.matrix(IRanges::values(normalizedData(mopsres)))
	
	if (!is.null(mopsres@params$L[,mainCN])){
		r <- mopsres@params$L[,mainCN]
		#cat("lambda.\n")
	} else {
		r <- rowMedians(X)
	}
	
	R <- X/r
	R[which(is.na(R))] <- 1
	
	R[which(R==Inf)] <- max(R[which(is.finite(R))],na.rm=TRUE)
	
	R[which(R==0)] <- min(R[which(R>0)],na.rm=TRUE)
	R <- log2(R)
	gr <- normalizedData(mopsres)
	colnames(R) <-  unique(as.character(
					IRanges::values(segmentation(mopsres))$sampleName))
	values(gr) <- R
	return(gr)
}



#' Plots read counts, call values and CNV calls in an identified CNV region.
#' 
#' @param x An instance of "CNVDetectionResult" 
#' @param which The index of the CNV region to be plotted.
#' @param margin Vector of two positive integers that states how many segments 
#' left and right of the CNV region should be included in the plot. Default = 
#' c(10,10).
#' @param toFile Logical value whether the output should be plotted to a file.
#' Default = FALSE.
#' @examples
#' data(cn.mops)
#' r <- cn.mops(X[1:200, ])
#' plot(r)
#' @return Generates a CNV calling plot.
#' @author Guenter Klambauer \email{klambauer@@bioinf.jku.at}
#' @noRd
#' @export
#' @importFrom graphics plot

setMethod("plot", signature(x="CNVDetectionResult",y="missing"),
		function(x,which,margin=c(10,10),toFile=FALSE){
			r <- x
			
			if (missing(which)){
				message("Missing argument \"which\". Plotting the first CNVR.")
				which <- 1
			}
			
			MMstart <- IRanges::match(GRanges(seqnames(r@cnvr),
							IRanges(start(r@cnvr),
									start(r@cnvr))), r@normalizedData )
			MMend <- IRanges::match(GRanges(seqnames(r@cnvr),
							IRanges(end(r@cnvr),
									end(r@cnvr))),r@normalizedData )
			
			for (select in which){
				if (!toFile){
					dev.new()
				}
				
				if (select>length(r@cnvr)){
					stop("Selected unknown CNVR for plotting.")}
				refSeq <- as.character(seqnames(r@cnvr)[select])
				start <- MMstart[select]
				end <- MMend[select]
				plotStart <- max(start-margin[1],1)
				plotEnd <- min(end+margin[2],length(r@normalizedData))
				
				if (plotEnd==plotStart){plotStart <- plotStart-1}
				#plot.new()
				layout(matrix(1:3,nrow=1))
				ND <- (do.call("cbind",
									(values(r@normalizedData)@listData)))
				LA <- (do.call("cbind",
									(values(r@localAssessments)@listData)))
				IC <- (do.call("cbind",
									(values(r@individualCall)@listData)))
				
				rn <- paste(start(r@normalizedData),
						end(r@normalizedData),sep=" - ")
				
				rownames(ND) <- rownames(LA) <- rownames(IC) <- rn
				
				refSeqName <- unique(as.character(seqnames(r@cnvr)))
				xlab <- paste(refSeq,": ",
						unlist(start(r@normalizedData))[plotStart],
						" - ",unlist(end(r@normalizedData))[plotEnd],sep="")
				
				lt <- r@params$lowerThreshold
				ut <- r@params$upperThreshold
				col <- rep("black",ncol(ND))
				if (is.numeric(lt) & is.numeric(ut)){
					col[apply(IC[start:end, ,drop=FALSE] >= ut,2,any)] <- "red"
					col[apply(IC[start:end, ,drop=FALSE] <= lt,2,any)] <- "blue"
				}
				
				lty <- sample(2:6,replace=TRUE,size=ncol(ND))
				matplot(ND[plotStart:plotEnd, ],type="l",lty=lty,lwd=2,
						main="Normalized Read Counts",ylab="Read Count",
						xlab=xlab,xaxt="n",col=col)
				axis(1,at=1:length(plotStart:plotEnd),
						labels=FALSE)
				matplot(LA[plotStart:plotEnd, ],type="l",lty=lty,lwd=2,
						main="Local Assessments",ylab="Local Assessment Score",
						xlab=xlab,xaxt="n",col=col)
				axis(1,at=1:length(plotStart:plotEnd),
						labels=FALSE)
				matplot(IC[plotStart:plotEnd, ],type="l",lty=lty,lwd=2,
						main="CNV Call",ylab="CNV Call Value",
						xlab=xlab,xaxt="n",col=col)
				axis(1,at=1:length(plotStart:plotEnd),
						labels=FALSE)
			}
		})


#setMethod("segPlot", signature(r="CNVDetectionResult"),
setGeneric("segplot",
		function(r,mainCN="CN2", sampleIdx, seqnames, plot.type="chrombysample", 
				altcol=TRUE, sbyc.layout, cbys.nchrom=1,
				cbys.layout, include.means=TRUE, zeroline=TRUE,
				pt.pch=".", pt.cex=3, pt.cols=c("green","black"),segcol, 
				zlcol="grey", ylim, lwd=3, ...) {
			standardGeneric("segplot")
		})

#' Plots the log normalized read counts and the detected segments as a 
#' segmentation plot.
#' 
#' @param r An instance of "CNVDetectionResult" 
#' @param mainCN The name of the main copy number. That is "CN2" for diploid
#' individuals. For haplocn.mops this should be set to "CN1".
#' @param sampleIdx The index of the samples to be plotted. (Default = missing)
#' @param seqnames The names of the reference sequence (chromosomes) to
#' be plotted. (Default = missing)
#' @param plot.type the type of plot. (Default = "s").
#' @param altcol logical flag to indicate if chromosomes should be
#'   plotted in alternating colors in the whole genome plot. (Default = TRUE).
#' @param sbyc.layout \code{layout} settings for the multifigure grid layout
#'    for the `samplebychrom' type.  It should be specified as a vector of
#'    two integers which are the number of rows and columns.  The default
#'    values are chosen based on the number of chromosomes to produce a
#'    near square graph.   For normal genome it is 4x6 (24 chromosomes)
#'    plotted by rows. (Default = NULL).
#' @param cbys.layout \code{layout} settings for the multifigure grid layout
#'   for the `chrombysample' type.  As above it should be specified as
#'    number of rows and columns and the default chosen based on the
#'    number of samples. (Default = NULL).
#' @param cbys.nchrom the number of chromosomes per page in the layout.
#'(Default = 1).
#' @param include.means logical flag to indicate whether segment means
#'   are to be drawn. (Default = TRUE).
#' @param zeroline logical flag to indicate whether a horizontal line at
#'    y=0 is to be drawn. (Default = TRUE).
#' @param pt.pch the plotting character used for plotting the log-ratio
#'    values. (Default = ".")
#' @param pt.cex the size of plotting character used for the log-ratio
#'    values (Default = 3).
#' @param pt.cols the color list for the points. The colors alternate
#'    between chromosomes. (Default = c("green","black").)
#' @param segcol the color of the lines indicating the segment means.
#' (Default = "red").
#' @param zlcol the color of the zeroline. (Default = "grey").
#' @param ylim this argument is present to override the default limits
#'    which is the range of symmetrized log-ratios. (Default = NULL).
#' @param lwd line weight of lines for segment mean and zeroline. (Default = 3).
#' @param ... other arguments which will be passed to \code{plot}
#'   commands.
#' @examples
#' data(cn.mops)
#' r <- cn.mops(X[1:200, ])
#' segplot(r,sampleIdx=1)
#' @return Generates a segmentation plot.
#' @author Guenter Klambauer \email{klambauer@@bioinf.jku.at}
#' @export


setMethod("segplot",
		signature(r="CNVDetectionResult"),
		function(r, mainCN="CN2",sampleIdx, seqnames, 
				plot.type="chrombysample", 
				altcol=TRUE, sbyc.layout, cbys.nchrom=1,
				cbys.layout, include.means=TRUE, zeroline=TRUE,
				pt.pch=".", pt.cex=3, pt.cols=c("green","black"),segcol, 
				zlcol="grey", ylim, lwd=3, ...) {
			
			if (!missing(sampleIdx)){
				sn <- sampleNames(r)
				idx <- which(as.character(
								values(segmentation(r))
										$sampleName)==sn[sampleIdx])
				
				r@segmentation <- segmentation(r)[idx]
				nd <- normalizedData(r)
				IRanges::values(nd) <- IRanges::values(nd)[,sampleIdx]
				IRanges::colnames(IRanges::values(nd)) <- sn[sampleIdx]
				r@normalizedData <- nd
				
				
			} 
			
			if (!missing(seqnames)){
				#browser()
				r@segmentation <- segmentation(r)[which(seqnames(
										segmentation(r)) %in% seqnames)]
				nd <- normalizedData(r)
				idx2 <- which(seqnames(nd) %in% seqnames)
				if (length(idx2)==0){
					stop(paste("Given \"seqnames\" do not appear in result",
					"object. Try to exchange \"chr1\" <--> \"1\"."))
				}
				nd <- nd[idx2]
				
				r@normalizedData <- nd
				if (!is.null(r@params$L)){
					r@params$L <- r@params$L[idx2, ]
				}
				
				
			} 

			
			.segPlot(x=cn.mops:::.makeLogRatios(r,mainCN),
					res=cn.mops:::.convertToFastSegRes(r),
					plot.type=plot.type, 
					altcol=altcol, 
					cbys.nchrom=cbys.nchrom,
					include.means=include.means,zeroline=zeroline, 
					pt.pch=pt.pch,pt.cex=pt.cex, pt.cols=pt.cols,
					zlcol=zlcol, ylim=ylim, lwd=lwd, ...)
			
		})
