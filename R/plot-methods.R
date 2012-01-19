# Copyright (C) 2011 Klambauer Guenter 
# <klambauer@bioinf.jku.at>

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
		function(r, seqname, sampleIdx, ...) {
			standardGeneric("segplot")
		})

#' Plots the log normalized read counts and the detected segments for 
#' one sample on one reference sequence.
#' 
#' @param r An instance of "CNVDetectionResult" 
#' @param seqname The name of the reference sequence or chromosome to be
#' plotted.
#' @param sampleIdx The index of the sample as it appears in the read count
#' matrix.
#' @examples
#' data(cn.mops)
#' r <- cn.mops(X[1:200, ])
#' segplot(r,sampleIdx=1)
#' @return Generates a segmentation plot.
#' @author Guenter Klambauer \email{klambauer@@bioinf.jku.at}
#' @export

setMethod("segplot",
		signature(r="CNVDetectionResult"),
		function(r, seqname, sampleIdx) {
			if (!"DNAcopy" %in% rownames(installed.packages())){
				stop("DNAcopy must be installed for these plot types.")
			}
			
			library(DNAcopy)
			X <- r@normalizedData
			if (all(IRanges::width(IRanges::ranges(X))==1) ){
				# idx mode
				WL <- 1
			} else {
				wir <- width(IRanges::ranges(X))
				WL <- wir[1]
				wir <- wir[-length(wir)]
				if (!all(wir==WL)){
					stop("Plot function only for equally spaced segments.")
				}
			}
			genomdat <-	do.call("cbind",values(X)@listData)
			if (is.null(colnames(X))){colnames(genomdat) <- as.character(
						paste("Sample",1:ncol(
										genomdat),sep="_"))} else {
				colnames(genomdat) <- paste("Sample",colnames(X),sep="_")
			}
			if (missing(sampleIdx)){
				sampleIdx <- 1
				message(paste("Missing \"sampleIdx\" argument. Selecting",
								sampleIdx,".\n"))
			}
			
			if (!is.null(r@params$lambda)){
				genomdat <- (genomdat/r@params$lambda[,"CN2"])
				#cat("lambda.\n")
			} else {
				genomdat <- (genomdat/rowMedians(genomdat))	
			}
			genomdat[which(is.na(genomdat))] <- 1
			genomdat[which(genomdat==Inf)] <- max(
					genomdat[which(is.finite(genomdat))])
			genomdat[which(genomdat==-Inf)] <- min(
					genomdat[which(is.finite(genomdat))])
			genomdat <- log2(genomdat)
			genomdat <- pmax(genomdat,min(
							genomdat[which(is.finite(genomdat))]))
			genomdat <- genomdat[,sampleIdx,drop=FALSE]
			
			#colnames(genomdat) <- sampleNames
			if (missing(seqname)){
				seqname <- rep(X@seqnames@values,X@seqnames@lengths)[1]
				message(paste("Missing \"seqname\" argument. Selecting"
								,seqname,".\n"))
			}
			#browser()
			chrIdx <- which(rep(X@seqnames@values,X@seqnames@lengths)==seqname)
			maploc <- as.integer(
					(start(ranges(X[chrIdx]))+end(ranges(X[chrIdx])))/2)
			
			data.type <- "logratio"
			zzz <- data.frame(chrom=as.character(seqname), maploc=maploc, 
					genomdat[chrIdx, ,drop=FALSE])
			attr(zzz, "data.type") <- data.type
			class(zzz) <- c("CNA", "data.frame")
			
			segDataTmp <- IRanges::as.data.frame(segmentation(r),as.is=TRUE)
			## segDataTmp$sampleName <- paste("S",
			##         as.character(segDataTmp$sampleName),sep="_")
			segDataTmp <- segDataTmp[which(segDataTmp$sampleName==
									colnames(genomdat)), ]
			segDataTmp$sampleName <- as.character(segDataTmp$sampleName)
			colnames(segDataTmp) <- c("chrom", "loc.start", "loc.end",
					"num.mark","strand","ID","seg.median","seg.mean","CN")
			segDataTmp$num.mark <- segDataTmp$num.mark/WL
			segDataTmp$chrom <- as.character(segDataTmp$chrom)
			#segDataTmp$ID <- as.character(colnames(genomdat)[1])
			
			
			segres <- list()
			segres$data <- zzz
			segres$output <- segDataTmp[,c("ID","chrom","loc.start","loc.end",
							"num.mark","seg.mean")]
			#segres$segRows <- segDataTmp[, 9:10]
			segres$segRows <- segDataTmp[, 2:3]
			segres$call <- "unknown"    
			class(segres) <- "DNAcopy"
			plot(segres,xmaploc=TRUE)
			#plot(segres,xmaploc=FALSE)
		})

