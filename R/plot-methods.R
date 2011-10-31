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
