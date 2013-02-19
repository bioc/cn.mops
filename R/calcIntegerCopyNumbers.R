#' @title Calculation of integer copy numbers for the CNVs and CNV regions.
#' @description This generic function calculates the integer copy numbers of
#'  a CNV detection method stored in an instance of 
#' \code{\link{CNVDetectionResult-class}}. 
#' 
#' @param object An instance of "CNVDetectionResult" 
#' @examples
#' data(cn.mops)
#' r <- cn.mops(X[1:100,1:5])
#' calcIntegerCopyNumbers(r)
#' @return \code{calcIntegerCopyNumbers} returns an 
#' instance of "CNVDetectionResult".
#' @author Guenter Klambauer \email{klambauer@@bioinf.jku.at}
#' @export
#' @importFrom IRanges findOverlaps
#' @importFrom IRanges as.list
#' @importFrom IRanges as.matrix
#' @importFrom GenomicRanges values


setMethod("calcIntegerCopyNumbers", signature="CNVDetectionResult",
		definition = function(object){
			
			
			priorImpact <- object@params$priorImpact
			cyc <- object@params$cyc
			classes <- object@params$classes
			I <- object@params$folds
			minReadCount <- object@params$minReadCount
			X <- object@normalizedData
			cnvr <- object@cnvr
			segmentation <- segmentation(object)
			cnvs <- cnvs(object)
			gr <- object@gr
			cov <- rep(1,ncol(X))
			uT <- object@params$upperThreshold
			lT <- object@params$lowerThreshold
			mainClass <- object@params$mainClass
			
			if (length(cnvr)==0 | length(cnvs)==0)
				stop(paste("No CNV regions in result object. Rerun cn.mops",
								"with different parameters!"))
			
			
			# for CNV regions
			M <- IRanges::as.list(IRanges::findOverlaps(cnvr,object@gr))
			XX <- lapply(M,function(i){ 
						if (length(i)>=3) ii <- i[-c(1,length(i))]
						apply(X[ii, ,drop=FALSE],2,mean) })
			CN <-t(sapply(XX,function(x) .cn.mopsC(x,I=I,
										classes=classes,
										cov=cov,priorImpact=priorImpact,
										cyc=cyc,
										minReadCount=minReadCount)$expectedCN))
			
			colnames(CN) <- colnames(X)
			resObject <- object
			GenomicRanges::values(cnvr) <- CN
			resObject@cnvr <- cnvr
			
			## now for individual CNVs and segmentation
			
			iCN <- rep(mainClass,length(segmentation))
			#mapping from CNVs to segmentation
			csM <- IRanges::as.matrix(IRanges::findOverlaps(segmentation,cnvs))
			tmpIdx <- which(values(segmentation)$sampleName[csM[,1]]==values(cnvs)$sampleName[csM[,2]])
			csM <- csM[tmpIdx, ,drop=FALSE]
			idx <- csM[,1]
			
			M2 <- IRanges::as.list(IRanges::findOverlaps(segmentation[idx],object@gr))
			XX2 <- lapply(M2,function(i){ 
						if (length(i)>=3) ii <- i[-c(1,length(i))]
						apply(X[ii, ,drop=FALSE],2,mean) })
			CN2 <-t(sapply(XX2,function(x) .cn.mopsC(x,I=I,
										classes=classes,
										cov=cov,priorImpact=priorImpact,
										cyc=cyc,
										minReadCount=minReadCount)$expectedCN))
			colnames(CN) <- colnames(X)
			extractedCN <- CN2[cbind(1:length(idx),
							match(as.character(values(segmentation[idx])$sampleName),
									colnames(X)))]
			iCN[idx] <- extractedCN 
			GenomicRanges::values(segmentation)$CN <- iCN
			GenomicRanges::values(cnvs)$CN <- extractedCN[csM[,2]]
			resObject@cnvs <- cnvs
			resObject@segmentation <- segmentation
			
			return(resObject)							
		})

