# Copyright (C) 2011 Klambauer Guenter 
# <klambauer@bioinf.jku.at>

.statmod <- function(x,na.rm=FALSE) {
	if (na.rm){
		z <- table(as.vector(x[!is.na(x)]))
		r <- names(z)[z == max(z)]
		return(as.numeric(r))
	} else {
		if (any(is.na(x))){return(NA)
		} else {
			z <- table(as.vector(x[!is.na(x)]))
			r <- names(z)[z == max(z)]
			return(as.numeric(r))
		}
	}
} 


#' Normalize quantitative NGS data in order to make counts comparable over
#' samples. Scales each samples' reads such that the coverage is even for
#' all samples after normalization. 
#' @param X Matrix of positive real values, where
#' columns are interpreted as samples and rows as genomic regions. An entry is
#' the read count of a sample in the genomic region. 
#' @param chr Character vector that has as many elements as "X" has rows. The
#' vector assigns each genomic segment to a reference sequence (chromosome).
#' @param normType Type of the normalization technique. Each samples'
#' read counts
#' are scaled such that the total number of reads is equal after normalization.
#' By this parameter one can decide to which coverage (i.e. total reads) the 
#' read counts should be normalized. Possible choices are the minimal coverage 
#' ("min"), the mean or median coverage ("mean", "median") or any quantile 
#' ("quant"). If this parameter is set to the value "mode", 
#' the read counts are scaled such that each samples'
#' most frequent value (the "mode") is equal after normalization. If the
#' parameter is set to "poisson" the values are scaled such that the 
#' distribution is (rowwise) close to a Poisson distribution.
#' Possible values are "mean","min","median","quant","poisson, and "mode". 
#' Default = "poisson".
#' @param qu Real value between 0 and 1. Default = 0.25.
#' @examples 
#' data(cn.mops)
#' X.norm <- normalizeChromosomes(X)
#' @return A data matrix of normalized read counts with the same dimensions 
#' as the input matrix X.
#' @author Guenter Klambauer \email{klambauer@@bioinf.jku.at}
#' @export


normalizeChromosomes <-
		function(X,chr,normType="poisson",qu=0.25){
	if (!(normType %in% c("mean","min","median","quant","mode","poisson"))){
		stop(paste("Set TO of normalization to \"mean\"",
						"\"min\", \"median\", \"quant\" or \"mode\"."))
	}
	if (is.vector(X)){X <- matrix(X,nrow=1)}
	
	if (missing(chr)){
		chr <- rep("undef",nrow(X))
	}
	# Sequencing data matrix
	# vector of chromosome - length equal to rows of X
	if (length(chr)!=nrow(X)){
		stop("Length of \"chr\" must be equal to number of rows of \"X\".")}
	chr <- (as.character(chr))
	Y <- matrix(0,nrow=nrow(X),ncol=ncol(X))
	
	for (l in (unique(chr))){
		chrIdx <- which(chr==l)
		Ytmp <- X[chrIdx, ,drop=FALSE]
		idxSG <- apply(Ytmp,1,function(x) all(x<1))
		Ytmp[idxSG, ] <- NA
		
		if (nrow(Ytmp) > 1){
			
			mappedreads <- colSums(Ytmp,na.rm=TRUE)
			if (any(mappedreads==0)){
				warning(paste("There exists a reference sequence with zero reads"
								,"for some samples."))
				mappedreads <- pmax(mappedreads, 1)
			}
			
			if (normType == "min"){
				correctwiththis <-  min(mappedreads,na.rm=TRUE)/mappedreads
				
			} else if (normType=="mean"){
				correctwiththis <-  mean(mappedreads,na.rm=TRUE)/mappedreads
			} else if (normType=="median"){
				cat("normalizing to median \n")
				correctwiththis <-  median(mappedreads,na.rm=TRUE)/mappedreads
			} else if (normType=="quant"){
				correctwiththis <-  quantile(mappedreads,probs=qu,
						na.rm=TRUE)/mappedreads
			} else if (normType=="mode"){
				correctwiththis <-  .statmod(as.vector(Ytmp),na.rm=TRUE)/
						apply(Ytmp,2,.statmod,na.rm=TRUE)
			} else if (normType=="poisson"){
				
				correctwiththis <-  median(mappedreads,na.rm=TRUE)/mappedreads
				YYtmp <- t(t(Ytmp)*correctwiththis)
				v2m <- apply(YYtmp,1,var)/rowMeans(YYtmp)
				mv2m <- median(v2m,na.rm=TRUE)
				if (is.finite(mv2m)) {correctwiththis <- correctwiththis*1/mv2m}
			} else {
				stop("normType not known.")
			}
			if (any(!is.finite(correctwiththis))){
				warning(paste("Normalization for reference sequence ",l,"not", 
								"applicable, because at least one sample has zero",
								"reads."))
				correctwiththis <- rep(1,ncol(X))
			}
			Ytmp <- t(t(Ytmp)*correctwiththis)
			Ytmp[idxSG, ] <- 0
		}
		
		Y[chrIdx, ] <- Ytmp
		
	}
	rownames(Y) <- rownames(X)
	colnames(Y) <- colnames(X)
	return(Y)	
}

