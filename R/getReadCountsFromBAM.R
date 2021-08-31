# Copyright (C) 2012 Klambauer Guenter 
# <klambauer@bioinf.jku.at>

#' @title Calculation of read counts from BAM files.
#' @description Generates the read counts from BAM Files. 
#' These counts are necessary for CNV detection methods based
#' on depth of coverage information.
#' 
#' This function can also be run in a parallel version.
#' 
#' @param BAMFiles BAMFiles
#' @param sampleNames The corresponding sample names to the BAM Files. 
#' @param refSeqNames Names of the reference sequence that should be analyzed.
#' The name must appear in the header of the BAM file. If it is not given
#' the function will select the first reference sequence that appears in the
#' header of the BAM files. Can be set to analyze multipe chromosomes at once, 
#' e.g. refSeqNames=c("chr1","chr2")
#' @param WL Windowlength. Length of the initial segmentation of the genome in
#' basepairs. Should be chosen such that on the average 100 reads are contained
#' in each segment. 
#' @param parallel The number of parallel processes to be used for this function.
#' Default=0.
#' @param ... Additional parameters passed to the function "countBamInGRanges" 
#' of the Bioconductor package "exomeCopy". Quality filters for read counts can 
#' be adjusted there. Please see "??countBamInGRanges" for more information.
#' @examples 
#' BAMFiles <- list.files(system.file("extdata", package="cn.mops"),pattern=".bam$",full.names=TRUE)
#' bamDataRanges <- getReadCountsFromBAM(BAMFiles,
#' 					sampleNames=paste("Sample",1:3),WL=5000)
#' X <- getReadCountsFromBAM(BAMFiles,
#' 					sampleNames=paste("Sample",1:3),WL=5000,parallel=2)
#' @return An instance of "GRanges", that contains the breakpoints of the 
#' initial segments and the raw read counts that were extracted from the BAM
#' files. This object can be used as input for cn.mops and other CNV detection
#' methods.
#' @importFrom Rsamtools scanBamHeader
#' @importFrom Rsamtools ScanBamParam
#' @importFrom Rsamtools scanBamFlag
#' @importFrom Rsamtools scanBam
#' @importFrom parallel makeCluster
#' @importFrom parallel clusterEvalQ
#' @importFrom parallel parApply
#' @importFrom parallel stopCluster 
#' @importFrom exomeCopy countBamInGRanges
#' @author Guenter Klambauer \email{klambauer@@bioinf.jku.at}
#' @export


getReadCountsFromBAM <- function(BAMFiles,sampleNames,refSeqNames, WL=25000,parallel=0,...){
	
	if (missing(sampleNames)){
		sampleNames <- basename(BAMFiles)	
	}
	
	headerInfo <- Rsamtools::scanBamHeader(BAMFiles)
	
	targets <- lapply(headerInfo,.subset2,1)
	sn <- names(targets[[1]])
	sl <- as.integer(targets[[1]])
	
	
	message(paste("Identified the following reference sequences: ",
					paste(unique(unlist(sn)),collapse=",")   ))
	
	if (missing(refSeqNames)){
		#refSeqNames <- unique(unlist(sn))[1]
		refSeqNames <- unique(unlist(sn))
		message(paste("Missing \"refSeqNames\"! Selecting all identified reference sequence for analysis."))
	} else{
		message("Using ",paste(refSeqNames,collapse=", ")," as reference.")
	}
	
	message("\n PLEASE BE PATIENT... this might take a while. Consider using the parallel version of this function\n")
	
	if (!(all(refSeqNames %in% unique(unlist(sn))))){
		stop("RefSeqNames does not match identified reference sequences.")
	}
	
	if (any(is.na(sn))){
		stop(paste(refSeqNames,"does not appear in header!"))}
	
	if (length(targets)>1){
		for (i in 2:length(targets)){
			if (!(all(sn==names(targets[[i]])) & 
						all(sl==as.integer(targets[[i]]))  )){
				stop(paste("Targets in Header file of ",BAMFiles[i]," are not identical", 
								"to the header of the file",BAMFiles[1],"!"))
			} 
		}
	}
	
	
	nidx <- match(refSeqNames,sn)
	sn <- sn[nidx]  # sequence names
	sl <- sl[nidx]  # sequence lengths
	sl <- as.integer(unique(sl))
	m <- length(sn)  # number of sequences
	
	# assembling GRanges
	
	GR <- GRanges()
	for (i in 1:m){
		if (WL >= sl[i]){
			tmpGr <- GRanges(sn[i],IRanges(0,sl[i]))
			GenomeInfoDb::seqlevels(tmpGr) <- sn
			
		} else {
			if (sl[i] %% WL ==0) sl[i] <- sl[i] - 1
			tmpBrkpts1 <- seq(0,sl[i],WL)+1
			tmpBrkpts2 <- c(seq(WL,sl[i],WL),sl[i])
			tmpGr <- GRanges(rep(sn[i],length(tmpBrkpts1)),IRanges(tmpBrkpts1,tmpBrkpts2)) 
			GenomeInfoDb::seqlevels(tmpGr) <- sn
		}
		GR <- c(GR,tmpGr)
	}
	
	
	
	# counting 
	
	if (parallel==0){
		XL <- list()
		for (i in 1:length(BAMFiles)){
			XL[[i]] <- exomeCopy::countBamInGRanges(bam.file=BAMFiles[i],granges=GR,...)
		}
		
	} else {
		message("Using parallel version of this function.")
		cl <- parallel::makeCluster(as.integer(parallel),type="SOCK")
		parallel::clusterEvalQ(cl,"exomeCopy::countBamInGRanges")
		XL <- parallel::parLapply(cl,BAMFiles,exomeCopy::countBamInGRanges,granges=GR,...)
	}
	
	if (length(BAMFiles)==1){
		X <- as.matrix(unlist(XL),ncol=1)
	} else	{X <- do.call("cbind",XL)}
	colnames(X) <- BAMFiles
	
		
	mode(X) <- "integer"
	
	colnames(X) <- sampleNames
	mcols(GR) <- X
	
	
	GR <- sortSeqlevels(GR)
	GenomeInfoDb::seqlengths(GR) <- targets[[1]][match(GenomeInfoDb::seqlevels(GR),names(targets[[1]]))]
	
	return(GR)
}
