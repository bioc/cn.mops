# Copyright (C) 2012 Klambauer Guenter 
# <klambauer@bioinf.jku.at>

#' @title Calculation of read counts from BAM files for predefined segments.
#' @description Generates the read counts from BAM Files for predefined segments. 
#' This is the appropiate choice for exome sequencing data, where the
#' bait regions, target regions or exons are the predefined segments.
#' These counts are necessary for CNV detection methods based
#' on depth of coverage information.
#' 
#' This function can also be run in a parallel version.
#' 
#' @param BAMFiles BAMFiles
#' @param sampleNames The corresponding sample names to the BAM Files. 
#' @param GR A genomic ranges object that contains the genomic coordinates of
#' the segments. 
#' @param parallel The number of parallel processes to be used for this function.
#' Default=0.
#' @param ... Additional parameters passed to the function "countBamInGRanges" 
#' of the Bioconductor package "exomeCopy". Quality filters for read counts can 
#' be adjusted there. Please see "??countBamInGRanges" for more information.
#' @examples 
#' BAMFiles <- list.files(system.file("extdata", package="cn.mops"),pattern=".bam$", full.names=TRUE)
#' gr <- GRanges(c("20","20"),IRanges(c(60000,70000),c(70000,80000)))
#' bamDataRanges <- getSegmentReadCountsFromBAM(BAMFiles,GR=gr)
#' bamDataRanges <- getSegmentReadCountsFromBAM(BAMFiles,GR=gr,parallel=2)
#' @return An instance of "GRanges", that contains the breakpoints of the 
#' initial segments and the raw read counts that were extracted from the BAM
#' files. This object can be used as input for cn.mops and other CNV detection
#' methods.
#' @importFrom parallel makeCluster
#' @importFrom parallel clusterEvalQ
#' @importFrom parallel parApply
#' @importFrom parallel stopCluster 
#' @author Guenter Klambauer \email{klambauer@@bioinf.jku.at}
#' @export


getSegmentReadCountsFromBAM <- function(BAMFiles,GR,sampleNames,parallel=0,...){
	
	if (missing(sampleNames)){
		sampleNames <- basename(BAMFiles)	
	}
	
	if (missing(GR) | !inherits(GR,"GRanges")){
		stop("You must submit the coordinates as GRanges object.")
	}
	
	
	message("This may take a couple of minutes per BAM file.",
			"Please be patient.\n\n")
	
	
	X <- matrix(NA,nrow=length(GR),ncol=length(BAMFiles))
	
	if (parallel==0){
		for (i in 1:length(BAMFiles)){
			message("Processing ",BAMFiles[i])
			X[,i] <- countBamInGRanges(bam.file=BAMFiles[i],granges=GR,...)
			
		}	
	} else {
		message("Using parallel version of this function.")
		cl <- parallel::makeCluster(as.integer(parallel),type="SOCK")
		parallel::clusterEvalQ(cl,"countBamInGRanges")
		XL <- parallel::parLapply(cl,BAMFiles,countBamInGRanges,granges=GR,...)
		
		parallel::stopCluster(cl)
		for (i in 1:length(BAMFiles)){
			X[,i] <- XL[[i]]
		}
		rm("XL")
	}
	
	
	
	colnames(X) <- BAMFiles
	
	mode(X) <- "integer"
	colnames(X) <- sampleNames
	
	mcols(GR) <- X
#names(gr@elementMetadata@listData) <- sampleNames
#IRanges::colnames(IRanges::elementMetadata(GR)) <- sampleNames
#gr <- GRanges(GenomicRanges::seqnames(GR),IRanges::ranges(GR),
#		sampleNames=X)
	
	return(GR)
}
