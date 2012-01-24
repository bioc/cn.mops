# Copyright (C) 2012 Klambauer Guenter 
# <klambauer@bioinf.jku.at>

#' Generates the read counts from BAM Files for predefined segments. 
#' This is the appropiate choice for exome sequencing data, where the
#' bait regions, target regions or exons are the predefined segments.
#' These counts are necessary for CNV detection methods based
#' on depth of coverage information.
#' Note that the function is much faster, if the BAM files have an index file.
#' The index file is assumed to be in the same folder and have an identical
#' file name except that ".bai" is appended.
#' 
#' @param BAMFiles BAMFiles
#' @param sampleNames The corresponding sample names to the BAM Files. 
#' @param GR A genomic ranges object that contains the genomic coordinates of
#' the segments. 
#' @param mode Possible values are "paired" and "unpaired", whether the mapping 
#' algorithm was using a "paired" or "unpaired" strategy. Default = "unpaired".
#' @examples 
#' BAMFiles <- list.files(system.file("extdata", package="cn.mops"),pattern=".bam$",
#' 	full.names=TRUE)
#' gr <- GRanges(c("20","20"),IRanges(c(60000,70000),c(70000,80000)))
#' bamDataRanges <- getSegmentReadCountsFromBAM(BAMFiles,GR=gr)
#' @return An instance of "GRanges", that contains the breakpoints of the 
#' initial segments and the raw read counts that were extracted from the BAM
#' files. This object can be used as input for cn.mops and other CNV detection
#' methods.
#' @importFrom Rsamtools countBam
#' @author Guenter Klambauer \email{klambauer@@bioinf.jku.at}
#' @export


getSegmentReadCountsFromBAM <- function(BAMFiles,GR,sampleNames,
		mode="unpaired"){
	
	if (!(mode %in% c("paired","unpaired"))){
		stop("Mode parameter must be \"paired\" or \"unpaired\"!")
	}
	if (missing(sampleNames)){
		sampleNames <- as.character(BAMFiles)	
	}
	
	if (missing(GR) | !inherits(GR,"GRanges")){
		stop("You must submit the coordinates as GRanges object.")
	}
	
	
	if (!all(file.exists(paste(BAMFiles,".bai",sep="")))){
		stop("The indices of the BAM files must be present.\n",
				"They are supposed to have the same file name with ",
				".bai appended.")
	} 
	
	message("This may take a couple of minutes per BAM file.",
			"Please be patient.\n\n")
	
	if (mode=="unpaired"){
		param <- Rsamtools::ScanBamParam(Rsamtools::scanBamFlag(
						isPaired = FALSE),which=GR)
	} else {
		param <- Rsamtools::ScanBamParam(Rsamtools::scanBamFlag(
						isPaired = TRUE),which=GR)
	}
	
	X <- matrix(NA,nrow=length(GR),ncol=length(BAMFiles))
	
	for (i in 1:length(BAMFiles)){
		message("Processing ",BAMFiles[i])
		X[,i] <- Rsamtools::countBam(BAMFiles[i],param=param)$records
	}
	
	colnames(X) <- BAMFiles
		
	mode(X) <- "integer"
	values(GR) <- X
	#names(gr@elementMetadata@listData) <- sampleNames
	colnames(elementMetadata(GR)) <- sampleNames
	
	
	return(GR)
}
