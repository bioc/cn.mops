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
#' be adjusted there. 
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
#' @importFrom GenomicRanges "strand<-"
#' @importFrom parallel makeCluster
#' @importFrom parallel clusterEvalQ
#' @importFrom parallel parApply
#' @importFrom parallel stopCluster 
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
			XL[[i]] <- .countBamInGRanges(bam.file=BAMFiles[i],granges=GR,...)
		}
		
	} else {
		message("Using parallel version of this function.")
		cl <- parallel::makeCluster(as.integer(parallel),type="SOCK")
		parallel::clusterEvalQ(cl,".countBamInGRanges")
		XL <- parallel::parLapply(cl,BAMFiles,.countBamInGRanges,granges=GR,...)
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


#' countBamInGRanges from package exomeCopy
.countBamInGRanges <- function (bam.file, granges, min.mapq = 1, read.width = 1, stranded.start = FALSE, get.width = FALSE, remove.dup = FALSE) {
	rds.counts <- integer(length(granges))
	seq.names <- unique(as.character(GenomicRanges::seqnames(granges)))
	seq.names.in.bam <- names(Rsamtools::scanBamHeader(bam.file)[[1]]$targets)
	for (seq.name in seq.names) {
		if (seq.name %in% seq.names.in.bam) {
			granges.subset <- granges[GenomicRanges::seqnames(granges) == seq.name]
			strand(granges.subset) <- "*"
			scan.what <- c("pos","mapq")
			if (stranded.start | get.width) {
				scan.what <- c(scan.what, "qwidth")
			}
			if (stranded.start | remove.dup) {
				scan.what <- c(scan.what, "strand")
			}
			rds <- Rsamtools::scanBam(bam.file, param = Rsamtools::ScanBamParam(what = scan.what, which = range(granges.subset)))
			if (min.mapq > 0) {
				mapq.test <- rds[[1]]$mapq >= min.mapq & !is.na(rds[[1]]$mapq)
			} else {
				mapq.test <- rep(TRUE,length(rds[[1]]$mapq))
			}
			if (remove.dup) {
				if (get.width) {
					# this check is fast and accurate, assuming that read widths
					# are less than 1e4 bp and positions are less than 1e12 bp
					# the precision for doubles is 16 significant digits
					mapq.test <- mapq.test & !duplicated(rds[[1]]$pos + as.numeric(rds[[1]]$strand)/10 + rds[[1]]$qwidth/1e5)
				} else {
					mapq.test <- mapq.test & !duplicated(rds[[1]]$pos + as.numeric(rds[[1]]$strand)/10)
				}
			}
			if (sum(mapq.test) > 0) {
				if (stranded.start) {
					rds.ranges <- GenomicRanges::GRanges(seq.name, IRanges::IRanges(start = ifelse(rds[[1]]$strand[mapq.test] %in% c("+","*"), rds[[1]]$pos[mapq.test] , rds[[1]]$pos[mapq.test] + rds[[1]]$qwidth[mapq.test] - read.width), end = ifelse(rds[[1]]$strand[mapq.test] %in% c("+","*"), rds[[1]]$pos[mapq.test] + read.width - 1, rds[[1]]$pos[mapq.test] + rds[[1]]$qwidth[mapq.test] - 1)))
				}
				else if (get.width) {
					rds.ranges <- GenomicRanges::GRanges(seq.name, IRanges::IRanges(start = rds[[1]]$pos[mapq.test], width = rds[[1]]$qwidth[mapq.test]))
				}
				else {
					rds.ranges <- GenomicRanges::GRanges(seq.name, IRanges::IRanges(start = rds[[1]]$pos[mapq.test], width = read.width))
				}
				rds.counts.seq.name <- GenomicRanges::countOverlaps(granges.subset, rds.ranges)
				rds.counts[as.logical(GenomicRanges::seqnames(granges) == seq.name)] <- rds.counts.seq.name
			}
			else {
				rds.counts[as.logical(GenomicRanges::seqnames(granges) == seq.name)] <- 0
			}
		}
		else {
			rds.counts[as.logical(GenomicRanges::seqnames(granges) == seq.name)] <- 0
		}
	}
	if (sum(rds.counts) == 0) {
		warning("No reads found with minimum mapping quality")
	}
  rds.counts
}

