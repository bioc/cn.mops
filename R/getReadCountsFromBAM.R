# Copyright (C) 2011 Klambauer Guenter 
# <klambauer@bioinf.jku.at>
.countBAM <- function(bamFile,sl,WL,mode,refSeqName,quiet=FALSE){		
	if (!quiet){message("Reading file: ",bamFile)}
	if (mode=="paired"){
		param <- Rsamtools::ScanBamParam(
				Rsamtools::scanBamFlag(isPaired = TRUE,
						isFirstMateRead=TRUE,
						isProperPair=TRUE),
				what=c("rname","pos","mpos"))
				#which=RangesList(refSeqName))
		readPos <-Rsamtools::scanBam(bamFile,param=param)[[1]]
		readPosIdx <- (readPos$rname==refSeqName)
		readPos$pos <- readPos$pos[readPosIdx]
		readPos$mpos <- readPos$mpos[readPosIdx]
		pp <- ((readPos$pos+readPos$mpos)/2)
		gapSize <- abs(readPos$pos-readPos$mpos)
		rm("readPos")
		ppFiltered <- pp[!(gapSize>5*median(gapSize,na.rm=TRUE)|
							is.na(gapSize))]
		ppFiltered <- ppFiltered[!is.na(ppFiltered)]
	} else{
		param <- Rsamtools::ScanBamParam(Rsamtools::scanBamFlag(
						isPaired = FALSE),
				what=c("rname","pos"))
				#which=RangesList(refSeqName))
		readPos <- Rsamtools::scanBam(bamFile,param=param)[[1]]
		readPosIdx <- (readPos$rname==refSeqName)
		readPos$pos <- readPos$pos[readPosIdx] 
		ppFiltered <- readPos$pos[!is.na(readPos$pos)]
		rm("readPos")
	}
	if (length(ppFiltered)==0){
		warning("No reads found in file: ",bamFile)
	}
	if (sl%%WL!=1){
		brkpts <- c(seq(1,sl,WL),sl)
	} else{
		brkpts <- seq(1,sl,WL)
	}
	if (any(ppFiltered>=sl)){
		warning(paste("Some read positions are greater than",
						"length of reference sequence! File: ",bamFile,"\n"))
		ppFiltered <- ppFiltered[ppFiltered<sl]
	}
	
	x <- hist(ppFiltered,breaks=brkpts,plot=FALSE)$counts

	return(x)
}


#' Generates the read counts from BAM Files. 
#' These counts are necessary for CNV detection methods based
#' on depth of coverage information.
#' 
#' 
#' @param BAMFiles BAMFiles
#' @param sampleNames The corresponding sample names to the BAM Files. 
#' @param refSeqName Name of the reference sequence that should be analyzed.
#' The name must appear in the header of the BAM file. If it is not given
#' the function will select the first reference sequence that appears in the
#' header of the BAM files.
#' @param WL Windowlength. Length of the initial segmentation of the genome in
#' basepairs. Should be chosen such that on the average 100 reads are contained
#' in each segment. If not given, cn.mops will try to find an appropiate window 
#' length.
#' @param mode Possible values are "paired" and "unpaired", whether the mapping 
#' algorithm was using a "paired" or "unpaired" strategy. Default = "unpaired".
#' @examples 
#' BAMFiles <- list.files(system.file("extdata", package="cn.mops"),pattern=".bam",
#' 	full.names=TRUE)
#' bamDataRanges <- getReadCountsFromBAM(BAMFiles,
#' 					sampleNames=paste("Sample",1:3),WL=5000)
#' @return An instance of "GRanges", that contains the breakpoints of the 
#' initial segments and the raw read counts that were extracted from the BAM
#' files. This object can be used as input for cn.mops and other CNV detection
#' methods.
#' @importFrom Rsamtools scanBamHeader
#' @importFrom Rsamtools ScanBamParam
#' @importFrom Rsamtools scanBamFlag
#' @importFrom Rsamtools scanBam
#' @author Guenter Klambauer \email{klambauer@@bioinf.jku.at}
#' @export


getReadCountsFromBAM <- function(BAMFiles,sampleNames,refSeqName,WL,
		mode="unpaired"){

	if (!(mode %in% c("paired","unpaired"))){
		stop("Mode parameter muste be \"paired\" or \"unpaired\"!")
	}
	if (missing(sampleNames)){
		sampleNames <- as.character(BAMFiles)	
	}
	
	headerInfo <- Rsamtools::scanBamHeader(BAMFiles)
	sn <- lapply(headerInfo,function(h){ 
				sn <- sapply(h$text,.subset2,1)				
				sn <- sn[names(sn)=="@SQ"]
				sn <- sapply(strsplit(sn,"SN:"),.subset2,2)
				return(sn)
			}
	)
	message(paste("Identified the following reference sequences: ",
					paste(unique(unlist(sn)),collapse=",")   ))
	if (!all(table(unlist(sn))==length(BAMFiles))){
		stop("Headers are not identical in the given BAM files.")
	}
	if (missing(refSeqName)){
		refSeqName <- unique(unlist(sn))[1]
		message(paste("Missing \"refSeqName\"! Selecting",refSeqName,
						"as reference sequence."))
	} else{
		message("Using ",refSeqName," as reference.")
	}
	
	if (!(all(refSeqName %in% unique(unlist(sn))))){
		stop("RefSeqName does not match identified reference sequences.")
	}
	
	ln <- sapply(headerInfo,function(h){ 
				sn <- sapply(h$text,.subset2,1)
				ln <- sapply(h$text,.subset2,2)
				ln <- ln[grep(paste("SN:",refSeqName,"$",sep=""),sn)]
				if (length(ln)!=1){return(NA)} else{return(ln)}
			}
	)
	
	if (any(is.na(sn))){stop(paste(refSeqName,"does not appear in header!"))}
	
	sl <- sapply(strsplit(ln,"LN:"),.subset2,2)
	
	if (length(unique(sn))!=1 | length(unique(sl))!=1){
		stop("\"refSeqName\" not consistent in given bamFiles!")
	}
	sl <- as.integer(unique(sl))
	
	if (missing(WL)){
		message(paste("Missing \"WL\"! cn.mops will suggest an",
						"appropiate value for the window length."))
		BAMFiles <- BAMFiles[order(file.info(BAMFiles)$size)]
		xs <- sum(.countBAM(BAMFiles[1],sl=sl,WL=100000,mode=mode,
						refSeqName=refSeqName,quiet=TRUE))
		if (xs==0){
			message(paste("It is not possible calculate an appropiate",
							"windowlength."))
			WL <- 25000
		} else{
			WL <- max(round(100*sl/(xs),-3),100)
		}
		message("Window length set to: ",WL)
	}
	
	XL <- lapply(BAMFiles,.countBAM,sl=sl,WL=WL,mode=mode,
			refSeqName=refSeqName)
	
	if (sl%%WL!=1){
		brkpts <- c(seq(1,sl,WL),sl)
	} else{
		brkpts <- seq(1,sl,WL)
	}
#	browser()

	nSegm <- length(brkpts)
	if (length(BAMFiles)==1){
		X <- as.matrix(unlist(XL),ncol=1)
	} else	{X <- do.call("cbind",XL)}
	
	colnames(X) <- BAMFiles
	rownames(X) <- paste(refSeqName,"_",brkpts[1:(nSegm-1)],"_",
			brkpts[2:(nSegm)]-1,sep="")
	
	ir <- IRanges(start=brkpts[1:(nSegm-1)],
			end=brkpts[2:(nSegm)]-1)
	
	
	gr <- GenomicRanges::GRanges(seqnames=refSeqName, ranges = ir)
	mode(X) <- "integer"
	values(gr) <- X
	#names(gr@elementMetadata@listData) <- sampleNames
	colnames(elementMetadata(gr)) <- sampleNames
	
	
	return(gr)
}
