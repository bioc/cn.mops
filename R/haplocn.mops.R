# Copyright (C) 2011 Klambauer Guenter 
# <klambauer@bioinf.jku.at>

haplocn.mopsC <- function(x,I = c(0.025,1,2,3,4,5,6,7,8), 
		classes=c("CN0","CN1","CN2","CN3","CN4","CN5","CN6","CN7","CN8"), cov,
		priorimpact = 1,cyc = 20) {
	
	version <- packageDescription("cn.mops")$Version
	
	N <- length(x)
	n <- length(I)
	
	if (missing(cov)){cov <- rep(1,N)}
	
	idxCN1 <- which(classes=="CN1")
	
	alpha.init <- rep(0.05,n)
	alpha.init[idxCN1] <- 0.6
	alpha.init <- alpha.init/ sum(alpha.init)
	alpha.est <- alpha.init
	alpha.prior <- rep(0,n)
	alpha.prior[idxCN1] <- 1
	alpha.prior <- alpha.prior*priorimpact
	
	if (all(x<=1)) {
		lambda.est <- rep(0,n)
		alpha.est <- rep(0,n)
		alpha.est[idxCN1] <- 1
		expCN <- rep(classes[idxCN1],N)
		ini <- 0
		ExpLogFoldChange <- rep(0,N)
		post.ik <- matrix(0,nrow=n,ncol=N)
		post.ik[idxCN1, ] <- 1
		
		
		params <- list(n,classes,I,priorimpact,cyc)
		names(params) <- c("nclasses","classes","I","priorimpact","cyc")
		l <-  list ("lambda"=lambda.est, "alpha"=alpha.est, "expectedCN"=expCN,
				"sini"=ExpLogFoldChange,"ini"=ini,"post"=post.ik, 
				"params"=params,"version"=version)
		return(l)
		
	} else {
		lambda.est <- median(x*1/cov,na.rm=TRUE)
		if (lambda.est < 1e-10){lambda.est <- max(mean(x*1/cov,na.rm=TRUE),1.0)}
		lambda.init <- I*lambda.est
		ret=.Call("cnmops", as.numeric(x), I,cov, as.integer(cyc), alpha.init, lambda.init,
				alpha.prior)
		alpha.ik=ret$alpha.ik
		alpha.i=ret$alpha.i
		alpha.est=ret$alpha.est
		lambda.est=ret$lambda.est
		
		#Posterior
		post.ik <- alpha.ik
		if (is.null(names(x))){
			colnames(post.ik) <- paste("x_",1:N,sep="")
		} else {
			colnames(post.ik) <- names(x)
		}
		rownames(post.ik) <- classes
		expCN <- classes[apply(post.ik,2,function(x) which(x==max(x))[1] )]
		
		ini <- mean(abs(log2(I)) %*% alpha.ik)
		ExpLogFoldChange <-  log2(I) %*%  post.ik
		params <- list(n,classes,I,priorimpact,cyc)
		names(params) <- c("nclasses","classes","I","priorimpact","cyc")
		l <-  list ("lambda"=lambda.est, "alpha"=alpha.est, "expectedCN"=expCN, 
				"sini"=ExpLogFoldChange, "ini"=ini, "post"=post.ik,
				"params"=params,"version"=version)
		return(l)
	}
}


#' Performs the cn.mops algorithm for copy number detection in
#' NGS data adjusted to haploid genomes. It is assumed that the normal state
#' is copy number 1. This is an experimental method at the moment.
#' 
#' @param input Either an instance of "GRanges" or a raw data matrix, where
#' columns are interpreted as samples and rows as genomic regions. An entry is
#' the read count of a sample in the genomic region.
#' @param I Vector positive real values that contain the expected fold change
#' of the copy number classes.  Length of this vector must be equal to the 
#' length of the "classes" parameter vector. For copy number polymorphisms 
#' in haploid organisms we suggest to use the default 
#' I = c(0.025,1,2,3,4,5,6,7,8).
#' @param classes Vector of characters of the same length as the parameter
#' vector "I". One vector element must be named "CN1". The names reflect the 
#' labels of the copy number classes. 
#' Default = c("CN0","CN1","CN2","CN3","CN4","CN5","CN6","CN7","CN8").
#' @param priorImpact Positive real value that reflects how strong the prior
#' assumption affects the result. The higher the value the more samples will
#' be assumed to have copy number 1. Default = 1.
#' @param cyc Positive integer that sets the number of cycles for the algorithm.
#' Usually after less than 15 cycles convergence is reached. Default = 20.
#' @param parallel How many cores are used for the computation. If set to zero
#' than no parallelization is applied. The package "snow" has to be installed
#' for this option. Default = 0.
#' @param normType Mode of the normalization technique. Possible values are 
#' "mean","min","median","quant", "poisson" and "mode". 
#' Read counts will be scaled sample-wise. Default = "poisson".
#' @param normQu Real value between 0 and 1.  
#' If the "normType" parameter is set to "quant" then this parameter sets the 
#' quantile that is used for the normalization. Default = 0.25. 
#' @param norm Logical that indicates whether normalization should be 
#' applied or not. Default = TRUE.
#' @param upperThreshold Positive real value that sets the cut-off for copy
#' number gains. All CNV calling values above this value will be called as 
#' "gain". The value should be set close to the log2 of the expected foldchange
#' for copy number 3 or 4. Default = 0.5.
#' @param lowerThreshold Negative real value that sets the cut-off for copy
#' number losses. All CNV calling values below this value will be called as 
#' "loss". The value should be set close to the log2 of the expected foldchange
#' for copy number 1 or 0. Default = -0.9.
#' @param minWidth Positive integer that is exactly the parameter "min.width"
#' of the "segment" function of "DNAcopy". minWidth is the minimum number 
#' of segments a CNV should span. Default = 4.
#' @param segAlgorithm Which segmentation algorithm should be used. If set to
#' "DNAcopy" circular binary segmentation is performed. Any other value will
#' initiate the use of our fast segmentation algorithm. Default = "fast".
#' @param ... Additional parameters will be passed to the "DNAcopy"
#' or the standard segmentation algorithm.
#' @examples 
#' data(cn.mops)
#' haplocn.mops(XRanges[1:200, ])
#'
#' @useDynLib cn.mops
#' @return An instance of "CNVDetectionResult".
#' @author Guenter Klambauer \email{klambauer@@bioinf.jku.at}
#' @export


haplocn.mops <- function(input,I = c(0.025,1,2,3,4,5,6,7,8),
		classes=c("CN0","CN1","CN2","CN3","CN4","CN5","CN6","CN7","CN8"),
		priorImpact = 1,cyc = 20,parallel=0,
		normType="poisson",normQu=0.25,norm=TRUE,
		upperThreshold=0.6,lowerThreshold=-0.9,
		minWidth=3,segAlgorithm="fast",...){
	
	version <- packageDescription("cn.mops")$Version
	
	
	if(class(input)=="GRanges"){
		inputType <- "GRanges"
		input <- IRanges::sort(input)
		#X <- (do.call("cbind",(values(input)@unlistData@listData)))
		#X <- do.call("cbind",values(input)@listData)
		#X <- matrix(as.numeric(X),nrow=nrow(X))
		X <- IRanges::as.matrix(IRanges::values(input))
		
		if (ncol(X)==1){
			stop("It is not possible to run cn.mops on only ONE sample.\n")
		}
		chr <- as.character(seqnames(input))
		start <- start(input)
		end <- end(input)
		
		irAllRegions <- IRanges(start,end)
		names(irAllRegions) <- NULL
		
		#maploc <- (start+end)/2
	} else if (is.matrix(input)){
		if (nrow(input)> 1){
			inputType <- "DataMatrix"
			X <- input
			X <- matrix(as.numeric(X),nrow=nrow(X))
			colnames(X) <- colnames(input)	
			chr <- rep("undef",nrow(X))
			irAllRegions <- IRanges(start=1:nrow(X),end=1:nrow(X))
		} else{
			inputType <- "DataMatrix"
			chr <- "undef"
			irAllRegions <- IRanges(start=1:nrow(X),end=1:nrow(X))
			parallel <- 0
		}
	} else if (is.vector(input)) {
		inputType <- "DataMatrix"
		X <- matrix(input,nrow=1)
		X <- matrix(as.numeric(X),nrow=nrow(X))
		chr <- "undef"
		irAllRegions <- IRanges(start=1:nrow(X),end=1:nrow(X))
		parallel <- 0
	}else{
		stop("GRanges object or read count matrix needed as input.")
	}
	
	
	if (any(X<0) | any(!is.finite(X))){
		stop("All values must be greater or equal zero and finite.\n")
	}	
	if (length(I)!=length(classes)){
		stop("I and classes must have same length!")
	}
	if (!("CN1" %in% classes)){stop("One element of classes must be CN1 .\n")}
	
	m <- nrow(X)
	N <- ncol(X)
	n <- length(I)
	
	if (is.null(colnames(X))){
		colnames(X) <- paste("Sample",1:N,sep="_")
	}
	
	
	params <- list("cn.mops",I,
			classes,
			priorImpact,cyc,
			normType,normQu,
			upperThreshold,lowerThreshold,
			minWidth,paste(...))
	names(params) <- c("method","folds",
			"classes",
			"priorimpact","cyc",
			"normType","normQu",
			"upperThreshold","lowerThreshold",
			"minWidth","SegmentationParams")
	
	if (m < 100){
		warning(paste("For this small number of segments normalization",
						"might not be applicable."))
	}
	
	if (norm) {
		message("Normalizing...")
		X.norm <- normalizeChromosomes(X,chr=chr,normType=normType,qu=normQu)
	} else {
		#X.norm <- normalizeChromosomes(X,chr=chr,normType=normType,qu=normQu)
		X.norm <- X
	}
	message("Starting local modeling, please be patient...  ")
	
	if (parallel > 0){
		library(snow)
		cl <- makeCluster(as.integer(parallel),type="SOCK")
		clusterEvalQ(cl,"haplocn.mopsC")
	} 
	
	res <- list()
	Xchr <- list()
	for (chrom in unique(chr)){
		message(paste("Reference sequence: ",chrom))
		chrIdx <- which(chr==chrom)
		if (norm) {Xchr[[chrom]] <- X.norm[chrIdx, ]} else {
			Xchr[[chrom]] <- X[chrIdx, ]
		}
		
#		if(m > 1 & !norm){
#			cov <- colSums(X[chrIdx, ])
#			if (median(cov) > 0){
#				cov <- cov/median(cov)
#			} else {
#				stop(paste("Median of total reads is zero.",
#								"Too many samples with zero reads."))
#			}
#		} else {
#			#X.norm <- X
#			cov <- rep(1,N)
#		}
		cov <- rep(1,N)
		
		#cat("Coverage: ",cov ," Norm: ",norm)
		
		#cov[which(cov < 1e-15)] <- 1e-15
		
		if (norm & m > 1){
			if (parallel==0){
				resChr <-apply(X.norm[chrIdx, ,drop=FALSE],1,haplocn.mopsC, I=I,
						classes=classes,
						cov=cov,priorimpact=priorImpact,cyc=cyc)
			} else {
				resChr <- parApply(cl,X.norm[chrIdx, ,drop=FALSE],1,haplocn.mopsC, 
						I=I,classes=classes,cov=cov,priorimpact=priorImpact,
						cyc=cyc)				
			}
		} else {
			if (parallel==0){
				resChr <-apply(X[chrIdx, ,drop=FALSE],1,haplocn.mopsC, I=I,
						classes=classes,
						cov=cov,priorimpact=priorImpact,cyc=cyc)
			} else {
				resChr <- parApply(cl,X[chrIdx, ,drop=FALSE],1,haplocn.mopsC, 
						I=I,classes=classes,cov=cov,priorimpact=priorImpact,
						cyc=cyc)				
			}
		}
		
		#cat(paste(Sys.time(),"\n"))
		res <- c(res, resChr)
	}
	if (parallel > 0){
		stopCluster(cl)
	} 
	
	message("Postprocessing result...")
	
	## Postprocess result
	L <- t(sapply(res,.subset2,1))
	rownames(L) <- rownames(X)
	colnames(L) <- classes
	params$L <- L
	A <- t(sapply(res,.subset2,2))
	rownames(A) <- rownames(X)
	colnames(A) <- classes
	params$A <- A
	CN <- t(sapply(res,.subset2,3))
	rownames(CN) <- rownames(X)
	colnames(CN) <- colnames(X)
	sINI <- t(sapply(res,.subset2,4))
	rownames(sINI) <- rownames(X)
	colnames(sINI) <- colnames(X)
	#write.table(sINI,file="sINIbefore.txt")
	INI <- (sapply(res,.subset2,5))
	names(INI) <- rownames(X)
	post <- array(dim=c(m,n,N))
	post.tmp <- t(lapply(res,.subset2,6))
	for (i in 1:m){
		post[i, ,] <- post.tmp[[i]]
	}
	dimnames(post) <- list(NULL,classes,colnames(X))
	rm("post.tmp")	
	
	
	if (m>1){
		message("Starting segmentation algorithm...")
		
		if (segAlgorithm=="DNAcopy"){
			message("Using \"DNAcopy\" for segmentation.")
			library(DNAcopy)
			if (!exists("eta")){eta <- 0.05}
			if (!exists("nperm")){nperm <- 10000}
			if (!exists("alpha")){alpha <- 0.01}
			if (minWidth > 5){
				message("For DNAcopy the maximum \"minWidth\" is 5.")
				message("Resetting \"minWidth\" to 5.")
				minWidth <- 5
			}
			DNAcopyBdry <- DNAcopy::getbdry(eta=eta,nperm=nperm,tol=alpha,
					max.ones=floor(nperm*alpha)+1)
			
			
			if (parallel==0){
				resSegm <- apply(sINI,2,.segmentation,
						chr=chr,minWidth=minWidth,DNAcopyBdry=DNAcopyBdry,...)
			} else {
				cl <- makeCluster(as.integer(parallel),type="SOCK")
				clusterEvalQ(cl,".segmentation")
				resSegm <- parApply(cl,sINI,2,.segmentation,
						chr=chr,minWidth=minWidth,DNAcopyBdry=DNAcopyBdry,...)
				stopCluster(cl)
			}
			
		
			
			
			segDf <- cbind(do.call(rbind,resSegm),
					rep(colnames(X),sapply(resSegm,nrow)))
			rm("resSegm")
			
			segDf <- data.frame("chr"=as.character(segDf[,1]),"from"=segDf[,2],
					"to"=segDf[,3],"value"=segDf[,4],"sample"=segDf[,5],
					stringsAsFactors=FALSE)
			#segDf <- segDf[order(segDf$chr,segDf$sample,segDf$from), ]
			
			#colnames(segDf) <- c("chr","from","to","value","sample")
			#segDf$chr <- as.character(segDf$chr)
			
			segCN <- apply(segDf,1,function(x){
						sIdx <- as.integer(x["from"]):as.integer(x["to"])
						if (length(sIdx)>=3){
							sIdx2 <- sIdx[-c(1,length(sIdx))]
						} else{
							sIdx2 <- sIdx
						}
						CN <- haplocn.mopsC(apply(Xchr[[x["chr"]]][sIdx2,
												,drop=FALSE],2,
										median),
								I=I,
								classes=classes,
								cov=cov,priorimpact=priorImpact,
								cyc=cyc)$expectedCN
						
						#Median
						#segValue <- quantile(sINI[sIdx,x["sample"]],probs=0.5,
						#		na.rm=TRUE)
						
						#Mean
						#segValue <- mean(sINI[sIdx,x["sample"]],na.rm=TRUE)
						#segValue2 <- median(sINI[sIdx,x["sample"]],na.rm=TRUE)
						#segValue <- x["value"]
						#segValue2 <- x["value"]
						
						return(CN[match(x["sample"],colnames(X))])
						
					})
			
			
			## value <- sapply(segCN,.subset2,2)
			## value2 <- sapply(segCN,.subset2,3)
			## segDf <- data.frame(segDf,"CN"=sapply(segCN,.subset2,1),
			##         "mean"=value,"median"=value2,
			##         stringsAsFactors=FALSE)
			segDf$CN <- segCN
			
			# we do not compute the seg median:
			segDf$value2 <- segDf$value 
			segDf <- segDf[,c("chr","from","to","sample","value","value2","CN")]
			colnames(segDf) <-
					c("chr","from","to","sample","mean","median","CN")
			
			callsS <- matrix(NA,nrow=m,ncol=N)
			for (chrom in unique(chr)){
				chrIdx <- which(chr==chrom)
				segDfTmp <- subset(segDf,chr==chrom)
				callsS[chrIdx, ] <- 
						matrix(rep(segDfTmp$mean,segDfTmp$to-segDfTmp$from+1),
								ncol=N)
			}
			
			colnames(callsS) <- colnames(X)
			
			
			segDfSubset <- segDf[which(
							segDf$mean >= upperThreshold |
									segDf$mean <= lowerThreshold), ]
			segDfSubset <- segDfSubset[which(
							(segDfSubset$to-segDfSubset$from+1) >= minWidth), ]
			
			
			
		} else {
			message("Using \"fastseg\" for segmentation.")
			resSegmList <- list()
			segDf <- data.frame(stringsAsFactors=FALSE)
			if (!exists("segMedianT")){
				segMedianT <- 0.01
			} else {
				stop(paste("Use \"upper-\" and \"lowerThreshold\" parameters",
								"instead of segMedianT."))
			}
			
			#if (!exists("alpha")){alpha <- 0.1}
			callsS <- matrix(NA,nrow=m,ncol=N)
			colnames(callsS) <- colnames(X)
			for (chrom in unique(chr)){
				chrIdx <- which(chr==chrom)
				
				if (parallel==0){
					resSegmList[[chrom]] <- apply(sINI[chrIdx, ],2,
							cn.mops:::segment,
							minSeg=minWidth,segMedianT=segMedianT,
							segPlot=FALSE,...)
				} else {
					cl <- makeCluster(as.integer(parallel),type="SOCK")
					clusterEvalQ(cl,"segment")
					resSegmList[[chrom]] <- parApply(cl,sINI[chrIdx, ],2,
							segment,minSeg=minWidth, segMedianT=segMedianT,
							segPlot=FALSE,...)
					stopCluster(cl)
				}
				
				segDfTmp <- cbind(do.call(rbind,resSegmList[[chrom]]),
						"sample"=rep(colnames(X),
								sapply(resSegmList[[chrom]],nrow)))
				segDfTmp$chr <- chrom
				
				
				callsS[chrIdx, ] <- 
						matrix(rep(segDfTmp$mean,segDfTmp$end-segDfTmp$start+1),
								ncol=N)
				
				segDf <- rbind(segDf,segDfTmp)
			}
			
			
			#segDf <- segDf[,c("chr","start","end","sample")]
			#colnames(segDf) <- c("chr","from","to","sample")
			
			message("Calculating integer copy numbers...")
			#browser()
			
			segCN <- apply(segDf[,c("chr","start","end","sample")],1,
					function(x){
						
						sIdx <- as.integer(x["start"]):as.integer(x["end"])
						if (length(sIdx)>=3){
							sIdx <- sIdx[-c(1,length(sIdx))]
						}
						
						CN <- haplocn.mopsC(apply(Xchr[[x["chr"]]][sIdx,
												,drop=FALSE],2,
										median),
								I=I,
								classes=classes,
								cov=cov,priorimpact=priorImpact,
								cyc=cyc)$expectedCN
						
#						cnProbs <- post[sIdx , , x["sample"]]
#						if (length(sIdx)>1){
#							cnProbs <- colSums(-log(cnProbs))
#						} else {
#							cnProbs <- -log(cnProbs)
#						}
#						mIdx <- which(cnProbs==min(cnProbs))
#						
#						return(paste(classes[mIdx],collapse="/"))
						return(CN[which(x["sample"]==colnames(X))])
					})
			
			segDf <- data.frame(segDf,"CN"=segCN,stringsAsFactors=FALSE)
		
			
			colnames(segDf) <- c("from","to","mean","median","sample",
					"chr","CN")
			
			#value <- segDf$value
			
			#median
			#callsS <- matrix(rep(segDf$value,segDf$to-segDf$from+1),ncol=N)
			#segDf <- segDf[which(segDf$value >= upperThreshold
			#						| segDf$value <= lowerThreshold), ]
			
			#mean for AUC
			
			
			
			#mean for the table
			segDfSubset <- segDf[which(segDf$mean >= upperThreshold
									| segDf$mean <= lowerThreshold), ]
			
			
			#segDf$value <- segDf$value2
			
			
			segDfSubset <- 
					segDfSubset[which((segDfSubset$to-segDfSubset$from+1)
											>= minWidth), ]
			
			
			
		}
		
		
		
		
		if (nrow(segDfSubset)>0){
			
			# Assembly of result object
			r <- new("CNVDetectionResult")
			cnvrR <- reduce(GRanges(seqnames=segDfSubset$chr,
							IRanges(segDfSubset$from,segDfSubset$to)))
			cnvrCN <- matrix(NA,ncol=N,nrow=length(cnvrR))
			for (jj in 1:length(cnvrR)){
				sIdx3 <- (start(cnvrR)[jj]:end(cnvrR)[jj])
				chrIdx2 <- as.character(seqnames(cnvrR))[jj]
				if (length(sIdx3) >= 3){sIdx3 <- sIdx3[2:(length(sIdx3)-1)]}
				cnvrCN[jj, ] <- haplocn.mopsC(apply(Xchr[[chrIdx2]][sIdx3,
										,drop=FALSE],
								2,median),
						I=I,
						classes=classes,
						cov=cov,priorimpact=priorImpact,
						cyc=cyc)$expectedCN
				
				
			}
			colnames(cnvrCN) <- colnames(X) 
			
			sampleNames <- segDfSubset$sample
			
			if (inputType=="GRanges"){
				ir <- IRanges()
				irCNVR <- IRanges()
				for (chrom in unique(chr)){
					inputChr <- input[which(as.character(
											seqnames(input))==chrom)]
					segDfSubsetChr <- subset(segDfSubset,chr==chrom)
					cnvrRChr <- cnvrR[which(as.character(
											seqnames(cnvrR))==chrom)]
					if (nrow(segDfSubsetChr) >0){
						ir <- c(ir,IRanges(start(inputChr)[
												segDfSubsetChr$from],
										end(inputChr)[segDfSubsetChr$to]))
						
						irCNVR <- c(irCNVR,IRanges(start(inputChr)[
												start(cnvrRChr)],
										end(inputChr)[end(cnvrRChr)]))
					}
				}
			} else if (inputType=="DataMatrix"){
				ir <- IRanges(start=segDfSubset$from,end=segDfSubset$to)
				irCNVR <- IRanges(start=start(cnvrR),end=end(cnvrR))
			}
			
			
			rd <- GRanges(seqnames=segDfSubset$chr,ir,"sampleName"=sampleNames,
					"median"=segDfSubset$median,"mean"=segDfSubset$mean,
					"CN"=segDfSubset$CN)
			
			#cnvr <- reduce(rd)
			
			cnvr <- GRanges(seqnames=seqnames(cnvrR),irCNVR,CN=cnvrCN)
			#values(cnvrR) <- cnvrCN
			#colnames(elementMetadata(cnvrR)) <- colnames(X)
			
			
			r@normalizedData    <- X.norm
			r@localAssessments  <- sINI
			#write.table(sINI,file="sINIafter.txt")
			
			r@individualCall   	<- callsS
			r@iniCall        	<- INI
			r@cnvs				<- rd
			r@cnvr				<- cnvr
			
			
			
			if (inputType=="GRanges"){
				irS <- IRanges()
				for (chrom in unique(chr)){
					inputChr <- input[which(as.character(
											seqnames(input))==chrom)]
					
					segDfChr <- subset(segDf,chr==chrom)
					if (nrow(segDfChr) >0 ){
						irS <- c(irS,IRanges(start(inputChr)[segDfChr$from],
										end(inputChr)[segDfChr$to]))
					}
				}
				r@segmentation 			<- GRanges(seqnames=segDf$chr,
						irS,
						"sampleName"=segDf$sample,"median"=segDf$median,
						"mean"=segDf$mean,"CN"=segDf$CN)
			} else if (inputType=="DataMatrix"){
				r@segmentation 			<- GRanges(seqnames=segDf$chr,
						IRanges(segDf$from,segDf$to),"sampleName"=segDf$sample,
						"median"=segDf$median,
						"mean"=segDf$mean,"CN"=segDf$CN)
			}
			
			r@gr <- GRanges(seqnames=chr,irAllRegions)
			r@posteriorProbs 	<- post
			r@params			<- params
			r@integerCopyNumber	<- CN
			r@sampleNames		<- colnames(X)
			
			return(r)	
		} else {
			message(paste("No CNVs detected. Try changing \"normalization\",", 
							"\"priorimpact\" or \"thresholds\"."))
			# Assembly of result object
			r <- new("CNVDetectionResult")
			if (inputType=="GRanges"){
				irS <- IRanges()
				for (chrom in unique(chr)){
					inputChr <- input[which(as.character(
											seqnames(input))==chrom)]
					
					segDfChr <- subset(segDf,chr==chrom)
					if (nrow(segDfChr) >0 ){
						irS <- c(irS,IRanges(start(inputChr)[segDfChr$from],
										end(inputChr)[segDfChr$to]))
					}
				}
				r@segmentation 			<- GRanges(seqnames=segDf$chr,
						irS,
						"sampleName"=segDf$sample,"median"=segDf$median,
						"mean"=segDf$mean,"CN"=segDf$CN)
			} else if (inputType=="DataMatrix"){
				r@segmentation 			<- GRanges(seqnames=segDf$chr,
						IRanges(segDf$from,segDf$to),
						"sampleName"=segDf$sample,"median"=segDf$median,
						"mean"=segDf$mean,"CN"=segDf$CN)
			}
			
			r@gr <- GRanges(seqnames=chr,irAllRegions)
			r@normalizedData    <- X.norm
			r@localAssessments  <- sINI
			r@gr <- GRanges(seqnames=chr,irAllRegions)
			r@individualCall   	<- callsS
			r@iniCall        	<- INI
			r@posteriorProbs 	<- post
			r@params			<- params
			r@integerCopyNumber	<- CN
			r@sampleNames		<- colnames(X)
			return(r)	
		}
		
		
	} else {
		message(paste("Only one genomic segment considered, therefore no ",
						"segmentation."))	
		# Assembly of result object
		r <- new("CNVDetectionResult")	#
		
		r@gr <- GRanges(seqnames=chr,irAllRegions)
		r@normalizedData    <- X.norm
		r@localAssessments  <- sINI
		r@individualCall   	<- sINI
		
		r@params			<- params
		
		r@integerCopyNumber	<- CN
		r@sampleNames		<- colnames(X)
		
		
		return(r)	
		
	}
}
