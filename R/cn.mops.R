# Copyright (C) 2011 Klambauer Guenter 
# <klambauer@bioinf.jku.at>

.cn.mopsC <- function(x,I = c(0.025,0.5,1,1.5,2,2.5,3,3.5,4), 
		classes=c("CN0","CN1","CN2","CN3","CN4","CN5","CN6","CN7","CN8"), cov,
		priorimpact = 1,cyc = 15) {
	
	version <- "0.99"
	
	N <- length(x)
	n <- length(I)
	
	if (missing(cov)){cov <- rep(1,N)}
	
	idxCN2 <- which(classes=="CN2")
	
	alpha.init <- rep(0.05,n)
	alpha.init[idxCN2] <- 0.6
	alpha.init <- alpha.init/ sum(alpha.init)
	alpha.est <- alpha.init
	alpha.prior <- rep(0,n)
	alpha.prior[idxCN2] <- 1
	alpha.prior <- alpha.prior*priorimpact
	
	if (all(x<=10)) {
		lambda.est <- rep(0,n)
		alpha.est <- rep(0,n)
		alpha.est[idxCN2] <- 1
		expCN <- rep(classes[idxCN2],N)
		ini <- 0
		ExpLogFoldChange <- rep(0,N)
		post.ik <- matrix(0,nrow=N,ncol=n)
		post.ik[idxCN2, ] <- 1
		
		
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
		ret=.Call("cnmops", x, I,cov, as.integer(cyc), alpha.init, lambda.init,
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

.segmentation <- function(x,chr,minWidth,DNAcopyBdry,...){
	m <- length(x)
	xx <- x+rnorm(mean=0,sd=0.00001,n=m)
	
	if (length(chr)!=m){
		stop("Vector \"chr\" must have the same length as \"x\"")
	}
	
	CNA.object <- DNAcopy::CNA(xx,
			chr,1:m,
			data.type="logratio")
	
	
	segment.CNA.object <- DNAcopy::segment(CNA.object,
			min.width=minWidth, sbdry=DNAcopyBdry,...)
	
	segDf <- segment.CNA.object$output
	names(segDf) <- c("sample","chr","from","to","idx","value")
	
	return(segDf[,c("chr","from","to")])
}

#' Performs the cn.mops algorithm for copy number detection in
#' NGS data.
#' 
#' @param input Either an instance of "GRanges" or a raw data matrix, where
#' columns are interpreted as samples and rows as genomic regions. An entry is
#' the read count of a sample in the genomic region.
#' @param I Vector positive real values that contain the expected fold change
#' of the copy number classes.  Length of this vector must be equal to the 
#' length of the "classes" parameter vector. For human copy number polymorphisms 
#' we suggest to use the default I = c(0.05,0.5,1,1.5,2,2.5,3,3.5,4).
#' @param classes Vector of characters of the same length as the parameter
#' vector "I". One vector element must be names "CN2". The names reflect the 
#' labels of the copy number classes. 
#' Default = c("CN0","CN1","CN2","CN3","CN4","CN5","CN6","CN7","CN8").
#' @param priorImpact Positive real value that reflects how strong the prior
#' assumption affects the result. The higher the value the more samples will
#' be assumed to have copy number 2. Default = 1.0.
#' @param cyc Positive integer that sets the number of cycles for the algorithm.
#' Usually after less than 15 cycles convergence is reached. Default = 15.
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
#' initiate the use of our fast segmentation algorithm.
#' @param ... Additional parameters will be passed to the "DNAcopy"
#' or the standart segmentation algorithm.
#' @examples 
#' data(cn.mops)
#' cn.mops(XRanges)
#'
#' @useDynLib cn.mops
#' @return An instance of "CNVDetectionResult".
#' @author Guenter Klambauer \email{klambauer@@bioinf.jku.at}
#' @export


cn.mops <- function(input,I = c(0.025,0.5,1,1.5,2,2.5,3,3.5,4),
		classes=c("CN0","CN1","CN2","CN3","CN4","CN5","CN6","CN7","CN8"),
		priorImpact = 0.1,cyc = 20,parallel=0,
		normType="poisson",normQu=0.25,norm=TRUE,
		upperThreshold=0.5,lowerThreshold=-0.9,
		minWidth=3,segAlgorithm="fast",...){
	
	version <- "0.99"
	
	
	if(class(input)=="GRanges"){
		inputType <- "GRanges"
		#X <- (do.call("cbind",(values(input)@unlistData@listData)))
		X <- do.call("cbind",values(input)@listData)
		X <- matrix(as.numeric(X),nrow=nrow(X))
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
	if (!("CN2" %in% classes)){stop("One element of classes must be CN2 .\n")}
	
	m <- nrow(X)
	N <- ncol(X)
	n <- length(I)
	
	if (is.null(colnames(X))){
		colnames(X) <- as.character(1:N)
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
		# Normalization only for visualization purposes
		X.norm <- normalizeChromosomes(X,chr=chr,normType=normType,qu=normQu)
	}
	message("Starting local modeling, please be patient...  ")
	
	if (parallel > 0){
		library(snow)
		cl <- makeCluster(as.integer(parallel),type="SOCK")
		clusterEvalQ(cl,".cn.mopsC")
	} 
	
	res <- list()
	for (chrom in unique(chr)){
		message(paste("Reference sequence: ",chrom))
		chrIdx <- which(chr==chrom)
		
		if(m > 1 & !norm){
			cov <- colSums(X[chrIdx, ])
			if (median(cov) > 0){
				cov <- cov/median(cov)
			} else {
				stop(paste("Median of total reads is zero.",
								"Too many samples with zero reads."))
			}
		} else {
			#X.norm <- X
			cov <- rep(1,N)
		}
		
		#cat("Coverage: ",cov ," Norm: ",norm)
		
		cov[which(cov < 1e-15)] <- 1e-15
		
		if (norm & m > 1){
			if (parallel==0){
				resChr <-apply(X.norm[chrIdx, ,drop=FALSE],1,.cn.mopsC, I=I,
						classes=classes,
						cov=cov,priorimpact=priorImpact,cyc=cyc)
			} else {
				resChr <- parApply(cl,X.norm[chrIdx, ,drop=FALSE],1,.cn.mopsC, 
						I=I,classes=classes,cov=cov,priorimpact=priorImpact,
						cyc=cyc)				
			}
		} else {
			if (parallel==0){
				resChr <-apply(X[chrIdx, ,drop=FALSE],1,.cn.mopsC, I=I,
						classes=classes,
						cov=cov,priorimpact=priorImpact,cyc=cyc)
			} else {
				resChr <- parApply(cl,X[chrIdx, ,drop=FALSE],1,.cn.mopsC, 
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
	A <- t(sapply(res,.subset2,2))
	rownames(A) <- rownames(X)
	colnames(A) <- classes
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
			colnames(segDf) <- c("chr","from","to","sample")
			
			segCN <- apply(segDf,1,function(x){
						sIdx <- as.integer(x["from"]):as.integer(x["to"])
						if (length(sIdx)>=3){
							sIdx2 <- sIdx[-c(1,length(sIdx))]
						} else{
							sIdx2 <- sIdx
						}
						if (norm){
							CN <- .cn.mopsC(colSums(X.norm[sIdx2, ,drop=FALSE]),
									I=I,
									classes=classes,
									cov=cov,priorimpact=priorImpact,
									cyc=cyc)$expectedCN
						}	else {
							CN <- .cn.mopsC(colSums(X.norm[sIdx2, ,drop=FALSE]), 
									I=I,
									classes=classes,
									cov=rep(1,N),priorimpact=priorImpact,
									cyc=cyc)$expectedCN
						}
						
						#Median
						#segValue <- quantile(sINI[sIdx,x["sample"]],probs=0.5,
						#		na.rm=TRUE)
						
						#Mean
						segValue <- mean(sINI[sIdx,x["sample"]],na.rm=TRUE)
						
						return(data.frame(CN[which(x["sample"]==colnames(X))],
										segValue,length(sIdx),
										stringsAsFactors=FALSE))
						
					})
			
			value <- sapply(segCN,.subset2,2)
			
			callsS <- matrix(rep(value,sapply(segCN,.subset2,3)),ncol=N)
			colnames(callsS) <- colnames(X)
			
			segDf <- data.frame(segDf,"CN"=sapply(segCN,.subset2,1),
					"value"=value,
					stringsAsFactors=FALSE)
			
			
			segDf <- segDf[which(segDf$value >= upperThreshold
									| segDf$value <= lowerThreshold), ]
			segDf <- segDf[which((segDf$to-segDf$from+1) >= minWidth), ]
			
			
		} else {
			resSegmList <- list()
			segDf <- data.frame(stringsAsFactors=FALSE)
			if (!exists("segMedianT")){
				segMedianT <- 0.05
			} else {
				stop(paste("Use \"upper-\" and \"lowerThreshold\" parameters",
								"instead of segMedianT."))
			}
			
			#if (!exists("alpha")){alpha <- 0.1}
			
			for (chrom in unique(chr)){
				chrIdx <- which(chr==chrom)
				
				if (parallel==0){
					resSegmList[[chrom]] <- apply(sINI[chrIdx, ],2,segment,
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
				segDf <- rbind(segDf,segDfTmp)
			}
			#browser()
			
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
						if (norm){
							CN <- .cn.mopsC(colSums(X.norm[sIdx, ,drop=FALSE]),
									I=I,
									classes=classes,
									cov=cov,priorimpact=priorImpact,
									cyc=cyc)$expectedCN
						}	else {
							CN <- .cn.mopsC(colSums(X.norm[sIdx, ,drop=FALSE]), 
									I=I,
									classes=classes,
									cov=rep(1,N),priorimpact=priorImpact,
									cyc=cyc)$expectedCN
						}
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
			#browser()
			
			colnames(segDf) <- c("from","to","mean","median","sample",
					"chr","CN")
			
			#value <- segDf$value
			
			#median
			#callsS <- matrix(rep(segDf$value,segDf$to-segDf$from+1),ncol=N)
			#segDf <- segDf[which(segDf$value >= upperThreshold
			#						| segDf$value <= lowerThreshold), ]
			
			#mean for AUC
			callsS <- matrix(rep(segDf$mean,segDf$to-segDf$from+1),ncol=N)
			colnames(callsS) <- colnames(X)
			
			#median for the table
			segDfSubset <- segDf[which(segDf$median >= upperThreshold
									| segDf$median <= lowerThreshold), ]
			
			
			#segDf$value <- segDf$value2
			
			
			segDfSubset <- 
					segDfSubset[which((segDfSubset$to-segDfSubset$from+1)
											>= minWidth), ]
			
		}
		
		#browser()
		
		
		
		if (nrow(segDfSubset)>0){
			#browser()
			# Assembly of result object
			r <- new("CNVDetectionResult")
			
			sampleNames <- segDfSubset$sample
			
			if (inputType=="GRanges"){
				ir <- IRanges(start(input)[segDfSubset$from],
						end(input)[segDfSubset$to])
			} else if (inputType=="DataMatrix"){
				ir <- IRanges(start=segDfSubset$from,end=segDfSubset$to)
			}
			rd <- GRanges(seqnames=segDfSubset$chr,ir,"sampleName"=sampleNames,
					"median"=segDfSubset$median,"mean"=segDfSubset$mean,
					"CN"=segDfSubset$CN)
			
			cnvr <- reduce(GRanges(seqnames=segDfSubset$chr,ir))
			
			r@normalizedData    <- GRanges(seqnames=chr,irAllRegions,
					normalizedData=X.norm)
			r@localAssessments  <- GRanges(seqnames=chr,irAllRegions,
					localAssessments=sINI)
			#write.table(sINI,file="sINIafter.txt")
			
			r@individualCall   	<- GRanges(seqnames=chr,irAllRegions,
					individualCall=callsS)
			r@iniCall        	<- GRanges(seqnames=chr,irAllRegions,
					iniCall=INI)
			r@cnvs				<- rd
			r@cnvr				<- cnvr
			
			
			if (inputType=="GRanges"){
				r@segmentation 			<- GRanges(seqnames=segDf$chr,
						IRanges(start(input)[segDf$from],end(input)[segDf$to]),
						"sampleName"=segDf$sample,"median"=segDf$median,
						"mean"=segDf$mean,"CN"=segDf$CN)
			} else if (inputType=="DataMatrix"){
				r@segmentation 			<- GRanges(seqnames=segDf$chr,
						IRanges(segDf$from,segDf$to),"sampleName"=segDf$sample,
						"median"=segDf$median,
						"mean"=segDf$mean,"CN"=segDf$CN)
			}
			
			r@posteriorProbs 	<- post
			r@params			<- params
			r@integerCopyNumber	<- GRanges(seqnames=chr,irAllRegions,
					integerCopyNumber=CN)
			return(r)	
		} else {
			message(paste("No CNVs detected. Try changing \"normalization\",", 
							"\"priorimpact\" or \"thresholds\"."))
			# Assembly of result object
			r <- new("CNVDetectionResult")
			if (inputType=="GRanges"){
				r@segmentation 			<- GRanges(seqnames=segDf$chr,
						IRanges(start(input)[segDf$from],end(input)[segDf$to]),
						"sampleName"=segDf$sample,"median"=segDf$median,
						"mean"=segDf$mean,"CN"=segDf$CN)
			} else if (inputType=="DataMatrix"){
				r@segmentation 			<- GRanges(seqnames=segDf$chr,
						IRanges(segDf$from,segDf$to),
						"sampleName"=segDf$sample,"median"=segDf$median,
						"mean"=segDf$mean,"CN"=segDf$CN)
			}
			
			r@normalizedData    <- GRanges(seqnames=chr,irAllRegions,
					normalizedData=X.norm)
			r@localAssessments  <- GRanges(seqnames=chr,irAllRegions,
					localAssessments=sINI)
			
			r@individualCall   	<- GRanges(seqnames=chr,irAllRegions,
					individualCall=callsS)
			r@iniCall        	<- GRanges(seqnames=chr,irAllRegions,
					iniCall=INI)
			r@posteriorProbs 	<- post
			r@params			<- params
			r@integerCopyNumber	<- GRanges(seqnames=chr,irAllRegions,
					integerCopyNumber=CN)
			
			return(r)	
		}
		
		
	} else {
		message(paste("Only one genomic segment considered, therefore no ",
						"segmentation."))	
		# Assembly of result object
		r <- new("CNVDetectionResult")	#
		
		r@normalizedData    <- GRanges(seqnames=chr,irAllRegions,
				normalizedData=X.norm)
		r@localAssessments  <- GRanges(seqnames=chr,irAllRegions,
				localAssessments=sINI)
		r@individualCall   	<- GRanges(seqnames=chr,irAllRegions,
				individualCall=sINI)
		
		r@params			<- params
		
		r@integerCopyNumber	<- GRanges(seqnames=chr,irAllRegions,
				integerCopyNumber=CN)
		
		return(r)	
		
	}
}
