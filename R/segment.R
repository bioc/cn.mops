# Copyright (C) 2011 Klambauer Guenter 
# <klambauer@bioinf.jku.at>

#' Performs a fast segmentation algorithm based on the cyber t test
#' and the t statistics.
#' 
#' @param x Values to be segmented.
#' @param alpha Real value between 0 and 1 is interpreted as the percentage of
#' total points that are considered as initial breakpoints. An integer greater 
#' than 1 is interpreted as number of initial breakpoints.
#' @param segMedianT Vector of length 2. Thresholds on the segment's median. 
#' Segments' medians above the first element are considered as gains and below
#' the second value as losses.
#' @param minSeg Minimum length of segments.
#' @param eps Real value greater or equal zero. A breakpoint is only possible 
#' between to consecutive values of x that have a distance of at least "eps".
#' @param delta Positive integer. A parameter to make the segmentation more 
#' efficient. If the statistics of a breakpoint lowers while extending the 
#' window, the algorithm extends the windows by "delta" more points until it 
#' stops.
#' @param maxInt The maximum length of a segment left of the breakpoint and
#' right of the breakpoint that is considered. 
#' @param squashing An experimental parameter that squashes the values "x" 
#' before segmentation. Should be left to zero, which means that squashing is 
#' not performed.
#' @param cyberWeight The "nu" parameter of the cyber t-test. 
#' @param segPlot Logical indicating whether the result of the segmentation a
#' algorithm should be plotted.
#' @param ... additional parameters passed to the plotting function.
#' @examples 
#' x <- rnorm(n=500,sd=0.5)
#' x[150:200] <- rnorm(n=51,mean=3,sd=0.5)
#' segment(x)
#' @author Guenter Klambauer \email{klambauer@@bioinf.jku.at}
#' @return A data frame containing the segments.
#' @export
#' @useDynLib cn.mops


segment <- function(x, alpha=.05, segMedianT, minSeg=3, 
		eps=0, delta=5, maxInt=40, squashing=0, cyberWeight=5,
		segPlot=TRUE, ...){
	
	if (missing("segMedianT")) {
		segMedianT <- c()
		segMedianT[1] <- mean(x, na.rm=TRUE)+2*sd(x, na.rm=TRUE)
		segMedianT[2] <- mean(x, na.rm=TRUE)-2*sd(x, na.rm=TRUE)
		
	} else {
		if (length(segMedianT)==1){
			segMedianT <- c(abs(segMedianT), -abs(segMedianT))
		}
	}
	if (any(is.na(x))){
		message("NA values detected. Replacing with median.")
		x[is.na(x)] <- median(x, na.rm=TRUE)
	}
	
	if (missing("eps")){
		#eps <- min(abs(quantile(x, 0.05)), abs(quantile(x, 0.95)))
		eps <- quantile(abs(diff(x)), probs=0.75)
	}
	
	res <-  .Call("segment", x, as.double(eps), as.integer(delta),
			as.integer(maxInt), as.integer(minSeg),
			as.integer(squashing), as.double(cyberWeight))
	
	#message("Finished C function.")
	
	if (alpha > 1){
		alpha <- as.integer(alpha)
		#message(paste("Number of initial breakpoints: ",alpha))
		brkptsInit <- sort(order(res$stat,decreasing=TRUE)[1:alpha])
		
	} else if (alpha < 1 & alpha > 0){
		#message(paste("Number of initial breakpoints: ", 
		#				as.integer(alpha*length(x)) ))
		pValT <- quantile(res$stat, probs=1-alpha)
		brkptsInit <- which(res$stat > pValT)
	} else{
		stop(paste("Alpha must be either between 0 and 1 or an integer",
						"greater than 1."))
	}
	
	brkptsInit <- unique(c(0, brkptsInit, length(x)))
	nbrOfBrkpts <- length(brkptsInit)-1
	
	#plot(x, pch=15, cex=0.5, ...)
	med <- vector(length=nbrOfBrkpts)
	m <- vector(length=nbrOfBrkpts)
	start <- vector(length=nbrOfBrkpts)
	end <- vector(length=nbrOfBrkpts)
	for (i in 1:nbrOfBrkpts){
		y <- x[(brkptsInit[i]+1):brkptsInit[i+1]]
		m[i] <- mean(y,na.rm=TRUE)
		med[i] <- quantile(y,probs=0.5,na.rm=TRUE)
		start[i] <- brkptsInit[i]+1
		end[i] <- brkptsInit[i+1]
		#   lines(c(start[i], end[i]), c(med[i], med[i]), lwd=5, col="green")
	}
	
	df <- data.frame("start"=start, "end"=end, "mean"=m, "median"=med)
	
	if (all(segMedianT==0)) {
		
		#message("No merging of segments.")
		ir <- IRanges(df$start, df$end)
		ir <- ir[which(ir@width>=minSeg)]
		
		
		irAll <- IRanges(1, length(x))
		segsFinal <- as.data.frame(sort(c(ir, setdiff(irAll, ir))))
		
		if (segPlot) plot(x, pch=15, cex=0.5, ...)
		nbrOfSegs <- nrow(segsFinal)
		med <- vector(length=nbrOfSegs)
		m <- vector(length=nbrOfSegs)
		start <- vector(length=nbrOfSegs)
		end <- vector(length=nbrOfSegs)
		for (i in 1:nbrOfSegs) {
			y <- x[segsFinal$start[i]:segsFinal$end[i]]
			m[i] <- mean(y)
			med[i] <- median(y)
			
			if (segPlot) lines(c(segsFinal$start[i], segsFinal$end[i]), 
						c(med[i], med[i]), lwd=5, col="green")
		}
		
		df2 <- data.frame("start"=segsFinal$start, "end"=segsFinal$end, 
				"mean"=m, "median"=med)
		
		
		return(df2)
		
		
	} else {
		dfAmp <- df[which(df$median > segMedianT[1]), ]
		irAmp <- IRanges(dfAmp$start, dfAmp$end)
		irAmp <- reduce(irAmp)
		
		dfLoss <- df[which(df$median < segMedianT[2]), ]
		irLoss <- IRanges(dfLoss$start, dfLoss$end)
		irLoss <- reduce(irLoss)
		
		ir <- sort(c(irAmp, irLoss))
		ir <- ir[which(ir@width>=minSeg)]
		
		rm(irAmp, irLoss, dfAmp, dfLoss)    
		
		irAll <- IRanges(1, length(x))
		segsFinal <- as.data.frame(sort(
						c(ir, setdiff(irAll, ir))))
		
		if (segPlot) plot(x, pch=15, cex=0.5, ...)
		nbrOfSegs <- nrow(segsFinal)
		med <- vector(length=nbrOfSegs)
		m <- vector(length=nbrOfSegs)
		start <- vector(length=nbrOfSegs)
		end <- vector(length=nbrOfSegs)
		for (i in 1:nbrOfSegs) {
			y <- x[segsFinal$start[i]:segsFinal$end[i]]
			m[i] <- mean(y)
			med[i] <- median(y)
			
			if (segPlot) lines(c(segsFinal$start[i], segsFinal$end[i]), 
						c(med[i], med[i]), lwd=5, col="green")
		}
		
		df2 <- data.frame("start"=segsFinal$start, "end"=segsFinal$end, 
				"mean"=m, "median"=med)
		
		
		return(df2)
	}
}


