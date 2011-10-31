# Copyright (C) 2011 Klambauer Guenter 
# <klambauer@bioinf.jku.at>

# S4 class definition for the result object of a CNV detection method

#' @export

setClass("CNVDetectionResult",
		representation = representation
				(
				normalizedData     	= "GRanges",
				localAssessments    = "GRanges",
				individualCall      = "GRanges",
				iniCall        		= "GRanges",
				posteriorProbs		= "array",
				cnvs				= "GRanges",
				cnvr				= "GRanges",
				segmentation		= "GRanges",
				integerCopyNumber	= "GRanges",
				params				= "list"
		),
		prototype = prototype
				(
				normalizedData     	= GRanges(),
				localAssessments     = GRanges(),
				individualCall      = GRanges(),
				iniCall        		= GRanges(),
				posteriorProbs		= array(NA,dim=c(1,1,1)),
				cnvs				= GRanges(),
				cnvr				= GRanges(),
				segmentation		= GRanges(),
				integerCopyNumber	= GRanges(),
				params				= list()
				)
)
