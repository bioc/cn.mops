# Copyright (C) 2011 Klambauer Guenter 
# <klambauer@bioinf.jku.at>


#' Displays the result object of a copy number detection method.
#'  
#' @param object An instance of a "CNVDetectionResult".
#' @examples 
#' data(cn.mops)
#' r <- cn.mops(XRanges[1:200])
#' show(r)
#' @return Displays the result object of a CNV detection method.
#' @author Guenter Klambauer \email{klambauer@@bioinf.jku.at}
#' @export
#' @noRd
#' @importFrom methods show

setMethod("show", "CNVDetectionResult",function(object){
			cat("\nIndividual CNVs: \n")
			show(object@cnvs)
			cat("\nCNV regions: \n")
			show(object@cnvr)
			
		})


