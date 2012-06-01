\name{gr}
\alias{gr}
\title{These generic function returns the genomic ranges
 of a CNV detection method stored in an instance of
\code{\link{CNVDetectionResult-class}}.}
\arguments{
  \item{object}{An instance of "CNVDetectionResult".}
}
\value{
  \code{normalizedData} returns a "GRanges" object
  containing the normalized data.
}
\description{
  These generic function returns the genomic ranges of a
  CNV detection method stored in an instance of
  \code{\link{CNVDetectionResult-class}}.
}
\examples{
data(cn.mops)
r <- cn.mops(X[1:100,1:5])
gr(r)
}
\author{
  Guenter Klambauer \email{klambauer@bioinf.jku.at}
}

