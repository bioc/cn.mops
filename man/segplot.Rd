\name{segplot}
\alias{segplot}
\title{Plots the log normalized read counts and the detected segments for
one sample on one reference sequence.}
\arguments{
  \item{r}{An instance of "CNVDetectionResult"}

  \item{seqname}{The name of the reference sequence or
  chromosome to be plotted.}

  \item{sampleIdx}{The index of the sample as it appears in
  the read count matrix.}
}
\value{
  Generates a segmentation plot.
}
\description{
  Plots the log normalized read counts and the detected
  segments for one sample on one reference sequence.
}
\examples{
data(cn.mops)
r <- cn.mops(X[1:200, ])
segplot(r,sampleIdx=1)
}
\author{
  Guenter Klambauer \email{klambauer@bioinf.jku.at}
}

