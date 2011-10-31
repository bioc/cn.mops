\name{normalizeChromosomes}
\alias{normalizeChromosomes}
\title{Normalize quantitative NGS data in order to make counts comparable over
samples. Scales each samples' reads such that the coverage is even for
all samples after normalization.}
\usage{
  normalizeChromosomes(X, chr, normType = "poisson", qu =
  0.25)
}
\arguments{
  \item{X}{Matrix of positive real values, where columns
  are interpreted as samples and rows as genomic regions.
  An entry is the read count of a sample in the genomic
  region.}

  \item{chr}{Character vector that has as many elements as
  "X" has rows. The vector assigns each genomic segment to
  a reference sequence (chromosome).}

  \item{normType}{Type of the normalization technique. Each
  samples' read counts are scaled such that the total
  number of reads is equal after normlization. By this
  parameter one can decide to which coverage (i.e. total
  reads) the read counts should be normalized. Possible
  choices are the minimal coverage ("min"), the mean or
  median coverage ("mean", "median") or any quantile
  ("quant"). If this parameter is set to the value "mode",
  the read counts are scaled such that each samples' most
  frequent value (the "mode") is equal after normalization.
  If the parameter is set to "poisson" the values are
  scaled such that the distribution is (rowwise) close to a
  Poisson distribution. Possible values are
  "mean","min","median","quant","poisson, and "mode".
  Default = "poisson".}

  \item{qu}{Real value between 0 and 1. Default = 0.25.}
}
\value{
  A data matrix of normalized read counts with the same
  dimensions as the input matrix X.
}
\description{
  Normalize quantitative NGS data in order to make counts
  comparable over samples. Scales each samples' reads such
  that the coverage is even for all samples after
  normalization.
}
\examples{
data(cn.mops)
X.norm <- normalizeChromosomes(X)
}
\author{
  Guenter Klambauer \email{klambauer@bioinf.jku.at}
}

