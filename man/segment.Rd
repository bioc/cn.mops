\name{segment}
\alias{segment}
\title{Fast segmentation of CNV calls.}
\usage{
  segment(x, alpha = 0.05, segMedianT = 0, minSeg = 3,
    eps = 0, delta = 20, maxInt = 40, squashing = 0,
    cyberWeight = 50)
}
\arguments{
  \item{x}{Values to be segmented.}

  \item{alpha}{Real value between 0 and 1 is interpreted as
  the percentage of total points that are considered as
  initial breakpoints. An integer greater than 1 is
  interpreted as number of initial breakpoints. Default =
  0.05.}

  \item{segMedianT}{Vector of length 2. Thresholds on the
  segment's median. Segments' medians above the first
  element are considered as gains and below the second
  value as losses.}

  \item{minSeg}{Minimum length of segments. Default = 3.}

  \item{eps}{Real value greater or equal zero. A breakpoint
  is only possible between to consecutive values of x that
  have a distance of at least "eps". Default = 0.}

  \item{delta}{Positive integer. A parameter to make the
  segmentation more efficient. If the statistics of a
  breakpoint lowers while extending the window, the
  algorithm extends the windows by "delta" more points
  until it stops. Default = 20.}

  \item{maxInt}{The maximum length of a segment left of the
  breakpoint and right of the breakpoint that is
  considered. Default = 40.}

  \item{squashing}{An experimental parameter that squashes
  the values "x" before segmentation. Should be left to
  zero, which means that squashing is not performed.
  Default = 0.}

  \item{cyberWeight}{The "nu" parameter of the cyber
  t-test. Default = 50.}
}
\value{
  A data frame containing the segments.
}
\description{
  Performs a fast segmentation algorithm based on the cyber
  t test and the t statistics.
}
\examples{
x <- rnorm(n=500,sd=0.5)
x[150:200] <- rnorm(n=51,mean=3,sd=0.5)
segment(x)
}
\author{
  Guenter Klambauer \email{klambauer@bioinf.jku.at}
}

