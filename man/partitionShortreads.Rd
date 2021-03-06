% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/all-generics.R
\name{partitionShortreads}
\alias{partitionShortreads}
\title{Assign short reads from mapInit to haplotypes}
\usage{
partitionShortreads(x, opts = list(), ...)
}
\arguments{
\item{x}{A \code{\link[=DR2S_]{DR2S}} object.}

\item{opts}{list with options passed to the mapper.}

\item{...}{Further arguments passed to methods.}
}
\value{
A \code{\link[=DR2S_]{DR2S}} object.
}
\description{
Assign short reads from mapInit to haplotypes
}
\examples{
\dontrun{
x <- DR2Smap(
  sample = "ID12300527",
  locus = "DPB1*04:01:01:01",
  datadir = "/path/to/data",
  outdir = "/path/to/output"
  ) \%>\%
  clear() \%>\%
  mapInit() \%>\%
  partitionLongreads() \%>\%
  mapIter() \%>\%
  partitionShortreads() \%>\%
  mapFinal() \%>\%
  report(blockWidth = 60)
}
}
\concept{DR2S partition functions}
