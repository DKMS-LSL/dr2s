% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/all-generics.R
\name{mapFinal}
\alias{mapFinal}
\title{Final mapping (long reads and short reads).}
\usage{
mapFinal(x, opts = list(), ...)
}
\arguments{
\item{x}{A \code{\link[=DR2S_]{DR2S}} object.}

\item{opts}{List with options passed to the mapper.}

\item{...}{Additional parameters passed to \code{\link[Rsamtools]{PileupParam}}.}
}
\value{
A \code{\link[=DR2S_]{DR2S}} object.
}
\description{
Partitioned long reads and partitioned or unpartitioned short reads are
mapped against the consensus sequences inferred at \code{mapIter}. Indels
are \strong{included} in pileup.
}
\examples{
\dontrun{
x <- DR2Smap(
  sample = "ID12300527",
  locus = "DPB1*04:01:01:01",
  datadir = "/path/to/data",
  outdir = "/path/to/output"
  ) \%>\%
  mapInit() \%>\%
  partitionLongreads() \%>\%
  mapIter() \%>\%
  partitionShortreads() \%>\%
  mapFinal() \%>\%
  report(blockWidth = 60)
}
}
\seealso{
Other DR2S mapper functions: 
\code{\link{cache}()},
\code{\link{clear}()},
\code{\link{createDR2SConf}()},
\code{\link{mapInit}()},
\code{\link{mapIter}()},
\code{\link{partitionLongreads}()},
\code{\link{remapAlignment}()},
\code{\link{reportCheckedConsensus}()},
\code{\link{report}()}
}
\concept{DR2S mapper functions}
