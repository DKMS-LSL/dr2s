% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pileup.R
\name{plotPileupBasecallFrequency}
\alias{plotPileupBasecallFrequency}
\title{Plot basecall frequency}
\usage{
plotPileupBasecallFrequency(
  x,
  threshold = 0.2,
  label = "",
  drop.indels = FALSE
)
}
\arguments{
\item{x}{A \code{pileup} object.}

\item{threshold}{At which frequency do we visualise polymorphisms.}

\item{label}{Optional plot label.}

\item{drop.indels}{Don't plot indels.}
}
\value{
A \code{ggplot} object.
}
\description{
Plot basecall frequency
}
\examples{
###
}
