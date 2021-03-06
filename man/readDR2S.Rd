% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DR2S-utils.R
\name{readDR2S}
\alias{readDR2S}
\title{Read a cached DR2S object}
\usage{
readDR2S(path)
}
\arguments{
\item{path}{the path to the folder storing the DR2S object. Usually the
\code{outdir} of a previous run.}
}
\description{
Read a cached DR2S object
}
\examples{
\dontrun{
library(DR2S)
path <- "path"
dr2s <- readDR2S(path)
}
}
