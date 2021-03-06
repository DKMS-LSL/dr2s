% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/config.R
\name{readDR2SConf}
\alias{readDR2SConf}
\title{Read a DR2S config file in yaml or json format}
\usage{
readDR2SConf(configFile, format = "auto")
}
\arguments{
\item{configFile}{The path to the valid DR2S config file or to a DR2S project
directrory containing a \file{config.json} or \file{config.yaml} file.}

\item{format}{Input format ("auto", "yaml" or "json")}
}
\value{
A \code{DR2Sconf} object or a list of \code{DR2Sconf} objects.
}
\description{
Read a DR2S config file in yaml or json format
}
\details{
DR2S config files can be created manually or by the
\code{\link{writeDR2SConf}} function.
}
