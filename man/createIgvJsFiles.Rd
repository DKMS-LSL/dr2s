% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DR2S-utils.R
\name{createIgvJsFiles}
\alias{createIgvJsFiles}
\title{Create subsampled bam files and fasta index files for igv.js instances.}
\usage{
createIgvJsFiles(
  reference,
  bamfile,
  outdir,
  paired = FALSE,
  sampleSize = 100,
  ...
)
}
\arguments{
\item{reference}{The reference fasta file used for mapping}

\item{bamfile}{The bamfile which should be viewed in IGV.js}

\item{...}{Additional parameters passed to \code{\link{subSampleBam}}.}
}
\description{
Create subsampled bam files and fasta index files for igv.js instances.
}
