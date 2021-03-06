% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/partitionLongreads.R
\name{HapPart}
\alias{HapPart}
\alias{mcoef.HapPart}
\alias{PWM.HapPart}
\alias{OC.HapPart}
\alias{SCR.HapPart}
\alias{K.HapPart}
\alias{SQS.HapPart}
\alias{CLS.HapPart}
\alias{SNP.HapPart}
\alias{PRT.HapPart}
\alias{HTR.HapPart}
\title{HapPart (Haplotype Partition) stores the haplotype information from
longread clustering.}
\usage{
HapPart(readNames, snpPos)

\method{mcoef}{HapPart}(x)

\method{PWM}{HapPart}(x)

\method{OC}{HapPart}(x)

\method{SCR}{HapPart}(x)

\method{K}{HapPart}(x)

\method{SQS}{HapPart}(x)

\method{CLS}{HapPart}(x)

\method{SNP}{HapPart}(x)

\method{PRT}{HapPart}(x)

\method{HTR}{HapPart}(x)
}
\arguments{
\item{readNames}{The names of reads in each cluster.}

\item{snpPos}{SNP positions used for clustering.}

\item{x}{A \code{\link{HapPart}} object.}
}
\value{
HapPart object.
}
\description{
HapPart (Haplotype Partition) stores the haplotype information from
longread clustering.
}
\section{Methods (by generic)}{
\itemize{
\item \code{mcoef}: The membership coefficient in the assigned clade.

\item \code{PWM}: Get the Position Weight Matrix

\item \code{OC}: The original clusters from \code{\link[stats]{hclust}}.

\item \code{SCR}: The membership scores for each read in each cluster.

\item \code{K}: The total number of polymorphic positions used.

\item \code{SQS}: Get the consensus sequences of the haplotypes.

\item \code{CLS}: Get the classification of the haplotype partitioning.

\item \code{SNP}: Get the SNP positions.

\item \code{PRT}: The haplotype assignment as a character vector <A>, <B>, ... for
each read.

\item \code{HTR}: The cluster tree produced by \code{\link[stats]{hclust}}.
}}

\section{Slots}{

\describe{
\item{\code{readNames}}{The names of reads in each cluster.}

\item{\code{snpPos}}{SNP positions used for clustering.}

\item{\code{mcoef}}{membership coefficient. Score for each read for belonging
to its cluster. Based on the sum of probabilities of a PWM of the haplotype
sequences.}

\item{\code{tree}}{The resulting tree from \code{\link[stats]{hclust}}.}

\item{\code{scores}}{The same scores as in mcoef, but for all clusters.}

\item{\code{mats}}{PWM matrices of the haplotype sequences.}

\item{\code{classification}}{... TODO}
}}

