# Distance functions --------------------------------------------------------

#' Construct a hamming distance matrix from an aligned DNAStringSet
#'
#' @param x DNAStringset with aligned sequences
#' @details
#' \code{consmat}: a \code{matrix} with positions row names and nucleotides as
#' column manes.
#' A \code{consmat} object includes the attributes:
#' \describe{
#'   \item{n}{<integer>; Number of reads per row.}
#'   \item{freq}{<logical>; Is it a frequency matrix or a count matrix.}
#'   \item{offset}{<integer>; Offset}
#'   \item{insertions}{<integer>; Insertions}
#' }
#'
#' @return A \code{consmat} object.
#' @export
#' @examples
#' ###
hammingDist <- function(x){
  stopifnot(is(x, "DNAStringSet"))
  xTmp <- as.matrix(x)
  dnaBaseMapping <- c("G" = 1, "A" = 2, "T" = 3, "C" = 4, "-" = 5, "+" = 6)
  xTat<- plyr::revalue(xTmp, dnaBaseMapping, warn_missing = FALSE)
  xMat<- sapply(xMat, as.numeric)
  dim(xMat) <- dim(xTmp)
  rm(xTmp)
  dist <- cpp_hamming(xMat)
  colnames(dist) <- names(x)
  rownames(dist) <- names(x)
  dist
}

# Todo Export and write man entry
PSDM <- function(x, consmat){
  stopifnot(is(x, "DNAStringSet"))

  ## Create seq matrix as input for cpp_PSDM
  xTmp <- as.matrix(x)
  dnaBaseMapping <- c("G" = 1, "A" = 2, "T" = 3, "C" = 4, "-" = 5, "+" = 6)
  xMat<- plyr::revalue(xTmp, dnaBaseMapping, warn_missing = FALSE)
  xMat<- sapply(xMat, as.numeric)
  dim(xMat) <- dim(xTmp)
  rm(xTmp)

  ## Get Position Specific Distance Matrix
  dist <- cpp_PSDM(consmat, xMat)
  colnames(dist) <- names(x)
  rownames(dist) <- names(x)
  dist
}
