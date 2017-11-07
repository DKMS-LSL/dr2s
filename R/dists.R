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
#'   \item{read_ids}{<character>; Read IDs}
#' }
#'
#' @return A \code{consmat} object.
#' @export
#' @examples
#' ###
hammingDist <- function(x){
  stopifnot(is(x, "DNAStringSet"))
  x_tmp <- as.matrix(x)
  x_mat<- plyr::revalue(x_tmp, c("G" = 1, "A" = 2, "T" = 3, "C" = 4, "-" = 5, "+" = 6), warn_missing = FALSE)
  x_mat<- sapply(x_mat, as.numeric)
  dim(x_mat) <- dim(x_tmp)
  rm(x_tmp)
  dist <- cpp_hamming(x_mat)
  colnames(dist) <- names(x)
  rownames(dist) <- names(x)
  dist
}

# Todo Export and write man entry
PSDM <- function(x, consmat){
  stopifnot(is(x, "DNAStringSet"))

  ## Create seq matrix as input for cpp_PSDM
  x_tmp <- as.matrix(x)
  x_mat<- plyr::revalue(x_tmp, c("G" = 1, "A" = 2, "T" = 3, "C" = 4, "-" = 5, "+" = 6), warn_missing = FALSE)
  x_mat<- sapply(x_mat, as.numeric)
  dim(x_mat) <- dim(x_tmp)
  rm(x_tmp)

  ## Get Position Specific Distance Matrix
  dist <- cpp_PSDM(consmat, x_mat)
  colnames(dist) <- names(x)
  rownames(dist) <- names(x)
  dist
}
