#' Construct a SNP matrix.
#'
#' @param bamfile A BAM file
#' @param polymorphicPositions SNP positions.
#' @param refseq Reference sequence as DNAString
#'
#' @return A \code{[read x position]} matrix.
#' @export
#' @examples
#' ##

# debug
#bamfile = self$mapInit$bamfile
#refseq = self$getRefSeq()
#polymorphicPositions = ppos
SNPmatrix <- function(bamfile,
                      refseq,
                      polymorphicPositions) {
  assertthat::assert_that(requireNamespace("readr", quietly = TRUE))
  if (is(polymorphicPositions, "tbl_df") &&
      all(colnames(polymorphicPositions) %in% c("position", "a1", "f1", "a2", "f2"))) {
    polymorphicPositions <- polymorphicPositions %>%
      dplyr::filter_(~a1 %in% VALID_DNA(), ~a2 %in% VALID_DNA()) %>%
      dplyr::select_(~position) %>%
      unlist(use.names = FALSE) %>%
      as.integer()
  }

  if (is.character(polymorphicPositions)) {
    polymorphicPositions <- as.integer(polymorphicPositions)
  }
  assertthat::assert_that(is.numeric(polymorphicPositions))
  msa <- .msaFromBam(bamfile, refseq)
  mat <- sapply(polymorphicPositions, function(x) as.matrix(Biostrings::subseq(msa, start = x, width = 1)))
  rownames(mat) <- names(msa)
  colnames(mat) <- polymorphicPositions
  mat
}
