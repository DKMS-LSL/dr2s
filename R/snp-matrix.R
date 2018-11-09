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
#bamfile = bampath(self$mapInit)
#refseq = self$getRefSeq()
#polymorphicPositions = ppos
SNPmatrix <- function(bamfile, polymorphicPositions) {
  assert_that(
    requireNamespace("readr", quietly = TRUE),
    file.exists(bamfile)
  )
  if (is(polymorphicPositions, "tbl_df") &&
      all(colnames(polymorphicPositions) %in% c("position", "a1", "f1", "a2", "f2"))) {
    polymorphicPositions <- polymorphicPositions %>%
      dplyr::filter(a1 %in% VALID_DNA(), a2 %in% VALID_DNA()) %>%
      dplyr::select(position) %>%
      unlist(use.names = FALSE) %>%
      as.integer()
  }
  if (is.character(polymorphicPositions)) {
    polymorphicPositions <- as.integer(polymorphicPositions)
  }
  assert_that(is.integer(polymorphicPositions))
  msa <- .msaFromBam(Rsamtools::BamFile(bamfile))
  mat <- vapply(polymorphicPositions, function(x, msa)
    as.matrix(Biostrings::subseq(msa, start = x, width = 1L)),
    msa = msa, FUN.VALUE = character(length(msa)))
  rownames(mat) <- names(msa)
  colnames(mat) <- polymorphicPositions
  mat
}
