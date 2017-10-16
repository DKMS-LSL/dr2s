#' Construct a SNP matrix.
#'
#' @param bamfile A BAM file
#' @param polymorphic_positions SNP positions.
#' @param refseq Reference sequence as DNAString
#'
#' @return A \code{[read x position]} matrix.
#' @export
#' @examples
#' ##

# debug
#bamfile = self$mapInit$bamfile
#refseq = self$getRefSeq()
#polymorphic_positions = ppos
SNPmatrix <- function(bamfile,
                      refseq,
                      polymorphic_positions) {

  assertthat::assert_that(requireNamespace("readr", quietly = TRUE))
  if (is(polymorphic_positions, "tbl_df") &&
      all(colnames(polymorphic_positions) %in% c("position", "a1", "f1", "a2", "f2"))) {
    DNA <- c("A", "C", "G", "T")
    polymorphic_positions <- polymorphic_positions %>%
      dplyr::filter_(~a1 %in% DNA, ~a2 %in% DNA) %>%
      dplyr::select_(~position) %>%
      unlist(use.names = FALSE) %>%
      as.integer()
  }

  if (is.character(polymorphic_positions)) {
    polymorphic_positions <- as.integer(polymorphic_positions)
  }
  assertthat::assert_that(is.numeric(polymorphic_positions))
  mat <- msa_from_bam(bamfile, refseq)
  mat <- as.matrix(mat)
  mat <- as.matrix(mat[,polymorphic_positions])
  colnames(mat) <- polymorphic_positions
  mat
}
