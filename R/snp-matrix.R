#' Construct a SNP matrix.
#'
#' @param bamfile A BAM file
#' @param max_depth How many reads do we consider.
#' @param polymorphic_positions SNP positions.
#' @param cleanup Remove temporary files.
#'
#' @return A \code{[read x position]} matrix.
#' @export
#' @examples
#' ##
SNPmatrix <- function(bamfile,
                      max_depth,
                      polymorphic_positions,
                      cleanup = TRUE) {
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
  pm <- as.integer(polymorphic_positions) - 1 # Python is 0-indexed!
  if (missing(max_depth)) {
    max_depth <- 8000L
  }
  outfile <- tempfile(fileext = ".csv")
  if (cleanup) {
    on.exit(unlink(outfile))
  }
  if (rPython::python.call("py_get_snp_matrix", bamfile, outfile, max_depth, pm, simplify = FALSE)) {
    ct <- paste0(c("c", rep("i", length(pm))), collapse = "")
    ans <- readr::read_csv(outfile, col_types = ct, trim_ws = TRUE)
    mat <- as.matrix(ans[, -1])
    ## Use directly a matrix of characters; Better handling and no need to convert always, because we need the SNPs as
    mat <- apply(mat, 2, function(x) factor(x, levels = 0:4, labels = c("-", "A", "C", "G", "T")))
    rownames(mat) <- ans$qname
    mat
  } else{
    stop("Failed to call SNP matrix from ", basename(bamfile))
  }
}
