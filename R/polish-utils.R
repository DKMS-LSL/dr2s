#' Print diagnostic summary of a variant
#'
#' @param x A \code{\link[=DR2S_]{DR2S}} mapper object
#' @param i Variant index for haplotype A
#' @param j Variant index for haplotype B
#'
#' @export
inspect_variant <- function(x, i, j = i) {
  cat("Haplotype A: ", sep = "")
  print(x$consensus$A$variants[[i]])
  cat("\n\nHaplotype B: ", sep = "")
  print(x$consensus$B$variants[[j]])
  invisible(NULL)
}
