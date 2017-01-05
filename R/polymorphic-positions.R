#' Calculate polymorphic positions from a consensus matrix.
#'
#' @param x A consensus matrix or pileup object.
#' @param theshold Threshold.
#' @param ... Additional arguments.
#'
#' @return A data.frame
#' @export
#' @examples
#' ##
polymorphic_positions <- function(x, theshold, ...) UseMethod("polymorphic_positions")

#' @export
polymorphic_positions.consmat <- function(x, threshold = 0.20) {
  if (!is.freq(x)) {
    x <- consmat(x, freq = TRUE)
  }
  nuc <- colnames(x)
  pos <- rownames(x)
  i <- cpp_polymorphic_positions(x, threshold)
  rs <- cpp_top2_cols(x[i, ])
  dplyr::data_frame(
    position = pos[i],
    a1 = nuc[rs$i1],
    f1 = rs$v1,
    a2 = nuc[rs$i2],
    f2 = rs$v2
  )
}

#' @export
polymorphic_positions.pileup <- function(x, threshold = NULL, ...) {
  if (is.null(threshold)) {
    threshold <- x$threshold
  }
  polymorphic_positions(consmat(x, freq = TRUE), threshold = threshold)
}

#' @keywords internal
#' @export
ambiguous_positions <- function(x, threshold, ...) UseMethod("ambiguous_positions")

#' @keywords internal
#' @export
ambiguous_positions.consmat <- function(x, threshold) {
  #f_ <- function(row) sum(row > threshold) > 1 || names(which.max(row)) == "-"
  f_ <- function(row) {
    sum(row > threshold) > 1L
  }
  x <- if (!is.freq(x)) {
    sweep(x, 1, n(x), `/`)
  } else x
  unname(which(apply(x, 1, f_)))
}

#' @keywords internal
#' @export
ambiguous_positions.pileup <- function(x, threshold = NULL) {
  if (is.null(threshold)) {
    threshold <- x$threshold
  }
  ambiguous_positions(consmat(x, freq = TRUE), threshold = threshold)
}

#' @keywords internal
#' @export
consensus_bases <- function(x, ...) UseMethod("consensus_bases")

#' @keywords internal
#' @export
consensus_bases.consmat <- function(x) {
  nucs <- x[, c("A", "C", "G", "T")]
  f <- sweep(nucs, 1, .rowSums(nucs, NROW(nucs), NCOL(nucs)), `/`)
  colnames(nucs)[apply(f, 1, which.max)]
}

n_polymorphic_positions <- function(cm, threshold = 0.2) {
  cm[, "-"] <- 0
  cm <- cm[rowSums(cm) != 0, ]
  length(cpp_polymorphic_positions(consmat(cm, freq = TRUE),threshold))
}
