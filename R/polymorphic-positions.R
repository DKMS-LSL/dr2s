#' Calculate polymorphic positions from a consensus matrix.
#'
#' @param x A consensus matrix or pileup object.
#' @param threshold The threshold when to call a position polymorphic.
#' @param ... Additional arguments.
#'
#' @return A data.frame
#' @export
#' @examples
#' ##
polymorphicPositions <- function(x, threshold)
  UseMethod("polymorphicPositions")
#' @export
polymorphicPositions.consmat <- function(x, threshold = 0.20) {
  ## make sure to ignore insertions!
  x[, "+"] <- 0
  if (!is.freq(x)) {
    x <- consmat(x, freq = TRUE)
  }
  nuc <- colnames(x)
  pos <- rownames(x)
  i <- cpp_polymorphicPositions(x, threshold)
  rs <- cpp_top2Cols(x[i, ])
  tibble::data_frame(
    position = pos[i],
    a1 = nuc[rs$i1],
    f1 = rs$v1,
    a2 = nuc[rs$i2],
    f2 = rs$v2
  )
}
#' @export
polymorphicPositions.pileup <- function(x, threshold = NULL) {
  if (is.null(threshold)) {
    threshold <- x$threshold
  }
  polymorphicPositions(consmat(x, freq = TRUE), threshold = threshold)
}

.ambiguousPositions <- function(x, threshold)
  UseMethod(".ambiguousPositions")

.ambiguousPositions.NULL <- function(x, threshold) {
  integer(0)
}

.ambiguousPositions.consmat <- function(x, threshold) {
  f_ <- function(row) {
    sum(row > threshold) > 1L
  }
  ## make sure to ignore insertions!
  x[, "+"] <- 0
  x <- if (!is.freq(x)) {
    sweep(x, 1, n(x), `/`)
  } else x
  unname(which(apply(x, 1, f_)))
}

.ambiguousPositions.pileup <- function(x, threshold = NULL) {
  if (is.null(threshold)) {
    threshold <- x$threshold
  }
  .ambiguousPositions(consmat(x, freq = TRUE), threshold = threshold)
}

.selectCorrelatedPolymorphicPositions <- function(mat, corr.method = "pearson", clust.method = "ward.D") {
  assert_that(is(mat, "matrix"))
  fmat <- as.numeric(factor(mat, levels = VALID_DNA("indel"), labels = VALID_DNA("indel")))
  dim(fmat) <- dim(mat)
  rownames(fmat) <- rownames(mat)
  colnames(fmat) <- colnames(mat)
  cmat <- abs(stats::cor(fmat, method = corr.method))
  cl <- stats::hclust(stats::as.dist(1 - cmat), method = clust.method)
  ctr <- stats::cutree(cl, 2)
  cl1 <- names(ctr)[ctr == 1]
  cl2 <- names(ctr)[ctr == 2]
  i <- which.max(c(mean(cmat[cl1, cl1], na.rm = TRUE), mean(cmat[cl2, cl2], na.rm = TRUE)))
  nm <- names(ctr)[ctr == i]
  ## get a nice plot
  ocmat <- cmat[cl$order, cl$order]
  ## set diag and upper tri to zero
  diag(ocmat) <- NA_real_
  ocmat[upper.tri(ocmat)] <- NA_real_
  cDf <- tibble::as_tibble(ocmat)
  ppos <- colnames(cDf)
  cDf$pos <- factor(ppos, levels = ppos, labels = ppos, ordered = TRUE)
  cDfLong <- tidyr::gather(cDf, pos1, corr, -pos, na.rm = TRUE) %>%
    dplyr::mutate(pos1 = factor(pos1, levels = ppos, labels = ppos, ordered = TRUE))
  p <- ggplot(cDfLong, aes(x = reorder(pos, dplyr::desc(pos)), y = pos1, fill = corr)) +
    geom_tile(aes(height = 0.95, width = 0.95)) + coord_flip() +
    scale_fill_gradient2(mid = "#e1e7fa", high = "#13209d", midpoint = 0, limit = c(0, 1)) +
    labs(x = "Polymorphic positions", y = "Polymorphic positions") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 90))

  structure(mat[, nm], snp.corr.mat = cmat, snp.clust = cl, snp.plot = p)
}

## alternative correlogram plot
corrgram <- function(x, toFile = FALSE, width = 6, height = 6) {
  assert_that(
    is(x, "DR2S"),
    requireNamespace("corrgram", quietly = TRUE))
  outf <- file.path(x$getOutdir(), "corrgram.png")
  cmat <- attr(x$lrpartition$prt, "snp.corr.mat")
  cl <- attr(x$lrpartition$prt, "snp.clust")
  if (!is.null(cmat)) {
    mat <- cmat[cl$order, cl$order]
    if (toFile) {
      png(outf, width = width, height = height, units = "in", res = 150)
      corrgram::corrgram(mat, order = FALSE,
                         lower.panel = corrgram::panel.fill,
                         upper.panel = NULL,
                         text.panel = corrgram::panel.txt)
      dev.off()
    } else {
      corrgram::corrgram(mat, order = FALSE,
                         lower.panel = corrgram::panel.fill,
                         upper.panel = NULL,
                         text.panel = corrgram::panel.txt)
    }
  }
}

