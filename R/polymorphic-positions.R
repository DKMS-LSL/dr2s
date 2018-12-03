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
  tibble::tibble(
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

.selectCorrelatedPolymorphicPositions <- function(mat,
                                                  method = "cramer.V",
                                                  dunnCutoff = 15,
                                                  minimumMeanAssociation = 0.1,
                                                  resample = NULL,
                                                  ...) {
  assert_that(is(mat, "matrix"))
  method <- match.arg(method, c("cramer.V", "spearman"))
  indent <- list(...)$indent %||% indentation()
  fmat <- factor(mat, levels = VALID_DNA("indel"), labels = VALID_DNA("indel"))
  dim(fmat) <- dim(mat)
  rownames(fmat) <- rownames(mat)
  colnames(fmat) <- colnames(mat)

  ## this is for purely exploratory purposess
  if (!is.null(resample)) {
    col <- sample(colnames(fmat), resample)
    for (i in seq_along(col)) {
      fmat[, col][, i] <- sample(fmat[, col][, i])
    }
  }

  if (method == "cramer.V") {
    ## Cramér's V association of nxm contingency tables
    cmat <- matrix(NA_real_, nrow = NCOL(fmat), ncol = NCOL(fmat))
    colnames(cmat) <- colnames(fmat)
    rownames(cmat) <- colnames(fmat)
    cCmb <- utils::combn(colnames(cmat), 2, simplify = FALSE)
    cTab <- purrr::map(cCmb, ~table(fmat[, .[1]], fmat[, .[2]], dnn = .))
    cV   <- purrr::map_dbl(cTab, cramerV)
    cmat[lower.tri(cmat, diag = FALSE)] <- cV
    cmat <- t(cmat)
    cmat[lower.tri(cmat, diag = FALSE)] <- cV
    y.lab <- "Mean Cramér's V"
    legend.label <- "V"
  } else if (method == "spearman") {
    nmat <- as.numeric(fmat)
    dim(nmat) <- dim(fmat)
    rownames(nmat) <- rownames(fmat)
    colnames(nmat) <- colnames(fmat)
    cmat <- abs(stats::cor(nmat, method = "spearman"))
    diag(cmat) <- NA_real_
    y.lab <- "Mean Spearman's rho"
    legend.label <- "rho"
  }

  ## calculate mean association of each PP too each other PP
  crm <- .rowMeans(cmat, m = NCOL(cmat), n = NROW(cmat), na.rm = TRUE)
  names(crm) <- colnames(cmat)
  ## try to split PP in a low association and a high association cluster.
  ## evaluate using a Dunn-like index, the ratio of the largest distance
  ## between the sorted average association (i.e., the smallest distance
  ## between the two putative clusters) to the mean average association
  ## (approximately the mean intra-cluster distance).
  scrm <- sort(crm)  # ranked mean association
  dcrm <- diff(scrm) # ranked difference in mean association
  c1 <- scrm[1:which.max(dcrm)] ## putative low association cluster
  c2 <- scrm[-(1:which.max(dcrm))] ## putative high association cluster
  dunnIndex <- max(dcrm)/mean(dcrm)
  if (dunnIndex > dunnCutoff) {
    nm <- names(c2)
    flog.info("%sInferring cluster of low association polymorphisms at Dunn index <%0.1f>",
              indent(), dunnIndex, name = "info")
    flog.info("%sRemoving %s low-association polymorphic positions ", indent(), length(c1), name = "info")
  } else {
    nm <- names(crm)
    flog.info("%sInferring no cluster of low association polymorphisms at Dunn index <%0.1f>",
              indent(), dunnIndex, name = "info")
  }
  ## use minimumMeanAssociation as a hard cutoff
  if (any(idx <- crm[nm] >= minimumMeanAssociation)) {
    flog.info("%sRemoving %s additional polymorphic positions below minimum mean association threshold <%0.2f>", indent(),
              sum(!idx), minimumMeanAssociation, name = "info")
    nm <- names(which(idx))
  }
  flog.info("%sRetaining %s high-association polymorphic positions", indent(),
            length(nm), name = "info")
  ## Correlogram
  cl <- stats::hclust(stats::as.dist(1 - cmat), method = "ward.D")
  ocmat <- cmat[cl$order, cl$order]
  ocmat[upper.tri(ocmat)] <- NA_real_
  cDf <- tibble::as_tibble(ocmat)
  ppos <- colnames(cDf)
  cDf$pos <- factor(ppos, levels = ppos, labels = ppos, ordered = TRUE)
  cDfLong <- tidyr::gather(cDf, pos1, Assoc, -pos, na.rm = TRUE) %>%
    dplyr::mutate(pos1 = factor(pos1, levels = ppos, labels = ppos, ordered = TRUE))
  p1 <- ggplot(cDfLong, aes(x = reorder(pos, dplyr::desc(pos)), y = pos1)) +
    geom_tile(aes(fill = Assoc, height = 0.95, width = 0.95)) + coord_flip() +
    scale_fill_gradient2(mid = "#e1e7fa", high = "#13209d", midpoint = 0, limit = c(0, 1)) +
    labs(x = "Polymorphic positions", y = "Polymorphic positions", fill = legend.label) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 90))
  rDf <- tibble::tibble(
    Position = factor(names(crm), levels = names(crm), labels = names(crm), ordered = TRUE),
    meanAssoc = unname(crm),
    picked = ifelse(names(crm) %in% nm, "yes", "no")
  )
  p2 <- ggplot(rDf, aes(x = Position, y = meanAssoc, colour =  picked)) +
    geom_point() +
    expand_limits(y = 0) +
    labs(x = "Polymorphic position", y = y.lab) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
          legend.position = "bottom")

  structure(nm, snp.corr.mat = cmat, dunn.index = dunnIndex,
            snp.correlogram = p1, snp.association = p2)
}

cramerV <- function(x) {
  ## x is an nxn contingency table of observations
  ## calc chi2 stats first
  n  <- sum(x)
  nr <- NROW(x)
  nc <- NCOL(x)
  sr <- rowSums(x)
  sc <- colSums(x)
  expected <- outer(sr, sc, "*")/n
  chi2 <- sum(abs(x - expected)^2/expected, na.rm = TRUE)
  ## Cramér's V
  k <- NCOL(x)
  V <- sqrt(chi2/(n * (k - 1)))
  V
}

