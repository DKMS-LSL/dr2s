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

.selectAssociatedPolymorphicPositions <- function(mat,
                                                  method.assoc = "cramer.V",
                                                  method.clust = "mclust",
                                                  expectedAbsDeviation = 0.05,
                                                  noSelect = FALSE,
                                                  resample = NULL,
                                                  ...) {
  method.assoc <- match.arg(method.assoc, c("cramer.V", "spearman", "kendall"))
  method.clust <- match.arg(method.clust, c("mclust", "dunn"))
  indent <- list(...)$indent %||% indentation()
  cmat <- .associationMatrix(mat, method = method.assoc, resample, ...)
  if (noSelect) {
    clustered.snps <- colnames(mat)
  } else {
    clustered.snps <- .clusterPolymorphicPositions(
      dist = cmat, method = method.clust, expectedAbsDeviation = expectedAbsDeviation,
      ..., indent = indent)
  }
  flog.info("%sRetaining %s high-association polymorphisms", indent(), length(clustered.snps), name = "info")
  ## create correlogram and association plots
  plts <- .correlogram(cmat, nm = clustered.snps)
  structure(clustered.snps, snp.corr.mat = cmat,
            snp.correlogram = plts$correlogram,
            snp.association = plts$association)
}

.associationMatrix <- function(mat, method = "cramer.V", resample = NULL, ...) {
  ## expect a read x snp matrix with elements {G, A, T, C, -, +}
  assert_that(is(mat, "matrix"))
  method <- match.arg(method, c("cramer.V", "spearman", "kendall"))
  indent <- list(...)$indent %||% indentation()
  fmat <- factor(mat, levels = VALID_DNA("indel"), labels = VALID_DNA("indel"))
  dim(fmat) <- dim(mat)
  rownames(fmat) <- rownames(mat)
  colnames(fmat) <- colnames(mat)

  ## this we do for purely exploratory purposess
  if (!is.null(resample)) {
    col <- sample(colnames(fmat), resample)
    for (i in seq_along(col)) {
      fmat[, col][, i] <- sample(fmat[, col][, i])
    }
  }

  if (method == "cramer.V") {
    ## Cramér's V association of nxn contingency tables
    cmat <- matrix(NA_real_, nrow = NCOL(fmat), ncol = NCOL(fmat))
    colnames(cmat) <- colnames(fmat)
    rownames(cmat) <- colnames(fmat)
    cCmb <- utils::combn(colnames(cmat), 2, simplify = FALSE)
    cTab <- purrr::map(cCmb, ~table(fmat[, .[1]], fmat[, .[2]], dnn = .))
    cV   <- purrr::map_dbl(cTab, cramerV)
    cmat[lower.tri(cmat, diag = FALSE)] <- cV
    cmat <- t(cmat)
    cmat[lower.tri(cmat, diag = FALSE)] <- cV
  }
  else if (method != "cramer.V") {
    nmat <- as.numeric(fmat)
    dim(nmat) <- dim(fmat)
    rownames(nmat) <- rownames(fmat)
    colnames(nmat) <- colnames(fmat)
    cmat <- abs(stats::cor(nmat, method = method))
  }

  diag(cmat) <- NA_real_
  attr(cmat, "method") <- method
  cmat
}

.clusterPolymorphicPositions <- function(dist, method = "mclust", expectedAbsDeviation = 0.05, ...) {
  ## @param expectedAbsDeviation How much do we expect 2 clusters to differ on
  ##   in mean Cramér's V. BIC-based clustering tends to split rather than lump
  ##   and this is a heuristical attempt to forestall this.
  method <- match.arg(method, c("mclust", "dunn"))
  indent <- list(...)$indent %||% indentation()
  assert_that(has_attr(dist, "method"))
  assoc.measure <- switch(attr(dist, "method"),
                          "cramer.V" = "Cramér's V",
                          "spearman" = "Spearman's Rho",
                          "kendall" = "Kendall's Tau")
  if (method == "mclust") {
    diag(dist) <- 1
    bic <- mclust::mclustBIC(dist, G = 1:2, verbose = FALSE)
    mc <- mclust::Mclust(dist, x = bic, verbose = FALSE)
    if (mc$G == 2) {
      diag(dist) <- NA_real_
      classification <- mc$classification
      i <- names(which(classification == 1))
      j <- names(which(classification == 2))
      if (length(i) == 1) {
        d <- mean(dist[i, ], na.rm = TRUE) - mean(dist[j, j], na.rm = TRUE)
      } else if (length(j) == 1) {
        d <- mean(dist[i, i], na.rm = TRUE) - mean(dist[j, ], na.rm = TRUE)
      } else {
        d <- mean(dist[i, i], na.rm = TRUE) - mean(dist[j, j], na.rm = TRUE)
      }
      flog.info("%s2 SNP clusters selected that differ by <%0.4f> mean %s",
                indent(), abs(d), assoc.measure, name = "info")
      if (abs(d) > expectedAbsDeviation) {
        if (d > expectedAbsDeviation) {
          clustered.snps <- i
        } else {
          clustered.snps <- j
        }
      } else {
        flog.info("%sRejecting clusters as mean difference does not exceed thhreshold <%s>",
                  indent(), expectedAbsDeviation, name = "info")
        clustered.snps <- c(i, j)
      }
    }
    else if (mc$G == 1) {
      clustered.snps <- mc$classification
    }
    attr(clustered.snps, "mclustBIC") <- bic
    attr(clustered.snps, "mclust") <- mc
  }
  else if (method == "dunn") {
    ## calculate mean association of each PP too each other PP
    dunnCutoff <- list(...)$dunnCutoff %||% 18
    minimumMeanAssociation <- list(...)$minimumMeanAssociation %||% 0.12
    cm <- .colMeans(dist, m = NCOL(dist), n = NROW(dist), na.rm = TRUE)
    names(cm) <- colnames(dist)
    ## try to split PP in a low association and a high association cluster.
    ## evaluate using a Dunn-like index, the ratio of the largest distance
    ## between the sorted average association (i.e., the smallest distance
    ## between the two putative clusters) to the mean average association
    ## (approximately the mean intra-cluster distance).
    scm <- sort(cm)  # ranked mean association
    dcm <- diff(scm) # ranked difference in mean association
    c1 <- scm[1:which.max(dcm)] ## putative low association cluster
    c2 <- scm[-(1:which.max(dcm))] ## putative high association cluster
    dunnIndex <- max(dcm)/mean(dcm)
    if (dunnIndex > dunnCutoff) {
      clustered.snps <- names(c2)
      flog.info("%sInferring cluster of low association polymorphisms at Dunn index <%0.1f>",
                indent(), dunnIndex, name = "info")
      flog.info("%sRemoving %s low-association polymorphic positions ", indent(), length(c1), name = "info")
    } else {
      clustered.snps <- names(cm)
      flog.info("%sInferring no cluster of low association polymorphisms at Dunn index <%0.1f>",
                indent(), dunnIndex, name = "info")
    }
    ## use minimumMeanAssociation as a hard cutoff
    if (any(idx <- cm[clustered.snps] >= minimumMeanAssociation)) {
      clustered.snps <- names(which(idx))
      flog.info("%sRemoving %s additional polymorphic positions below minimum mean association threshold <%0.2f>", indent(),
                sum(!idx), minimumMeanAssociation, name = "info")
    }
  }

  clustered.snps
}

.correlogram <- function(dist, nm) {
  ## @param nm Names of clustered SNPs to be highlighted in the association
  ##   plot as selected.
  method.assoc <- attr(dist, "method")
  ## Set up labels
  if (method.assoc == "cramer.V") {
    y.lab <- "Mean Cramér's V"
    legend.label <- "V"
  } else if (method.assoc == "spearman") {
    y.lab <- "Mean Spearman's Rho"
    legend.label <- "rho"
  } else if (method.assoc == "kendall") {
    y.lab <- "Mean Kendall's Tau"
    legend.label <- "tau"
  }
  ## Correlogram
  cl <- stats::hclust(stats::as.dist(1 - dist), method = "ward.D")
  ocmat <- dist[cl$order, cl$order]
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

  ## Association plot
  cm <- colMeans(dist, na.rm = TRUE)
  rDf <- tibble::tibble(
    Position = factor(names(cm), levels = names(cm), labels = names(cm), ordered = TRUE),
    meanAssoc = unname(cm),
    picked = ifelse(names(cm) %in% nm, "yes", "no")
  )
  p2 <- ggplot(rDf, aes(x = Position, y = meanAssoc, colour =  picked)) +
    geom_point() +
    expand_limits(y = 0) +
    labs(x = "Polymorphic position", y = y.lab) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
          legend.position = "bottom")

  invisible(list(correlogram = p1, association = p2))
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

findCutpoint <- function(d, m) {
  c1 <- m$posterior[, "comp.1"]
  c2 <- m$posterior[, "comp.2"]
  d1 <- d[which(c1 > c2)]
  d2 <- d[which(c1 < c2)]
  lower <- which.min(c(mean(d1, na.rm = TRUE), mean(d2, na.rm = TRUE)))
  if (lower == 1) {
    (max(d1) + min(d2))/2
  } else {
    (min(d1) + max(d2))/2
  }
}

