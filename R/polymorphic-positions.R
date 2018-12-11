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
                                                  minimumExpectedDifference = 0.06,
                                                  noSelect = FALSE,
                                                  resample = NULL,
                                                  ...) {
  method.assoc <- match.arg(method.assoc, c("cramer.V", "spearman", "kendall"))
  method.clust <- match.arg(method.clust, c("mclust", "dunn"))
  indent <- list(...)$indent %||% indentation()
  #resample <- floor(NCOL(mat)*0.1)
  cmat <- .associationMatrix(mat, method.assoc, resample, ...)
  if (noSelect) {
    selected.snps <- structure(colnames(mat), classification = "1")
  } else {
    selected.snps <- .clusterPolymorphicPositions(
      dist = cmat, method = method.clust, minimumExpectedDifference = minimumExpectedDifference,
      ..., indent = indent)
  }
  flog.info("%sRetaining %s high-association polymorphisms", indent(), length(selected.snps), name = "info")
  ## create correlogram and association plots
  plts <- .correlogram(dist = cmat, nm = selected.snps)
  structure(selected.snps, snp.corr.mat = cmat,
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

.clusterPolymorphicPositions <- function(dist, method = "mclust", minimumExpectedDifference = 0.05, ...) {
  ## @param minimumExpectedDifference How much do we expect 2 clusters to differ on
  ##   in mean Cramér's V. BIC-based clustering tends to split rather than lump
  ##   and this is a heuristical attempt to forestall this.
  method <- match.arg(method, c("mclust", "dunn"))
  indent <- list(...)$indent %||% indentation()

  if (method == "mclust") {
    diag(dist) <- 1
    bic <- mclust::mclustBIC(dist, G = 1:2, verbose = FALSE)
    mc <- mclust::Mclust(dist, x = bic, verbose = FALSE)
    if (mc$G == 2) {
      selected.snps <- .checkClusterMerge(dist, mc, minimumExpectedDifference, indent)
    }
    else if (mc$G == 1) {
      selected.snps <- structure(names(mc$classification), classification = as.character(mc$classification))
    }
    attr(selected.snps, "mclustBIC") <- bic
    attr(selected.snps, "mclust") <- mc
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
      selected.snps <- names(c2)
      flog.info("%sInferring cluster of low association polymorphisms at Dunn index <%0.1f>",
                indent(), dunnIndex, name = "info")
      flog.info("%sRemoving %s low-association polymorphic positions ", indent(), length(c1), name = "info")
    } else {
      selected.snps <- names(cm)
      flog.info("%sInferring no cluster of low association polymorphisms at Dunn index <%0.1f>",
                indent(), dunnIndex, name = "info")
    }
    ## use minimumMeanAssociation as a hard cutoff
    if (any(idx <- cm[selected.snps] >= minimumMeanAssociation)) {
      selected.snps <- names(which(idx))
      flog.info("%sRemoving %s additional polymorphic positions below minimum mean association threshold <%0.2f>", indent(),
                sum(!idx), minimumMeanAssociation, name = "info")
    }
  }

  selected.snps
}

.checkClusterMerge <- function(dist, mc, ...) {
  ## @param ...minimumExpectedDifference The absolute difference
  ##   in mean association between clusters that must be exceded
  ##   to accept the clusters as a secondary check performed after
  ##   the equivalence test
  assert_that(has_attr(dist, "method"))
  indent <- list(...)$indent %||% indentation()
  minimumExpectedDifference <- list(...)$minimumExpectedDifference %||% 0.05
  assoc.measure <- switch(attr(dist, "method"),
                          "cramer.V" = "Cramér's V",
                          "spearman" = "Spearman's Rho",
                          "kendall"  = "Kendall's Tau")
  ##
  diag(dist) <- NA_real_
  classification <- mc$classification
  i <- names(which(classification == 1))
  j <- names(which(classification == 2))
  ## infer mean and standard deviation for the high-association (h)
  ## and the low-association (l) cluster
  if (length(i) == 1) {
    mi  <- mean(dist[i, ], na.rm = TRUE)
    sdi <- sd(dist[i, ], na.rm = TRUE)
  } else {
    mi  <- mean(dist[i, i], na.rm = TRUE)
    sdi <- sd(dist[i, i], na.rm = TRUE)
  }
  if (length(j) == 1) {
    mj  <- mean(dist[j, ], na.rm = TRUE)
    sdj <- sd(dist[j, ], na.rm = TRUE)
  } else {
    mj  <- mean(dist[j, j], na.rm = TRUE)
    sdj <- sd(dist[j, j], na.rm = TRUE)
  }
  ## calculate average distance (dij) between clusters
  dij <- mi - mj
  ## perform an ad hoc equivalence test on the two clusters: calculate the
  ## lower confidence bound of the high-association cluster.
  ## If this bound overlaps the upper confidence bound of the low-association
  ## cluster, reject the clusters
  if (dij > 0) { # i is high-association, j is low-association
    #(mi - 1.96*sdi) <= (mj + 1.96*sdj)
    reject <- (mi - sdi) <= (mj + sdj)
    ovl <- (mj + sdj) - (mi - sdi)
  } else {
    #(mj - 1.96*sdj) <= (mi - 1.96*sdi)
    reject <- (mj - sdj) <= (mi - sdi)
    ovl <- (mi + sdi) - (mj - sdj)
  }
  if (reject) {
    flog.info("%sReject SNP clusters: mean difference <%0.3f> %s with <%0.4f> overlapping confidence bounds",
              indent(), abs(dij), assoc.measure, ovl, name = "info")
    selected.snps <- c(i, j)
  }
  else {
    ## perform a secondary test to see if dij is less than the minimum expected
    ## absolute difference between mi and mj
    if (abs(dij) < minimumExpectedDifference) {
      flog.info("%sReject SNP clusters: mean difference <%0.3f> %s does not exceed threshold <%s>",
                indent(), abs(dij), assoc.measure, minimumExpectedDifference, name = "info")
      selected.snps <- c(i, j)
    }
    else {
      ## perform yet another test to see if the putative high-association SNPs
      ## are likely to differ by a location shift within the gene. If this happens,
      ## it is quite likely that one cluster is a local high linkage block.
      if (.locationShift(i, j, p.value = 0.001)) {
        flog.info("%sReject SNP clusters: significant location shift.",
                  indent(), name = "info")
        selected.snps <- c(i, j)
      }
      else {
        flog.info("%sAccept SNP clusters with mean difference <%0.3f> %s",
                  indent(), abs(dij), assoc.measure, name = "info")
        if (dij > 0) {
          selected.snps <- i
        } else {
          selected.snps <- j
        }
      }
    }
  }

  o1 <- order(as.numeric(selected.snps))
  o2 <- order(as.numeric(c(i, j)))
  structure(selected.snps[o1], classification = c(rep("1", length(i)), rep("2", length(j)))[o2])
}

.locationShift <- function(i, j, p.value = 0.001) {
  i  <- as.numeric(i)
  j  <- as.numeric(j)
  ni <- length(i)
  nj <- length(j)
  ## we test the hypothesis that a randomly selected location from set i will be less
  ## or greater than a rondomly selected location from set j (Wilcoxon rank-sum test)
  rij <- rank(c(i, j))
  Ri <- sum(rij[1:ni])
  Rj <- sum(rij[(ni + 1):(ni + nj)])
  Wi <- Ri - (ni*(ni + 1))/2
  Wj <- Rj - (nj*(nj + 1))/2
  W  <- min(Wi, Wj)
  2*stats::pwilcox(W, ni, nj) < p.value
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
    selected = ifelse(names(cm) %in% nm, "yes", "no"),
    cluster = attr(nm, "classification") %||% 1
  )
  p2 <- ggplot(rDf, aes(x = Position, y = meanAssoc, colour =  selected, shape = cluster)) +
    geom_point() +
    expand_limits(y = 0) +
    labs(x = "Polymorphic position", y = y.lab) +
    theme_bw() +
    scale_colour_manual(values = c("#4c8cb5", "#e15a53")) +
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

