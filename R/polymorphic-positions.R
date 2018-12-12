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
                                                  measureOfAssociation = "cramer.V",
                                                  proportionOfOverlap = 1/3,
                                                  minimumExpectedDifference = 0.06,
                                                  noSelect = FALSE,
                                                  resample = NULL,
                                                  ...) {
  measureOfAssociation <- match.arg(measureOfAssociation, c("cramer.V", "spearman", "kendall"))
  indent <- list(...)$indent %||% indentation()
  #resample <- floor(NCOL(mat)*0.1)
  dist <- .associationMatrix(mat, measureOfAssociation, resample)#, ...)
  if (noSelect) {
    selected.snps <- structure(colnames(mat), classification = "1")
  } else {
    selected.snps <- .clusterPolymorphicPositions(
      dist, proportionOfOverlap, minimumExpectedDifference, indent = indent)#, ...)
  }
  flog.info("%sRetaining %s high-association polymorphisms", indent(), length(selected.snps), name = "info")
  ## create correlogram and association plots
  plts <- .correlogram(dist, nm = selected.snps)
  structure(selected.snps, snp.corr.mat = dist,
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

.clusterPolymorphicPositions <- function(dist,
                                         proportionOfOverlap = 0.1,
                                         minimumExpectedDifference = 0.05,
                                         ...) {
  ## @param proportionOfOverlap We perform an equivalence test on the two clusters:
  ##   calculate the lower 1-sigma bound of the high-association cluster i.
  ##   calculate the upper 1-sigma bound of the low-association cluster j.
  ##   reject clusters, if this bounds overlap by more than <proportionOfOverlap> of
  ##   the average distance (dij) between clusters.
  ## @param minimumExpectedDifference The absolute difference in mean association
  ##   between clusters that must be exceeded to accept the clusters as a secondary
  ##   check performed after the equivalence test
  indent <- list(...)$indent %||% indentation()
  diag(dist) <- 1
  bic <- mclust::mclustBIC(dist, G = 1:2, verbose = FALSE)
  mc  <- mclust::Mclust(dist, x = bic, verbose = FALSE)
  if (mc$G == 2) {
    selected.snps <- .checkClusterMerge(dist, mc, proportionOfOverlap, minimumExpectedDifference, indent)
  }
  else if (mc$G == 1) {
    selected.snps <- structure(
      names(mc$classification), classification = as.character(mc$classification))
  }
  attr(selected.snps, "mclustBIC") <- bic
  attr(selected.snps, "mclust") <- mc
  selected.snps
}

.checkClusterMerge <- function(dist, mc, ...) {
  ## @param proportionOfOverlap We perform an equivalence test on the two clusters:
  ##   calculate the lower 1-sigma bound of the high-association cluster i.
  ##   calculate the upper 1-sigma bound of the low-association cluster j.
  ##   reject clusters, if this bounds overlap by more than <proportionOfOverlap> of
  ##   the average distance (dij) between clusters.
  ## @param ...minimumExpectedDifference The absolute difference in mean association
  ##   between clusters that must be exceeded to accept the clusters as a secondary
  ##   check performed after the equivalence test
  assert_that(has_attr(dist, "method"))
  indent <- list(...)$indent %||% indentation()
  proportionOfOverlap <- list(...)$proportionOfOverlap %||% 0.1
  minimumExpectedDifference <- list(...)$minimumExpectedDifference %||% 0.05
  measureOfAssociation <- switch(attr(dist, "method"),
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
    mi  <- mean(dist[i, i][lower.tri(dist[i, i])])
    sdi <- sd(dist[i, i][lower.tri(dist[i, i])]) %|na|% 0.001
  }
  if (length(j) == 1) {
    mj  <- mean(dist[j, ], na.rm = TRUE)
    sdj <- sd(dist[j, ], na.rm = TRUE)
  } else {
    mj  <- mean(dist[j, j][lower.tri(dist[j, j])])
    sdj <- sd(dist[j, j][lower.tri(dist[j, j])]) %|na|% 0.001
  }
  ## mi: mean association of high-association cluster
  ## mj: mean association of low-association cluster
  if (mi < mj) { ## mi <-> mj; i <-> j
    tmp.mj  <- mi
    tmp.sdj <- sdi
    tmp.j   <- i
    mi  <- mj
    sdi <- sdj
    i   <- j
    mj  <- tmp.mj
    sdj <- tmp.sdj
    j   <- tmp.j
  }
  ## perform an ad hoc equivalence test on the two clusters
  ovl <- .OVL(mi, sdi, mj, sdj, limitLevel = 1, proportionOfOverlap = proportionOfOverlap)
  if (ovl$reject) {
    flog.info("%sReject SNP clusters: mean difference <%0.3f> %s with <%0.2f%%> overlapping 1-sigma limits",
              indent(), ovl$dij, measureOfAssociation, 100*ovl$ovl, name = "info")
    selected.snps <- c(i, j)
  } else {
    ## perform a secondary test to see if dij is less than the minimum expected
    ## absolute difference between mi and mj
    if (ovl$dij < minimumExpectedDifference) {
      flog.info("%sReject SNP clusters: mean difference <%0.3f> %s does not exceed threshold <%s>",
                indent(), ovl$dij, measureOfAssociation, minimumExpectedDifference, name = "info")
      selected.snps <- c(i, j)
    } else {
      ## perform yet another test to see if the putative high-association SNPs
      ## are likely to differ by a location shift within the gene. If this happens,
      ## it is quite likely that one cluster is a local high linkage block.
      if (.locationShift(i, j, p.value = 0.001)) {
        flog.info("%sReject SNP clusters: significant location shift.",
                  indent(), name = "info")
        selected.snps <- c(i, j)
      } else {
        flog.info("%sAccept SNP clusters with mean difference <%0.3f> %s",
                  indent(), ovl$dij, measureOfAssociation, name = "info")
        selected.snps <- i
      }
    }
  }

  o1 <- order(as.numeric(selected.snps))
  o2 <- order(as.numeric(c(i, j)))
  structure(selected.snps[o1],
            classification = c(rep("1", length(i)), rep("2", length(j)))[o2],
            dij = ovl$dij, ovl = ovl$ovl, ovl.coef = ovl$ovl.coef, ovl.plot = ovl$ovl.plot)
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
  measureOfAssociation <- attr(dist, "method")
  ## Set up labels
  if (measureOfAssociation == "cramer.V") {
    y.lab <- "Mean Cramér's V"
    legend.label <- "V"
  } else if (measureOfAssociation == "spearman") {
    y.lab <- "Mean Spearman's Rho"
    legend.label <- "rho"
  } else if (measureOfAssociation == "kendall") {
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

.OVL <- function(mi, sdi, mj, sdj, limitLevel = 1, proportionOfOverlap = 0.1) {
  min.fifj <- function(x, mi, mj, sdi, sdj) {
    fi <- dnorm(x, mean = mi, sd = sdi)
    fj <- dnorm(x, mean = mj, sd = sdj)
    pmin(fi, fj)
  }
  ## expected percentage by which each sample distribution overlaps the other
  rs <- stats::integrate(min.fifj, -Inf, Inf, mi, mj, sdi, sdj)
  ## calculate average distance (dij) between clusters
  dij <- mi - mj
  ## perform an ad hoc equivalence test on the two clusters:
  ## calculate the lower 1sigma bound of the high-association cluster.
  ## calculate the upper 1sigma bound of the low-association cluster.
  ## reject, if this bounds overlap by more than 10% of the average distance
  ## (dij) between clusters.
  ovl <- ((mj + limitLevel*sdj) - (mi - limitLevel*sdi))/dij
  reject <- ovl > proportionOfOverlap
  ## plot
  xs <- seq(mj - 3*sdj, mi + 3*sdi, 0.01)
  f1 <- dnorm(xs, mean = mi, sd = sdi)
  f2 <- dnorm(xs, mean = mj, sd = sdj)
  ys <- min.fifj(xs, mi, mj, sdi, sdj)
  xs2 <- c(xs, xs[1])
  ys2 <- c(ys, ys[1])
  grDevices::pdf(NULL, bg = "white")
  grDevices::dev.control(displaylist = "enable")
  plot(xs, f1, type = "n", ylim = c(0, max(f1, f2)),
       xlab = "association", ylab = "density", bg = "white")
  lines(xs, f1, lty = "dashed", lwd = 2)
  lines(xs, f2, lty = "dotted", lwd = 2)
  polygon(xs2, ys2, col = "gray80", density = 20)
  abline(v = c(mi - limitLevel*sdi), lty = "dashed", col = "red")
  abline(v = c(mj + limitLevel*sdj), lty = "dotted", col = "red")
  abline(v = mi, lty = "dashed", col = "green")
  abline(v = mj, lty = "dotted", col = "green")
  p.base <- grDevices::recordPlot()
  invisible(grDevices::dev.off())

  list(reject = reject, dij = dij, ovl = ovl, ovl.coef = rs$value, ovl.plot = p.base)
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

