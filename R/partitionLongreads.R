# partitionReads -----------------------------------------------------------

#' Cluster a SNP matrix into haplotype groups
#'
#' @param x A \code{[read x position]} SNP matrix.
#' @param distAlleles Number of distinct alleles in the sample.
#' @param sortBy sort by "distance" or "count". See
#' \code{\link{partitionLongreads}} for details
#' @param threshold  Only gaps above threshold will be treated as gaps.
#' Removes noisy positions in longreads.
#' @param clMethod clustering method passed to \code{\link[stats]{hclust}}.
#' @param minLen Minimal fraction of SNPs that need to be covered by a read
#' to be used
#' @param skipGapFreq Skip a badly behaved polymorphic position if the
#' frequency of gaps within a haplotype group exceeds \code{skipGapFreq}.
#' @param deepSplit sensitivity parameter passed to
#' \code{\link[dynamicTreeCut]{cutreeHybrid}}.
#' @param ...  Additional parameters.
#'
#' @return A \code{HapPart} object.
#' @export
#' @examples
#' ###
partitionReads <- function(x, distAlleles = 2, selectAllelesBy = "count", threshold = 0.2,
                           clMethod = "ward.D", minLen = 0.5, skipGapFreq = 0.8,
                           deepSplit = 1, minClusterSize = 15, ...) {

  indent <- list(...)$indent %||% indentation()
  selectAllelesBy <- match.arg(selectAllelesBy, c("count", "distance"))
  # get SNPs
  ppos <- colnames(x)
  badPpos <- c()
  xm <- x[order(rownames(x)), , drop = FALSE]
  ## if there is only one SNP for clustering, use it! If it does not match both
  ## sequencing types it will be reported
  if (length(ppos) > 1) {
    badPpos <- apply(xm, 2, function(col) {
      sum(col == "+" | col == "-")/NROW(col) > skipGapFreq
    })
    xm <- xm[, !badPpos, drop = FALSE]
    badPpos <- ppos[badPpos]
    flog.info("%s%s polymorphisms occupy less than %0.3g%% of the reads and are discarded",
              indent(), length(badPpos), 100*(1 - skipGapFreq), name = "info")
    flog.info("%sUsing the remaining %s polymorphisms for clustering",
              indent(), NCOL(xm), name = "info")
    if (NCOL(xm) == 0) {
      flog.warn("%s No polymorphisms remain. Aborting clustering." %<<%
                 " Check reads and reference and inspect the mapInit plot",
                 indent(), name = "info")
      part_ <- HapPart(readNames = row.names(xm), snpPos = colnames(xm))
      PRT(part_) <- "A"  ## <character>; the haplotype assignment for each read
      return(part_)
    }
  }

  ## Get the SNP matrix as sequences
  xseqs <- .getSeqsFromMat(xm)
  ## if only one SNP, remove every read without it
  if (NCOL(xm) == 1)
    xseqs <- xseqs[!xseqs == "+"]

  ## Get only the fraction of reads that contain at least minLen of total SNPs
  clusters <- .getClusts(xseqs, clMethod = clMethod, minLen = minLen,
                         deepSplit = deepSplit, threshold = threshold,
                         minClusterSize = minClusterSize, suppressGaps = FALSE,
                         indent = indent)
  if ("@" %in% levels(clusters$clades))
    flog.warn("%s Removed one clusters due to small cluster size of %s. To force processing run with option minClusterSize < %s.",
                 indent(), sum(clusters$clades == "@"), sum(clusters$clades == "@"), name = "info")
  subclades <- factor(clusters$clades[!clusters$clades == "@"])
  tree <- clusters$tree
  hptypes <- levels(subclades)

  if (length(subclades) == 0) {
    flog.error("%sTwo few longreads for clustering", indent(), name = "info")
    stop("To few longreads for clustering")
  }

  flog.info("%sInitial clustering results in %s haplotypes <%s>",
            indent(), length(hptypes), comma(hptypes), name = "info")

  ## Get scores and assign clades by freq mat
  ## Position Weight Matrix: Use frequency plus pseudocount/ basefrequency
  ## (here 0.25 for each).
  msa <- stats::setNames(lapply(levels(subclades), function(x) {
    xseqs[names(subclades[subclades == x])]
  }), hptypes)
  mats <- lapply(msa, function(x) createPWM(x))
  hpseqs <- Biostrings::DNAStringSet(vapply(msa, function(x) {
    conseq(t(Biostrings::consensusMatrix(x)[c(VALID_DNA(), "+"), ]), type = "simple")
  }, FUN.VALUE = character(1)))
  names(hpseqs) <- vapply(hptypes, function(hptype) {
    colon(c(hptype, sum(subclades == hptype)))
  }, FUN.VALUE = character(1))

  if (length(hptypes) > distAlleles) {
    if (distAlleles == 1) {
      rC <- "A"
    } else {
      flog.info("%sIdentify chimeric reads/haplotypes", indent(), name = "info")
      if (selectAllelesBy == "count") {
        rC <- names(sort(table(subclades), decreasing = TRUE)[seq_len(distAlleles)])
      }
      else if (selectAllelesBy == "distance") {
       #  rC <- sort(.findChimeric(seqs = hpseqs, distAlleles = distAlleles))
        rC <- .findMostDistantPwm(mats)
      }
    }
    flog.info("%sUse only clusters <%s>", indent(), comma(rC), name = "info")
    mats <- mats[rC]
  }
  if (length(hptypes) > 1) {
    maxDeltaScore <- .findMostDistantPwm(mats, onlyScore = TRUE)
    flog.info("Maximum score difference of the used haplotypes: %s", max(maxDeltaScore))
  }

  hptypes <- names(mats)
  scores <- .getScores(xseqs, mats)
  clades <- scores %>%
    dplyr::group_by(.data$read) %>%
    dplyr::slice(which.max(.data$score)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      clustclade = dplyr::if_else(.data$read %in% names(subclades),
                                  as.character(subclades[.data$read]),
                                  NA_character_)) %>%
    dplyr::mutate(
      correct = dplyr::if_else(.data$clustclade == .data$clade, TRUE, FALSE))

  ## Correctly classified clades in the initial clustering
  falseClassified <- sum(!clades$correct, na.rm = TRUE)/sum(clades$correct, na.rm = TRUE)
  flog.info("%sCorrected classification of %.2f%% of reads", indent(), 100*falseClassified, name = "info")
  lapply(hptypes, function(hp, clades) {
    invisible(flog.info("%s%s reads in haplotype <%s>", indent(), table(clades$clade)[hp], hp, name = "info"))
   }, clades = clades)

  # Create the partition table
  part_ <- HapPart(readNames = clades$read, snpPos = colnames(xm))
  PRT(part_)    <- clades$clade  ## <character>; the haplotype assignment for each read
  mcoef(part_)  <- clades$score  ## <numeric>; the membership coefficient for each read
  HTR(part_)    <- tree          ## <hclust>; hierarchical cluster tree
  K(part_)      <- length(ppos) - length(badPpos) ## <numeric>; The total number of polymorphic positions used.
  OC(part_)     <- levels(subclades) ## <character>; The original clusters
  SCR(part_)    <- scores
  SQS(part_)    <- hpseqs
  PWM(part_)    <- mats
  ##
  attr(part_, "snp.corr.mat") <- attr(x, "snp.corr.mat")
  attr(part_, "snp.clust") <- attr(x, "snp.clust")

  return(part_)
}



# Helpers -----------------------------------------------------------------

.getClusts <- function(xseqs, minLen = 0.80, clMethod = "ward.D",
                       deepSplit = 1, minReadsFrac = 1/3, threshold = 0.2,
                       minClusterSize = 15, suppressGaps = TRUE, ...) {
  assert_that(
    is.double(minLen),
    is.double(minReadsFrac),
    is.double(threshold)
  )
  indent <- list(...)$indent %||% indentation()
  # heuristic for how many reads should be used
  minLen <- .getMinLenClust(minReadsFrac, minLen, xseqs, indent = indent)
  flog.info("%sUse only longreads containing at least %s%% of all polymorphisms",
            indent(), minLen*100, name = "info")
  xSub <- xseqs[Biostrings::width(gsub("\\+", "", xseqs)) > minLen*unique(Biostrings::width(xseqs))]
  ## Consensus matrix with pseudocount
  flog.info("%sConstruct a Position Specific Distance Matrix" %<<%
            " from the remaining %s sequences", indent(), length(xSub), name = "info")
  consmat <- as.matrix(
    Biostrings::consensusMatrix(xSub, as.prob = TRUE, baseOnly = FALSE)[c(VALID_DNA(), "+" ), ] +
      1/length(xSub))
  if (suppressGaps) {
    # Remove gaps below a threshold as they are probably sequencing artifacts
    consmat <- apply(consmat, 2, function(a) {
      a["-"] <- ifelse(a["-"] < threshold, min(a), a["-"])
      a
    })
  }
  ## normalise to sum to 1
  consmat <- consmat/colSums(consmat)
  ## Get Position Specific Distance Matrix
  dist <- PSDM(x = xSub, consmat = as.matrix(consmat))
  dist <- stats::as.dist(dist)
  ## replace "na" with mean for a being able to cluster
  dist[is.na(dist)] <- mean(dist, na.rm = TRUE)

  ## Perform a hierarchical clustering
  hcc <-  stats::hclust(dist, method = clMethod)

  ## do a dynamic cut. Need to be evaluated
  clusts <- dynamicTreeCut::cutreeHybrid(hcc,
                                         minClusterSize = minClusterSize,
                                         distM = as.matrix(dist),
                                         deepSplit = deepSplit,
                                         verbose = FALSE)
  # extract the clusters for each sequence
  clades <- as.factor(vapply(clusts$labels, function(i) {
    rawToChar(as.raw(as.integer(i) + 64))
  }, FUN.VALUE = character(1), USE.NAMES = FALSE))
  names(clades) <- hcc$labels
  return(list(clades = clades, tree = hcc))
}

## Get minimal fraction of SNPs needed for a read to be used
.getMinLenClust <- function(minReadsFrac, minLen, xseqs, ...) {
  indent <- list(...)$indent %||% indentation()
  adjustedMinReadsFrac <- minReadsFrac +
    (0.03*Biostrings::width(xseqs[1])/50) * 1500^2/length(xseqs)^2
  adjustedMinReadsFrac <- min(0.5, adjustedMinReadsFrac)
  flog.info("%sAdjust the minimal fraction of reads used for clustering" %<<%
            " from %0.4g%% to %0.4g%%", indent(), minReadsFrac*100,
            adjustedMinReadsFrac*100, name = "info")
  # get long reads containing sufficient number of SNPs
  minLens <- seq(1, 0, -0.01)
  lenCounts <- vapply(minLens, function(minl) {
    length(xseqs[
      Biostrings::width(gsub("\\+", "", xseqs)) >
        minl*Biostrings::width(xseqs[1])
      ])/length(xseqs)
  }, FUN.VALUE = numeric(1))

  names(lenCounts) <- minLens
  minLen <- max(as.numeric(names(
    lenCounts[which(lenCounts > adjustedMinReadsFrac)][1])), minLen)
}

.findMostDistantPwm <- function(mats, onlyScore = FALSE) {
  ## Assert that all pwm are of same length
  assert_that(length(unique(vapply(mats, NCOL, numeric(1)))) == 1)
  combs <- combn(names(mats), 2, FUN = paste0)
  if (onlyScore) {
    dMats <- apply(combs, 2, function(comb) {
      .pwmDiff(mats[[comb[[1]]]], mats[[comb[[2]]]], max = TRUE)
    })
    names(dMats) <- apply(combs, 2, function(x) paste0(x, collapse = ""))
    return(dMats)
  }
  dMats <- apply(combs, 2, function(comb) {
    .pwmDiff(mats[[comb[[1]]]], mats[[comb[[2]]]], max = FALSE)
  })
  names(dMats) <- apply(combs, 2, function(x) paste0(x, collapse = ""))
  ## Use only pairs with an A.
  dMats <- dMats[grepl("A", names(dMats))]
  strsplit(names(which.max(dMats)), "")[[1]]
}
.pwmDiff <- function(m1, m2, max = FALSE) {
  ## Check if excluding gaps helps
  # m1 <- m1[VALID_DNA(include = "none"),]
  # m2 <- m2[VALID_DNA(include = "none"),]
  scores <- vapply(seq_len(NCOL(m1)), function(i, m1, m2) {
    sum((m1[,i] - m2[,i])^2)
  }, m1 = m1, m2 = m2, FUN.VALUE = numeric(1))
  if (max) 
    return(max(scores))
  sum(scores)/NCOL(m1)
}

.getScores <- function(reads, pwmlist) {
  assert_that(
    is(reads, "XStringSet"),
    is.list(pwmlist),
    !is.null(names(pwmlist))
  )
  rs <- do.call(dplyr::bind_rows, Map(function(pwm, hp) {
    b <- vapply(reads, function(read) .read_score(read, pwm), FUN.VALUE = double(1))
    tibble::tibble(read = names(b), score = b, clade = hp)
  }, pwm = pwmlist, hp = names(pwmlist))) %>%
    dplyr::arrange(.data$read, dplyr::desc(.data$score))
  rs
}

## Identify chimeric clusters
.findChimeric <- function(seqs, distAlleles, plotSeqs = FALSE) {
  # Use only non empty seqs
  seqs <- seqs[vapply(seqs, function(x) {
    !nchar(stripIndel(x)) == 0
  }, FUN.VALUE = logical(1))]

  forward <- lapply(seq_len(Biostrings::width(seqs[1])), function(x) {
    hammingDist(Biostrings::DNAStringSet(seqs, start = 1, width = x))})
  reverse <- lapply(seq_len(Biostrings::width(seqs[1])), function(x) {
    hammingDist(Biostrings::DNAStringSet(Biostrings::reverse(seqs),
                                         start = 1, width = x))})

  dm <- foreach(h = names(seqs), .combine = 'cbind') %:%
    foreach(hp = names(seqs)) %do% {
      f <- vapply(forward, function(x) x[h, hp], FUN.VALUE = numeric(1))
      r <- vapply(reverse, function(x) x[h, hp], FUN.VALUE = numeric(1))
      unname(stats::quantile(f + r)[2])
    }
  nms <- vapply(strsplit(names(seqs), ":"), function(x) x[1],
                FUN.VALUE = character(1))
  rownames(dm) <- colnames(dm) <- nms
  c("A", (names(sort(unlist(dm[1,]),
                     decreasing = TRUE)[seq_len(distAlleles - 1)])))
}

# Class: HapPart -----------------------------------------------------------

#' HapPart (Haplotype Partition) stores the haplotype information from
#' longread clustering.
#' @param readNames The names of reads in each cluster.
#' @param snpPos SNP positions used for clustering.
#' @param x A \code{\link{HapPart}} object.
#' @return HapPart object.
#' @slot readNames The names of reads in each cluster.
#' @slot snpPos SNP positions used for clustering.
#' @slot mcoef membership coefficient. Score for each read for belonging
#' to its cluster. Based on the sum of probabilities of a PWM of the haplotype
#' sequences.
#' @slot tree The resulting tree from \code{\link[stats]{hclust}}.
#' @slot scores The same scores as in mcoef, but for all clusters.
#' @slot mats PWM matrices of the haplotype sequences.
#' @slot classification ... TODO
#' @export
HapPart <- function(readNames, snpPos) {
  readNames <- as.character(readNames)
  snpPos   <- as.integer(snpPos)
  n <- length(readNames) # number of reads
  k <- length(snpPos)    # number of polymorphic positions
  structure(
    readNames,
    snpPos = snpPos,
    k      = 0L,            # total number of polymorphic positions
    mcoef  = rep(0, n),
    tree   = NULL,          # Add tree from hclust
    scores = NULL,
    mats   = NULL,
    classification = matrix( # running classification of polymorphic positions
      rep("", n * k),
      nrow = n,
      ncol = k,
      dimnames = list(NULL, snpPos)
    ),
    partition = rep("", n),  # final partitioning of reads into haplotypes
    class = c("HapPart", "character")
  )
}


## Methods: HapPart --------------------------------------------------------

#' @export
print.HapPart <- function(x, sortBy = "none", nrows = 8, ...) {
  sortBy <- match.arg(sortBy, c("none", "name", "mcoef"))
  n <- max(length(x), nrows)
  df0 <- data.frame(
    read      = substr(as.vector(x), 1, 12) %<<% "...",
    mcoef     = mcoef(x),
    partition = PRT(x)
  )[seq_len(n), ]
  cat("HapPart over", K(x), "polymorphic positions:\n")
  switch(
    sortBy,
    none  = print(df0, ...),
    name  = print(df0[order(df0$read),], ...),
    mcoef = print(df0[order(df0$mcoef),], ...)
  )
  invisible(x)
}

#' @export
`[.HapPart` <- function(x, i, j = NULL, drop = NULL) {
  i <- if (is.character(i)) x %in% i else i
  rs <- NextMethod()
  attr(rs, "snpPos") <- SNP(x)
  attr(rs, "k") <- K(x)
  # add mcoef attr
  attr(rs, "mcoef") <- mcoef(x)[i]
  attr(rs, "classification") <- CLS(x)[i, , drop = FALSE]
  attr(rs, "partition") <- PRT(x)[i]
  attr(rs, "tree") <- HTR(x)
  attr(rs, "scores") <- SCR(x)
  class(rs) <- c("HapPart", "character")
  rs
}

mcoef <- function(x) UseMethod("mcoef")
#' @describeIn HapPart
#' The membership coefficient in the assigned clade.
#' @export
mcoef.HapPart <- function(x) {
  attr(x, "mcoef")
}
`mcoef<-` <- function(x, value) UseMethod("mcoef<-")
`mcoef<-.HapPart` <- function(x, value) {
  attr(x, "mcoef") <- value
  x
}

## The PWM matrix of a cluster
PWM <- function(x) UseMethod("PWM")
#' @describeIn HapPart
#' Get the Position Weight Matrix
#' @export
PWM.HapPart <- function(x) {
  attr(x, "PWM")
}
`PWM<-` <- function(x, value) UseMethod("PWM<-")
`PWM<-.HapPart` <- function(x, value) {
  attr(x, "PWM") <- value
  x
}

## Get the original number of clusters for plotting the tree
OC <- function(x) UseMethod("OC")
#' @describeIn HapPart
#' The original clusters from \code{\link[stats]{hclust}}.
#' @export
OC.HapPart <- function(x) {
  attr(x, "oc")
}
`OC<-` <- function(x, value) UseMethod("OC<-")
`OC<-.HapPart` <- function(x, value) {
  attr(x, "oc") <- value
  x
}

## All vs all cluster distances for each read
SCR <- function(x) UseMethod("SCR")
#' @describeIn HapPart
#' The membership scores for each read in each cluster.
#' @export
SCR.HapPart <- function(x) {
  attr(x, "scores")
}

`SCR<-` <- function(x, value) UseMethod("SCR<-")
`SCR<-.HapPart` <- function(x, value) {
  attr(x, "scores") <- value
  x
}

## The actual partition table
partition <- function(x) UseMethod("partition")
partition.HapPart <- function(x) {
  dplyr::arrange(tibble::tibble(
    read      = as.vector(x),
    haplotype = as.factor(PRT(x)),
    mcoef     = mcoef(x)
  ),
  dplyr::desc(mcoef))
}

## Total number of polymorphic positions between haplotypes
K <- function(x) UseMethod("K")
#' @describeIn HapPart
#' The total number of polymorphic positions used.
#' @export
K.HapPart <- function(x) {
  attr(x, "k")
}

`K<-` <- function(x, value) UseMethod("K<-")
`K<-.HapPart` <- function(x, value) {
  attr(x, "k") <- value
  x
}

## The consensus sequence
SQS <- function(x) UseMethod("SQS")
#' @describeIn HapPart
#' Get the consensus sequences of the haplotypes.
#' @export
SQS.HapPart <- function(x) {
  attr(x, "seqs")
}

`SQS<-` <- function(x, value) UseMethod("SQS<-")
`SQS<-.HapPart` <- function(x, value) {
  attr(x, "seqs") <- value
  x
}

## The classification
CLS <- function(x) UseMethod("CLS")
#' @describeIn HapPart
#' Get the classification of the haplotype partitioning.
#' @export
CLS.HapPart <- function(x) {
  attr(x, "classification")
}

`CLS<-` <- function(x, value) UseMethod("CLS<-")
`CLS<-.HapPart` <- function(x, value) {
  attr(x, "classification")[, K(x)] <- value
  x
}

## All SNP positions
SNP <- function(x) UseMethod("SNP")
#' @describeIn HapPart
#' Get the SNP positions.
#' @export
SNP.HapPart <- function(x) {
  attr(x, "snpPos")
}

`SNP<-` <- function(x, value) UseMethod("SNP<-")
`SNP<-.HapPart` <- function(x, value) {
  attr(x, "snpPos") <- value
  x
}

## Get the partition vector
PRT <- function(x) UseMethod("PRT")
#' @describeIn HapPart
#' The haplotype assignment as a character vector <A>, <B>, ... for
#' each read.
#' @export
PRT.HapPart <- function(x) {
  attr(x, "partition")
}

`PRT<-` <- function(x, value) UseMethod("PRT<-")
`PRT<-.HapPart` <- function(x, value) {
  attr(x, "partition") <- value
  x
}

## get tree from hclust
HTR <- function(x) UseMethod("HTR")
#' @describeIn HapPart
#' The cluster tree produced by \code{\link[stats]{hclust}}.
#' @export
HTR.HapPart <- function(x) {
  attr(x, "tree")
}

`HTR<-` <- function(x, value) UseMethod("HTR<-")
`HTR<-.HapPart` <- function(x, value) {
  attr(x, "tree") <- value
  x
}

# Summarise and Plot ------------------------------------------------------

#' Plot distribution of haplotype partitioned reads
#'
#' @param x A \code{HapPart} object.
#' @param label Optional plot label.
#' @param limits Manually provided limits for plotting. Defaults to read limits
#' from the \code{HapPart} object.
#'
#' @return A \code{ggplot} object.
#' @export
#' @examples
#' ###
plotPartitionHistogram <- function(x, label = "", limits = NULL) {
  stopifnot(is(x, "HapPart"))
  if (is.null(limits)) {
    limits <- x$getLimits()
  }
  data <- partition(x) %>%
    dplyr::mutate(mcoef = ifelse(.data$haplotype == "A", mcoef, -1*mcoef))
  ggplot(data) +
    geom_histogram(aes_string(x = 'mcoef', fill = 'haplotype'), bins = 100) +
    scale_fill_manual(values = PARTCOL()) +
    geom_vline(xintercept = c(limits[1], -limits[2]), linetype = "dashed",
               colour = "grey80") +
    xlab("Haplotype membership coefficient") +
    ylab("Number of reads") +
    theme_bw()
}

#' Plot distribution of haplotype partitioned reads for more than two
#' haplotypes
#'
#' @param x A \code{HapPart} object.
#' @param label Optional plot label.
#' @param limits provided limits for each haplotype as a named list.
#'
#' @return A \code{ggplot} object.
#' @export
#' @examples
#' ###
#'
plotPartitionHistogramMulti <- function(x, limits, label = "") {
  stopifnot(is(x, "HapPart"))
  data <- partition(x) %>%
    dplyr::mutate(limit = unlist(limits)[.data$haplotype])

  ggplot(data) +
    geom_histogram(aes(x = mcoef, fill = haplotype), bins = 100) +
    facet_grid(~ haplotype) +
    geom_vline(aes(xintercept = limit), colour = "grey40",
                linetype = "dashed", size = 1) +
    scale_fill_manual(values = PARTCOL()) +
    xlab("Haplotype membership coefficient") +
    ylab("Number of reads") +
    ggtitle(label = label) +
    theme_bw()
}

#' Plot a radar chart for each haplotype. Membership values for each read
#' for each partition are shown as lines.
#'
#' @param x A \code{HapPart} object.
#'
#' @return A \code{ggplot} object.
#' @export
#'
#' @examples
#' ###
plotRadarPartition <- function(x){
  stopifnot(is(x, "HapPart"))

  df <- dplyr::full_join(SCR(x), partition(x), by = "read")

  # use bigger size for 2d radarplot.
  size <- ifelse(length(levels(df$haplotype)) == 2, 4, 0.05)
  ggplot(df, aes_string(x = 'clade', y = 'score')) +
    geom_polygon(aes_string(group = 'read', color = 'haplotype'), fill = NA,
                 size = size, show.legend = FALSE, alpha = 1) +
    facet_grid( ~ haplotype) +
    ggtitle("Similarity to clusters") +
    scale_color_manual(values = PARTCOL()) +
    theme_bw() +
    theme(
      axis.title = element_blank(),
      axis.text.y  = element_blank(),
      axis.ticks  = element_blank()
    ) +
    .coordRadar()
}

#' Plot tree of initial clustering
#'
#' @param x A \code{HapPart} object.
#'
#' @return A \code{ggplot} object.
#' @export
#' @examples
#' ###
plotPartitionTree <- function(x){
  assert_that(
    is(x, "HapPart"),
    requireNamespace("ggdendro", quietly = TRUE)
  )
  tree <- HTR(x)
  k <- length(OC(x))

  tryCatch({
    dendr <- ggdendro::dendro_data(tree, type = "rectangle")
    dendr$labels <- dendr$labels %>%
      dplyr::mutate(label = as.character(.data$label))

    clust <- stats::cutree(tree, k)
    clust <- tibble::tibble(label = names(clust), haplotype = clust)
    dendr$labels <- dplyr::left_join(dendr$labels, clust, by = "label")

    height <- unique(dendr$segments$y)[order(unique(dendr$segments$y),
                                             decreasing = TRUE)]
    cut.height <- mean(c(height[k], height[k - 1]))
    dendr$segments <- dendr$segments %>%
      dplyr::mutate(line =  dplyr::if_else(
        .data$y == .data$yend & .data$y > cut.height, 1, dplyr::if_else(
          .data$yend > cut.height,1, 2)))

    dendr$segments$cluster <- c(-1, diff(dendr$segments$line))
    change <- which(dendr$segments$cluster == 1)
    for (i in seq_len(k)) dendr$segments$cluster[change[i]] = i + 1
    dendr$segments <- dendr$segments %>%
      dplyr::mutate(cluster = dplyr::if_else(.data$line == 1, 1, ifelse(
        .data$cluster == 0, NA, .data$cluster)))
    dendr$segments$cluster <- vapply(seq_len(NROW(dendr$segments$cluster)),
                                     function(x) {
      getCl(x, dendr$segments$cluster, change)}, FUN.VALUE = numeric(1))
    # Correct order
    labs <- c("N", OC(x))
    dendr$segments$cluster <- factor(labs[dendr$segments$cluster])

    # Make plot
    ggplot() +
      geom_segment(data = ggdendro::segment(dendr), aes_string(x = 'x', y = 'y',
                                                               xend = 'xend',
                                                               yend = 'yend',
                                                               color = 'cluster'),
                   size = 1.25) +
      scale_color_manual(values = PARTCOL(),
                         name = "Haplotype",
                         breaks = LETTERS[seq_len(k)],
                         labels =  LETTERS[seq_len(k)]) +
      geom_hline(yintercept = cut.height, color = "blue") +
      ggtitle("Initial clustering of haplotypes") +
      theme_bw() +
      theme(legend.position = "bottom",
            axis.title = element_blank(),
            axis.text  = element_blank(),
            axis.ticks  = element_blank()
      )
  }, error = function(e){
    flog.warn("Too many or too nested nodes for plotting a nice tree.",
              name = "info")
    return(NULL)
  })
}


# Plot helper
getCl <- function(n, cluster, change) {
  ifelse(!is.na(cluster[n]),
         cluster[n],
         cluster[change[max(which(change < n))]]
  )
}
# Helper function for the plotting of radar charts
.coordRadar <- function(theta = "x", start = 0, direction = 1)  {
  theta <- match.arg(theta, c("x", "y"))
  r <- if (theta == "x")
    "y"
  else "x"
  ggproto("CoordRadar", CoordPolar, theta = theta, r = r, start = start,
          direction = sign(direction),
          is_linear = function(coord) TRUE)
}

# Optimal partitioning ----------------------------------------------------


.optimalPartitionLimits <- function(scores, pickiness = 0.8, lowerLimit = 60) {
  ## pickiness > 1: pick more reads
  ## pickiness < 1: pick less reads
  score_range <- range(unlist(scores))
  score_bins <- seq(score_range[1], score_range[2], length.out = 100)
  #coeffCuts <- seq(0, max(unlist(scores)), length.out = 100)
  rHap <- lapply(scores, function(x, score_bins) {
      vapply(score_bins, function(cutoff) sum(x >= cutoff), FUN.VALUE = double(1))
    }, score_bins = score_bins)
  #hp <- "A"
  df <- do.call(dplyr::bind_rows, Map(function(hp) {
    l <- pickiness*length(unlist(scores[names(rHap) != hp]))/length(scores[[hp]])
    tibble::tibble(haplotype = hp, nreads = rHap[[hp]], score = score_bins) %>%
      dplyr::mutate(benefit = ((.data$nreads^l)*.data$score))
  }, hp = names(rHap)))
  dfmax <- df %>%
    dplyr::group_by(haplotype) %>%
    dplyr::filter(.data$benefit == max(.data$benefit))

  if (any(dfmax$nreads < lowerLimit)) {
    hp <- dfmax$haplotype[which(dfmax$nreads < lowerLimit)]
    dfmax2 <- dplyr::filter(df, .data$haplotype %in% hp & .data$nreads >= lowerLimit) %>%
      dplyr::group_by(.data$haplotype) %>%
      dplyr::top_n(1, dplyr::desc(nreads)) %>%
      dplyr::slice(1) ## in case of ties take only the top row
    if (NROW(dfmax2) == 0) {
      flog.warn("No reads above threshold")
      lowerLimit <- dfmax$nreads
    } else {
      dfmax[dfmax$haplotype %in% hp, ] <- dfmax2
    }
  }

  list(
    limits = dfmax,
    plt = ggplot(df, aes_(x = ~nreads, y = ~score, colour = ~haplotype)) +
      geom_line() +
      geom_point(data = dfmax, size = 3) +
      ylab("Cluster similarity coefficient") +
      xlab("N reads") +
      theme_bw()
  )
}
