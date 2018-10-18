# partitionReads -----------------------------------------------------------

#' Cluster a SNP matrix into haplotype groups
#'
#' @param x A \code{[read x position]} SNP matrix.
#' @param distAlleles Number of distinct alleles in the sample.
#' @param sortBy sort by "distance" or "count". See
#' \code{\link{partitionLongReads}} for details
#' @param threshold  Only gaps above threshold will be treated as gaps.
#' Removes noisy positions in longreads.
#' @param clMethod clustering method passed to \code{\link[stats]{hclust}}.
#' @param minLen Minimal fraction of SNPs that need to be covered by a read
#' to be used
#' @param skipGapFreq Skip a badly behaved polymorphic position if the
#' frequency of gaps within a haplotype group exceeds \code{skipGapFreq}.
#' @param deepSplit sensitivity parameter passed to
#' \code{\link[dynamicTreeCut]{cutreeHybrid}}.
#'
#' @return A \code{HapPart} object.
#' @export
#' @examples
#' ###

partitionReads <- function(x, distAlleles = 2, sortBy = "count", threshold = 0.2,
                           clMethod = "ward.D", minLen = 0.5, skipGapFreq = 2/3,
                           deepSplit = 1, minClusterSize = 15) {
  sortBy <- match.arg(sortBy, c("count", "distance"))
  # get SNPs
  ppos <- colnames(x)
  badPpos <- c()
  xm <- as.matrix(x[order(rownames(x)), , drop = FALSE])
  ## if there is only one SNP for clustering, use it! If it does not match both
  ## sequencing types it will be reported
  if (length(ppos) > 1) {
    badPpos <- apply(xm, 2, function(x) {
      NROW(x[x == "+"])/NROW(x) > skipGapFreq
    })
    xm <- as.matrix(xm[, !badPpos])
    badPpos <- ppos[badPpos]
    flog.info("%s SNPs are covered by less than %g%% of sequences and" %<<%
               " discarded. Using the remaining %s SNPs for clustering ...",
              length(badPpos), 1 - skipGapFreq, NCOL(xm), name = "info")
    if (NCOL(xm) == 0) {
      flog.error("  Aborting. No SNP remaining for clustering!" %<<%
                  " Check your reads and reference and have a look" %<<%
                        " at mapInit plots!")
      stop("No SNPs remaining for clustering. Check mapInit plots!")
    }
  }

  ## Get the SNP matrix as sequences
  xseqs <- .getSeqsFromMat(xm)
  ## if only one SNP, remove every read without it
  if (NCOL(xm) == 1)
    xseqs <- xseqs[!xseqs == "+"]

  ## Get only the fraction of reads that contain at least minLen of total SNPs
  clustres <- .getClusts(xseqs, clMethod = clMethod, minLen = minLen,
                         deepSplit = deepSplit, threshold = threshold,
                         minClusterSize = minClusterSize)
  subclades <- factor(clustres$clades[!clustres$clades == "@"])
  tree <- clustres$tree
  hptypes <- levels(subclades)

  if (length(subclades) == 0) {
    flog.error("  Two few longreads for clustering", name = "info")
    stop("To few longreads for clustering")
  }

  flog.info("  Initial clustering results in %s haplotypes %s",
            length(hptypes), comma(hptypes), name = "info")

  ## Get scores and assign clades by freq mat
  ## Position Weight Matrix: Use frequency plus pseudocount/ basefrequency
  ## (here 0.25 for each).
  msa <- lapply(levels(subclades), function(x) {
    xseqs[names(subclades[subclades == x])]
  })
  names(msa) <- hptypes
  mats <- lapply(msa, function(x) createPWM(x))
  hpseqs <- Biostrings::DNAStringSet(vapply(msa, function(x) {
    unlist(conseq(
      t(Biostrings::consensusMatrix(x)[c(VALID_DNA(), "+"), ]),
      type = "simple"))
  }, FUN.VALUE = character(1)))
  names(hpseqs) <- vapply(hptypes, function(hptype)
    paste(hptype, sum(subclades == hptype), sep = ":"),
    FUN.VALUE = character(1))

  if (length(hptypes) > distAlleles) {
    flog.info("  Trying to identify chimeric reads/haplotypes ...",
              name = "info")
    if (sortBy == "count") {
      rC <- names(sort(table(subclades),
                       decreasing = TRUE)[seq_len(distAlleles)])
    } else if (sortBy == "distance") {
      rC <- sort(.findChimeric(seqs = hpseqs, distAlleles = distAlleles))
    }
    flog.info("  Use only clusters %s ...", comma(rC), name = "info")
    mats <- mats[rC]
  }
  hptypes <- names(mats)

  scores <- dplyr::bind_rows(
    bplapply(seq_along(xseqs), function(s) .getScores(s, xseqs, mats)))

  clades <- scores %>%
    dplyr::group_by(.data$read) %>%
    dplyr::slice(which.max(.data$score)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(clustclade = ifelse(.data$read %in% names(subclades),
                                      as.character(subclades[.data$read]),
                                      NA)) %>%
    dplyr::mutate(correct = dplyr::if_else(.data$clustclade == .data$clade,
                                           TRUE, FALSE))

  ## Correctly classified clades in the initial clustering
  falseClassified <- NROW(dplyr::filter(clades, .data$correct == FALSE)) /
    NROW(dplyr::filter(clades, .data$correct == TRUE))
  flog.info("  Corrected classification of %.2f%% of reads",
            100*falseClassified, name = "info")
  lapply(hptypes, function(hp, clades) {
    invisible(flog.info("  %s reads in haplotype %s", table(clades$clade)[hp],
                        hp, name = "info"))
   }, clades = clades)

  # Create the partition table
  part_ <- HapPart(readNames = clades$read, snpPos = colnames(xm))
  PRT(part_) <- clades$clade
  HTR(part_) <- tree
  K(part_)   <- length(ppos) - length(badPpos)
  oc(part_)  <- as.character(unique(subclades))
  mcoef(part_)  <- clades$score
  scores(part_) <- scores
  SQS(part_)    <- hpseqs
  PWM(part_)    <- mats

  return(part_)
}


# Helpers -----------------------------------------------------------------

.getClusts <- function(xseqs, minLen = 0.80, clMethod = "ward.D",
                       deepSplit = 1, minReadsFrac = 1/3, threshold = 0.2,
                       minClusterSize = 15) {

  assert_that(
    is.double(minLen),
    is.double(minReadsFrac),
    is.double(threshold)
  )

  # heuristic for how many reads should be used
  minLen <- .getMinLenClust(minReadsFrac, minLen, xseqs)
  flog.info("  Using only longreads containing at least %s%% of all SNPs ...",
            minLen*100, name = "info")
  xSub <- xseqs[Biostrings::width(gsub("\\+", "", xseqs)) >
                   minLen*Biostrings::width(xseqs[1])]

  ## Consensus matrix with pseudocount
  flog.info("  Constructing a Position Specific Distance Matrix" %<<%
            " of the remaining %s sequences ...",
            length(xSub), name = "info")
  consmat  <- as.matrix(
    Biostrings::consensusMatrix(xSub, as.prob = TRUE)[c(VALID_DNA(), "+" ), ] +
      1/length(xSub)
  )
  # Remove gaps below a threshold as they are probably sequencing artifacts
  con <- apply(consmat, 2, function(a) {
    a["-"] <- ifelse(a["-"] < threshold, min(a), a["-"])
    a
  })

  ## Get Position Specific Distance Matrix
  dist <- PSDM(xSub, as.matrix(consmat))
  dist <- stats::as.dist(dist)

  ## replace "na" with mean for a being able to cluster
  dist[is.na(dist)] <- mean(dist,na.rm = TRUE)

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
.getMinLenClust <- function(minReadsFrac, minLen, xseqs){
  adjustedMinReadsFrac <- minReadsFrac +
    (0.03*Biostrings::width(xseqs[1])/50) * 1500^2/length(xseqs)^2
  adjustedMinReadsFrac <- min(0.5, adjustedMinReadsFrac)
  flog.info("  Adjusting minimal fraction of reads used for clustering" %<<%
              " from %g%% to %g%% ...", minReadsFrac*100,
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

.getScores <- function(s, xseqs, mats) {
  seq <- as.character(xseqs[[s]])
  seq <- unlist(strsplit(seq, split = ""))
  read <- names(xseqs[s])

  b <- vapply(mats, function(t) {
    sum(vapply(seq_along(seq), function(x) t[seq[x],x], FUN.VALUE = numeric(1)))
  }, numeric(1))
  t <- data.frame(read, b, names(b))
  t$read <- as.character(t$read)
  names(t) <- c("read", "score", "clade")
  t
}

## Identify chimeric clusters
.findChimeric <- function(seqs, distAlleles, plotSeqs = FALSE) {
  # Use only non empty seqs
  seqs <- seqs[vapply(seqs, function(x) {
    !nchar(gsub("\\+|-","", as.character(x))) == 0}, FUN.VALUE = logical(1))]

  forward <- lapply(seq_len(Biostrings::width(seqs[1])), function(x) {
    hammingDist(Biostrings::DNAStringSet(seqs, start = 1, width = x))})
  reverse <- lapply(seq_len(Biostrings::width(seqs[1])), function(x) {
    hammingDist(Biostrings::DNAStringSet(Biostrings::reverse(seqs),
                                         start = 1, width = x))})

  dm <- foreach(h = names(seqs), .combine = 'cbind') %:%
    foreach(hp = names(seqs)) %do% {
      f <- vapply(forward, function(x) x[h,hp], FUN.VALUE = numeric(1))
      r <- vapply(reverse, function(x) x[h,hp], FUN.VALUE = numeric(1))
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
  k <- length(snpPos)   # number of polymorphic positions
  structure(
    readNames,
    snpPos = snpPos,
    k       = 0L,            # total number of polymorphic positions
    # add mcoef and tree
    mcoef   = rep(0, n),
    tree    = NULL,          # Add tree from hclust
    scores   = NULL,
    mats    = NULL,
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
  n <- maximum(length(x), nrows)
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
  attr(rs, "scores") <- scores(x)
  class(rs) <- c("HapPart", "character")
  rs
}

mcoef <- function(x) UseMethod("mcoef")
#' @describeIn HapPart
#' Get the membership coefficient mcoef.
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
oc <- function(x) UseMethod("oc")
#' @describeIn HapPart
#' Get the original clusters from \code{\link[stats]{hclust}}.
#' @export
oc.HapPart <- function(x) {
  attr(x, "oc")
}
`oc<-` <- function(x, value) UseMethod("oc<-")
`oc<-.HapPart` <- function(x, value){
  attr(x, "oc") <- value
  x
}

## All vs all cluster distances for each read
scores <- function(x) UseMethod("scores")
#' @describeIn HapPart
#' Get the membership scores for each read for each cluster.
#' @export
scores.HapPart <- function(x) {
  attr(x, "scores")
}
`scores<-` <- function(x, value) UseMethod("scores<-")
`scores<-.HapPart` <- function(x, value) {
  attr(x, "scores") <- value
  x
}

## The actual partition table
partition <- function(x) UseMethod("partition")
partition.HapPart <- function(x) {
  dplyr::arrange(dplyr::data_frame(
    read      = as.vector(x),
    haplotype = as.factor(PRT(x)),
    mcoef     = mcoef(x)
  ),
  dplyr::desc(mcoef))
}

## Total number of polymorphic positions between haplotypes
K <- function(x) UseMethod("K")
#' @describeIn HapPart
#' Get the total number of polymorphic positions.
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

## Get the partition table
PRT <- function(x) UseMethod("PRT")
#' @describeIn HapPart
#' Get the haplotype partition as a
#' \code{\link[base]{data.frame}} containing
#' columns of \code{read}, \code{haplotype} and \code{mcoef} for each read.
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
#' Get the resulting tree from \code{\link[stats]{hclust}}.
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
  if (is.null(limits)){
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

#' Tile plot of haplotyped reads.
#'
#' @param x A \code{HapPart} object.
#' @param thin Subsample reads. [1]
#' @param label Optional plot label.
#' @param sort Sort [TRUE]
#' @param nameReads Name reads [FALSE]
#'
#' @return A \code{ggplot} object.
#' @export
#' @examples
#' ###
plotPartitionHaplotypes <- function(x, thin = 1, label = "", sort = TRUE,
                                      nameReads = FALSE) {
  stopifnot(is(x, "HapPart"))
  q <- mcoef(x)
  snpPos <- SNP(x)
  rs <- CLS(x)[, order(snpPos), drop = FALSE]
  i <- if (thin < 1) {
    sample(NROW(rs), thin*NROW(rs))
  } else seq_len(NROW(rs))
  qOrder <- if (sort) {
    rev(order(q[i], decreasing = TRUE))
  } else i
  rs2 <- rs <- rs[i, , drop = FALSE][qOrder, , drop = FALSE]
  dim(rs2) <- NULL
  alphaSteps <- function(a) {
    a <- abs(a)
    a <- ifelse(a <= 0.25, 0L, ifelse(a <= 0.5, 1L, ifelse(a <= 0.75, 2L, 3L)))
    factor(a, levels = 0:3,
           labels = c("q <= 0.25", "0.25 < q <= 0.50",
                      "0.50 < q <= 0.75", "0.75 < q <= 1"))
  }

  df <- dplyr::data_frame(
    snp   = rep(seq_len(NCOL(rs)), each = NROW(rs)),
    read  = rep.int(seq_len(NROW(rs)), NCOL(rs)),
    hap   = factor(rs2, levels = c("A", "B", "-"), ordered = TRUE),
    trans = alphaSteps(a = rep(q[i][qOrder], NCOL(rs))),
    mcoef = rep(q[i][qOrder], NCOL(rs))
  )

  if (nameReads) {
    readnames <- as.vector(x[i])[qOrder]
    readnames <- factor(readnames, readnames, readnames, ordered = TRUE)
    df$read   <- rep(readnames, NCOL(rs))
  }

  ggplot(df) +
    geom_tile(aes(x =~ snp, y =~ read, fill =~ hap, alpha =~ trans)) +
    scale_alpha_discrete(range = c(0.25, 0.95)) +
    scale_fill_manual(values = PARTCOL()) +
    labs(x = "Polymorphic positions", y = "Reads", fill = "Haplotype",
         alpha = "Coefficient") +
    ggtitle(label = label) +
    theme_bw() +
    theme(legend.position = "top")
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

  df <- dplyr::full_join(scores(x), partition(x), by = "read")

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
  k <- length(oc(x))

  tryCatch({
    dendr <- ggdendro::dendro_data(tree, type = "rectangle")
    dendr$labels <- dendr$labels %>%
      dplyr::mutate(label = as.character(.data$label))

    clust <- stats::cutree(tree, k)
    clust <- dplyr::data_frame(label = names(clust), haplotype = clust)
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
    labs <- c("N", oc(x))
    dendr$segments$cluster <- factor(labs[dendr$segments$cluster])

    # Make plot
    ggplot() +
      geom_segment(data=ggdendro::segment(dendr), aes_string(x = 'x', y = 'y',
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

.optimalPartitionLimits <- function(scores, f = 0.8) {
  coeffCuts <- seq(0, max(unlist(scores)), length.out = 100)
  rHap <- lapply(scores, function(x, coeffCuts) {
      vapply(coeffCuts, function(cutoff) sum(x >= cutoff),
             FUN.VALUE = double(1))
    },
    coeffCuts = coeffCuts)
  df <- data.frame()
  for(hapType in names(rHap)) {
    l <- f * length(unlist(
      scores[!names(rHap) == hapType])) / length(scores[[hapType]])
    df <- dplyr::bind_rows(df,
                           dplyr::data_frame(c = coeffCuts,
                                             r = unlist(rHap[hapType])) %>%
                             dplyr::mutate(haplotype = hapType,
                                           score = ((.data$r^l)*c)) %>%
                             dplyr::mutate(max.score = .data$score ==
                                             max(.data$score)))
  }
  dfmax <- dplyr::filter(df, .data$max.score) %>%
    dplyr::select(.data$haplotype, .data$c, .data$r, .data$score)
  list(
    limits = dfmax,
    plt = ggplot(df, aes_(x =~ r, y =~ c, colour =~ haplotype)) +
      geom_line() +
      geom_point(data = dfmax, size = 3) +
      ylab("Cluster similarity coefficient") +
      xlab("N reads") +
      theme_bw()
  )
}
