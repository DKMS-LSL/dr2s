# partition_reads -----------------------------------------------------------

#' Cluster a SNP matrix into haplotype groups
#'
#' @param x A \code{[read x position]} SNP matrix.
#' @param skip_gap_freq Skip a badly behaved polymorphic position if the
#' frequency of gaps within a haplotype group exceeds \code{skip_gap_freq}.
#' @param min_len Minimal fraction of SNPs that need to be covered by a read
#' to be used
#' @param threshold  Only gaps above threshold will be treated as gaps.
#' Removes noisy positions in longreads.
#' @param cl_method clustering method passed to \code{\link[stats]{hclust}}.
#' @param deepSplit sensitivity parameter passed to
#' \code{\link[dynamicTreeCut]{cutreeHybrid}}.
#'
#' @return A \code{HapPart} object.
#' @export
#' @examples
#' ###

# debug
#x <- dpb1_3$partition$mat
#sort_by <- "distance"
# cl_method="ward.D"
# min_len = 0.5
# skip_gap_freq = 2/3
# deepSplit = 1
# x <- dpb1_3$partition$mat
# x <- mat

partition_reads <- function(x, cl_method = "ward.D", min_len = 0.5,
                            skip_gap_freq = 2/3, deepSplit = 1,
                            threshold = 0.2, dist_alleles = 2, sort_by = "count"){
  match.arg(sort_by, c("count", "distance"))
  # get SNPs
  ppos <- colnames(x)
  bad_ppos <- c()
  xm <- as.matrix(x[order(rownames(x)),])
  ## if there is only one SNP for clustering, use it! If it does not match both
  ## sequencing types it will be reported
  if (length(ppos) > 1) {
    bad_ppos <- apply(xm, 2, function(x) {
      NROW(x[x == "+"])/NROW(x) > skip_gap_freq
    })
    xm <- as.matrix(xm[,!bad_ppos])
    bad_ppos <- ppos[bad_ppos]
    flog.info(paste0("  %s SNPs are covered by less than %g%% of sequences and",
                     " discarded. Using the remaining %s SNPs for clustering ..."),
              length(bad_ppos), 1 - skip_gap_freq, ncol(xm), name = "info")
    if (NCOL(xm) == 0) {
      flog.error(paste0("  Aborting. No SNP remaining for clustering!",
                        " Check your reads and reference and have a look",
                        " at mapInit plots!"))
      stop("No SNPs remaining for clustering. Check mapInit plots!")
    }
  }


  ## Get the SNP matrix as sequences
  xseqs <- get_seqs_from_mat(xm)
  ## if only one SNP, remove every read without it
  if (NCOL(xm) == 1)
    xseqs <- xseqs[!xseqs == "+"]

  ## Get only the fraction of reads that contain at least min_len of total SNPs
  clustres <- get_clusts(xseqs, cl_method = cl_method, min_len = min_len,
                         deepSplit = deepSplit, threshold = threshold)
  subclades <- factor(clustres$clades[!clustres$clades == "@"])
  tree <- clustres$tree
  hptypes <- levels(subclades)

  if (length(subclades) == 0){
    flog.error("  Two few longreads for clustering", name = "info")
    stop("To few longreads for clustering")
  }

  flog.info("  Initial clustering results in %s haplotypes %s",
            length(hptypes), comma(hptypes), name = "info")
  # Get scores and assign clades by freq mat
  ## Position Weight Matrix: Use frequency plus pseudocount/ basefrequency
  ## (here 0.25 for each).
  msa <- lapply(levels(subclades), function(x) {
    xseqs[names(subclades[subclades == x])]
  })
  names(msa) <- hptypes
  mats <- lapply(msa, function(x) create_PWM(x))
  hpseqs <- Biostrings::DNAStringSet(sapply(msa, function(x) {
    unlist(simple_consensus(
      t(Biostrings::consensusMatrix(x)[c(VALID_DNA(), "+"),])))
  }))
  names(hpseqs) <- sapply(hptypes, function(hptype)
    paste(hptype, sum(subclades == hptype), sep = ":"))

  if (length(hptypes) > dist_alleles) {
    flog.info("  Trying to identify chimeric reads/haplotypes ...",
              name = "info")
    if (sort_by == "count") {
      rC <- names(sort(table(subclades),decreasing = TRUE)[1:dist_alleles])
      # ## !!!!! HACK!!! rm following line; uncomment previous one; removed
      # rC <- c("A", "C")
    } else if (sort_by == "distance"){
      rC <- sort(find_chimeric(seqs = hpseqs, dist_alleles = dist_alleles))
    }
    flog.info("  Use only clusters %s ...", comma(rC), name = "info")
    mats <- mats[rC]
  }
  hptypes <- names(mats)

  scores <- dplyr::bind_rows(
    lapply(1:length(xseqs),function(s) get_scores(s, xseqs, mats)))

  s <- 1
    lapply(1:length(xseqs),function(s) get_scores(s, xseqs, mats))

  clades <- scores %>%
    dplyr::group_by(read) %>%
    dplyr::slice(which.max(score)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(clustclade = ifelse(read %in% names(subclades),
                                      as.character(subclades[read]),
                                      NA)) %>%
    dplyr::mutate(correct = dplyr::if_else(clustclade == clade, TRUE, FALSE))

  ## Correctly classified clades in the initial clustering
  falseClassified <- NROW(dplyr::filter(clades, correct == FALSE)) /
    NROW(dplyr::filter(clades, correct == TRUE))
  flog.info("  Corrected classification of %.2f%% of reads",
            100*falseClassified, name = "info")
  foreach(hp = hptypes) %do% {
    flog.info("  %s reads in haplotype %s", table(clades$clade)[hp],
              hp, name = "info")
  }

  # Create the partition table
  part_ <- HapPart(read_name = clades$read, snp_pos = colnames(x))
  PRT(part_) <- clades$clade
  HTR(part_) <- tree
  SNP(part_) <- ppos
  K(part_)   <- length(ppos) - length(bad_ppos)
  oc(part_)  <- as.character(unique(subclades))
  mcoef(part_)  <- clades$score
  scores(part_) <- scores
  SQS(part_)    <- hpseqs
  PWM(part_)    <- mats

  return(part_)
}


# Helpers -----------------------------------------------------------------

get_clusts <- function(xseqs, xmat, min_len = 0.80, cl_method = "ward.D",
                       deepSplit = 1, min_reads_frac = 1/3, threshold = 0.2) {

  assertthat::assert_that(
    is.double(min_len),
    is.double(min_reads_frac),
    is.double(threshold)
  )

  # heuristic for how many reads should be used
  min_len <- .getMinLenClust(min_reads_frac, min_len, xseqs)
  flog.info("  Using only longreads containing at least %s%% of all SNPs ...",
            min_len*100, name = "info")
  x_sub <- xseqs[Biostrings::width(gsub("\\+", "", xseqs)) >
                   min_len*Biostrings::width(xseqs[1])]

  ## Consensus matrix with pseudocount
  flog.info(paste0("  Constructing a Position Specific Distance Matrix",
                   " of the reamaining %s sequences ..."),
            length(x_sub), name = "info")
  consmat  <- as.matrix(
    Biostrings::consensusMatrix(x_sub, as.prob = TRUE)[c(VALID_DNA(), "+" ), ] + 1/length(x_sub)
  )
  # Remove gaps below a threshold as they are probably sequencing artifacts
  consmat <- foreach(col = 1:ncol(consmat), .combine = cbind) %do% {
    a <- consmat[,col]
    a["-"] <- ifelse(a["-"] < threshold, min(a), a["-"])
    a
  }
  ## Get Position Specific Distance Matrix
  dist <- PSDM(x_sub, as.matrix(consmat))
  dist <- as.dist(dist)

  ## replace "na" with mean for a being able to cluster
  dist[is.na(dist)] <- mean(dist,na.rm = TRUE)

  ## Perform a hierarchical clustering
  hcc <-  hclust(dist, method = cl_method)
  ## do a dynamic cut. Need to be evaluated
  clusts <- dynamicTreeCut::cutreeHybrid(hcc,
                                         minClusterSize = 15,
                                         distM = as.matrix(dist),
                                         deepSplit = deepSplit,
                                         verbose = FALSE)
  # extract the clusters for each sequence
  clades <- as.factor(sapply(unname(clusts$labels), function(i) {
    rawToChar(as.raw(as.integer(i) + 64))
  }))
  clades
  names(clades) <- hcc$labels
  return(list(clades = clades, tree = hcc))
}

## Get minimal fraction of SNPs needed for a read to be used
.getMinLenClust <- function(min_reads_frac, min_len, xseqs){
  adjusted_min_reads_frac <- min_reads_frac +
    (0.03*Biostrings::width(xseqs[1])/50) * 1500^2/length(xseqs)^2
  adjusted_min_reads_frac <- min(0.5, adjusted_min_reads_frac)
  flog.info(paste0("  Adjusting minimal fraction of reads used for clustering",
                   " from %g%% to %g%% ..."), min_reads_frac*100,
            adjusted_min_reads_frac*100, name = "info")
  # get long reads containing sufficient number of SNPs
  min_lens <- seq(1, 0, -0.01)
  len_counts <- sapply(min_lens, function(minl) {
    length(xseqs[
      Biostrings::width(gsub("\\+", "", xseqs)) >
        minl*Biostrings::width(xseqs[1])
      ])/length(xseqs)
  })
  names(len_counts) <- min_lens
  min_len <- max(as.numeric(names(
    len_counts[which(len_counts > adjusted_min_reads_frac)][1])), min_len)
}

get_scores <- function(s, xseqs, mats){
  seq <- as.character(xseqs[[s]])
  seq <- unlist(strsplit(seq, split = ""))
  read <- names(xseqs[s])

  b <- sapply(mats, function(t) {
    sum( sapply(1:length(seq), function(x) t[seq[x],x]))
  })
  t <- data.frame(read, b, names(b))
  t$read <- as.character(t$read)
  names(t) <- c("read", "score", "clade")
  t
}

## Identify chimeric clusters
find_chimeric <- function(seqs, dist_alleles, plot_seqs = FALSE) {
  # Use only non empty seqs
  seqs <- seqs[sapply(seqs, function(x) {
    !nchar(gsub("\\+|-","", as.character(x))) == 0})]

  forward <- lapply(1:(Biostrings::width(seqs[1])), function(x) {
    hammingDist(Biostrings::DNAStringSet(seqs, start = 1, width = x))})
  reverse <- lapply(1:(Biostrings::width(seqs[1])), function(x) {
    hammingDist(Biostrings::DNAStringSet(Biostrings::reverse(seqs),
                                         start = 1, width = x))})

  dm <- foreach(h = names(seqs), .combine = 'cbind') %:%
    foreach(hp = names(seqs)) %do% {
      f <- sapply(forward, function(x) x[h,hp])
      r <- sapply(reverse, function(x) x[h,hp])
      # debug
      if (plot_seqs) {
        ggplot() +
          geom_line(aes(x = 1:length(r),y = f + r)) #+
        geom_line(aes(x = length(r):1, y = r))
      }
      unname(quantile(f + r)[2])
    }
  rownames(dm) <- sapply(names(seqs), function(x) strsplit1(x, ":")[1])
  colnames(dm) <- sapply(names(seqs), function(x) strsplit1(x, ":")[1])

  c("A", (names(sort(unlist(dm[1,]), decreasing = TRUE)[1:(dist_alleles - 1)])))
}

# Class: HapPart -----------------------------------------------------------


HapPart <- function(read_name, snp_pos) {
  read_name <- as.character(read_name)
  snp_pos   <- as.integer(snp_pos)
  n <- length(read_name) # number of reads
  k <- length(snp_pos)   # number of polymorphic positions
  structure(
    read_name,
    snp_pos = snp_pos,
    k       = 0L,            # total number of polymorphic positions
    q       = rep(0, n),     # cluster weight
    # add mcoef and tree
    mcoef   = rep(0, n),     # coefficient of membership to one cluster vs all other clusters
    tree    = NULL,          # Add tree from hclust
    scores   = NULL,
    mats    = NULL,
    classification = matrix( # running classification of polymorphic positions
      rep("", n * k),
      nrow = n,
      ncol = k,
      dimnames = list(NULL, snp_pos)
    ),
    partition = rep("", n),  # final partitioning of reads into haplotypes
    class = c("HapPart", "character")
  )
}


# Methods: HapPart --------------------------------------------------------


#' @export
print.HapPart <- function(x, sort_by = "none", nrows = 8, ...) {
  sort_by <- match.arg(sort_by, c("none", "name", "mcoef"))
  n <- maximum(length(x), nrows)
  df0 <- data.frame(
    read      = paste0(substr(as.vector(x), 1, 12), "..."),
    mcoef     = mcoef(x),
    partition = PRT(x)
  )[1:n, ]
  cat("HapPart over", K(x), "polymorphic positions:\n")
  switch(
    sort_by,
    none  = print(df0, ...),
    name  = print(df0[order(read),], ...),
    mcoef = print(df0[order(mcoef),], ...)
  )
  invisible(x)
}

#' @export
`[.HapPart` <- function(x, i, j = NULL, drop = NULL) {
  i <- if (is.character(i)) x %in% i else i
  rs <- NextMethod()
  attr(rs, "snp_pos") <- SNP(x)
  attr(rs, "k") <- K(x)
  attr(rs, "q") <- Q(x)[i]
  # add mcoef attr
  attr(rs, "mcoef") <- mcoef(x)[i]
  attr(rs, "classification") <- CLS(x)[i, , drop = FALSE]
  attr(rs, "partition") <- PRT(x)[i]
  attr(rs, "tree") <- HTR(x)
  attr(rs, "scores") <- scores(x)
  class(rs) <- c("HapPart", "character")
  rs
}

#' @export
mcoef <- function(x) UseMethod("mcoef")
#' @export
mcoef.HapPart <- function(x) {
  attr(x, "mcoef")
}
`mcoef<-` <- function(x, value) UseMethod("mcoef<-")
#' @export
`mcoef<-.HapPart` <- function(x, value) {
  attr(x, "mcoef") <- value
  x
}

#' The PWM matrix of a cluster
#' @export
PWM <- function(x) UseMethod("PWM")
#' @export
PWM.HapPart <- function(x) {
  attr(x, "PWM")
}
`PWM<-` <- function(x, value) UseMethod("PWM<-")
#' @export
`PWM<-.HapPart` <- function(x, value) {
  attr(x, "PWM") <- value
  x
}

#' Get the original number of clusters for plotting the tree
#' @export
oc <- function(x) UseMethod("oc")
#' @export
oc.HapPart <- function(x) {
  attr(x, "oc")
}
`oc<-` <- function(x, value) UseMethod("oc<-")
#' @export
`oc<-.HapPart` <- function(x, value){
  attr(x, "oc") <- value
  x
}

#' All vs all cluster distances for each read
#' @export
scores <- function(x) UseMethod("scores")
#' @export
scores.HapPart <- function(x) {
  attr(x, "scores")
}
`scores<-` <- function(x, value) UseMethod("scores<-")
#' @export
`scores<-.HapPart` <- function(x, value) {
  attr(x, "scores") <- value
  x
}

#' @export
partition <- function(x) UseMethod("partition")
#' @export
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
K.HapPart <- function(x) {
  attr(x, "k")
}

`K<-` <- function(x, value) UseMethod("K<-")
`K<-.HapPart` <- function(x, value) {
  attr(x, "k") <- value
  x
}

SQS <- function(x) UseMethod("SQS")
SQS.HapPart <- function(x) {
  attr(x, "seqs")
}

`SQS<-` <- function(x, value) UseMethod("SQS<-")
`SQS<-.HapPart` <- function(x, value) {
  attr(x, "seqs") <- value
  x
}

CLS <- function(x) UseMethod("CLS")
CLS.HapPart <- function(x) {
  attr(x, "classification")
}

`CLS<-` <- function(x, value) UseMethod("CLS<-")
`CLS<-.HapPart` <- function(x, value) {
  attr(x, "classification")[, K(x)] <- value
  x
}

SNP <- function(x) UseMethod("SNP")
SNP.HapPart <- function(x) {
  attr(x, "snp_pos")
}

`SNP<-` <- function(x, value) UseMethod("SNP<-")
`SNP<-.HapPart` <- function(x, value) {
  attr(x, "snp_pos") <- value
  x
}

PRT <- function(x) UseMethod("PRT")
PRT.HapPart <- function(x) {
  attr(x, "partition")
}

`PRT<-` <- function(x, value) UseMethod("PRT<-")
`PRT<-.HapPart` <- function(x, value) {
  attr(x, "partition") <- value
  x
}

# get tree from hclust
HTR <- function(x) UseMethod("HTR")
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
plot_partition_histogram <- function(x, label = "", limits = NULL) {
  stopifnot(is(x, "HapPart"))
  if (is.null(limits)){
    limits <- self$getLimits()
  }
  data <- partition(x) %>%
    dplyr::mutate(mcoef = ifelse(haplotype == "A", mcoef, -1*mcoef))

  colors <- c(x)
  ggplot(data) +
    geom_histogram(aes(x = mcoef, fill = haplotype), bins = 100) +
    scale_fill_manual(values = PARTCOL()) +
    geom_vline(xintercept = c(limits[1], -limits[2]), linetype = "dashed",
               colour = "grey80") +
    xlab("Haplotype membership coefficient") +
    ylab("Number of reads") +
    ggtitle(label = label) +
    theme_bw()
}

#' Plot distribution of haplotype partitioned reads for more than two haplotypes
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
plot_partition_histogram_multi <- function(x, limits, label = "") {
  stopifnot(is(x, "HapPart"))
  data <- partition(x) %>%
    dplyr::mutate(limit = unlist(limits)[haplotype])

  ggplot(data) +
    geom_histogram(aes(x = mcoef, fill = haplotype), bins = 100) +
    facet_grid(~ haplotype) +
    # geom_vline(xintercept = c(0.25, 0.5, 0.75),
    # linetype = "dashed", size = 0.2, colour = "grey80") +
    geom_vline( aes(xintercept = limit), colour = "grey40",
                linetype = "dashed", size = 1) +
    scale_fill_manual(values = PARTCOL()) +
    # xlim(c(0, 1.1)) +
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
#' @param name_reads Name reads [FALSE]
#'
#' @return A \code{ggplot} object.
#' @export
#' @examples
#' ###
plot_partition_haplotypes <- function(x, thin = 1, label = "", sort = TRUE,
                                      name_reads = FALSE) {
  stopifnot(is(x, "HapPart"))
  q <- mcoef(x)
  snp_pos <- SNP(x)
  rs <- CLS(x)[, order(snp_pos), drop = FALSE]
  i <- if (thin < 1) {
    sample(NROW(rs), thin*NROW(rs))
  } else seq_len(NROW(rs))
  q_order <- if (sort) {
    rev(order(q[i], decreasing = TRUE))
  } else i
  rs2 <- rs <- rs[i, , drop = FALSE][q_order, , drop = FALSE]
  dim(rs2) <- NULL
  alpha_steps <- function(a) {
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
    trans = alpha_steps(a = rep(q[i][q_order], NCOL(rs))),
    mcoef = rep(q[i][q_order], NCOL(rs))
  )

  if (name_reads) {
    readnames <- as.vector(x[i])[q_order]
    readnames <- factor(readnames, readnames, readnames, ordered = TRUE)
    df$read   <- rep(readnames, NCOL(rs))
  }

  ggplot(df) +
    geom_tile(aes(x = snp, y = read, fill = hap, alpha = trans)) +
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
plot_radar_partition <- function(x){
  stopifnot(is(x, "HapPart"))

  df <- dplyr::full_join(scores(x), partition(x), by = "read")

  # use bigger size for 2d radarplot.
  size <- ifelse(length(levels(df$haplotype)) == 2, 4, 0.05)
  ggplot(df, aes(x = clade, y = score)) +
    geom_polygon(aes(group = read, color = haplotype), fill = NA, size = size,
                 show.legend = FALSE, alpha = 1) +
    facet_grid( ~ haplotype) +
    # scale_fill_manual(PARTCOL()) +
    ggtitle("Similarity to clusters") +
    scale_color_manual(values = PARTCOL()) +
    theme_bw() +
    theme(
      axis.title = element_blank(),
      axis.text.y  = element_blank(),
      axis.ticks  = element_blank()
    ) +
    coord_radar()
}

#' Plot tree of initial clustering
#'
#' @param x A \code{HapPart} object.
#'
#' @return A \code{ggplot} object.
#' @export
#' @examples
#' ###
#x <- self$getPartition()
plot_partition_tree <- function(x){
  stopifnot(
    is(x, "HapPart"),
    requireNamespace("ggdendro", quietly = TRUE)
  )
  tree <- HTR(x)
  k <- length(oc(x))

  tryCatch({
    dendr <- ggdendro::dendro_data(tree, type = "rectangle")
    dendr$labels <- dendr$labels %>%
      dplyr::mutate(label = as.character(label))

    clust <- cutree(tree, k)
    clust <- dplyr::data_frame(label = names(clust), haplotype = clust)
    dendr$labels <- dplyr::left_join(dendr$labels, clust, by = "label")

    height <- unique(dendr$segments$y)[order(unique(dendr$segments$y),
                                             decreasing = TRUE)]
    cut.height <- mean(c(height[k], height[k - 1]))
    dendr$segments <- dendr$segments %>%
      dplyr::mutate(line =  dplyr::if_else(
        y == yend & y > cut.height, 1, dplyr::if_else(
          yend > cut.height,1, 2)))

    dendr$segments$cluster <- c(-1, diff(dendr$segments$line))
    change <- which(dendr$segments$cluster == 1)
    for (i in 1:k) dendr$segments$cluster[change[i]] = i + 1
    dendr$segments <- dendr$segments %>%
      dplyr::mutate(cluster = dplyr::if_else( line == 1, 1, ifelse(
        cluster == 0, NA, cluster)))
    dendr$segments$cluster <- sapply(1:NROW(dendr$segments$cluster), function(x) {
      getCl(x, dendr$segments$cluster, change)})

    # Correct order
    labs <- c("N", oc(x))
    dendr$segments$cluster <- labs[dendr$segments$cluster]

    # Make plot
    p <- ggplot() +
      geom_segment(data=ggdendro::segment(dendr), aes(x = x, y = y,
                                                      xend = xend, yend = yend,
                                                      color = factor(cluster)),
                   size = 1.25) +
      scale_color_manual(values = PARTCOL(),
                         name = "Haplotype",
                         breaks = LETTERS[1:k],
                         labels =  LETTERS[1:k]) +
      geom_hline(yintercept = cut.height, color = "blue") +
      ggtitle("Initial clustering of haplotypes") +
      theme_bw() +
      theme(legend.position = "bottom",
            axis.title = element_blank(),
            axis.text  = element_blank(),
            axis.ticks  = element_blank()
      )
    p
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
coord_radar <- function(theta = "x", start = 0, direction = 1)  {
  theta <- match.arg(theta, c("x", "y"))
  r <- if (theta == "x")
    "y"
  else "x"
  ggproto("CoordRadar", CoordPolar, theta = theta, r = r, start = start,
          direction = sign(direction),
          is_linear = function(coord) TRUE)
}

# Optimal partitioning ----------------------------------------------------

optimal_partition_limits <- function(scores, f = 0.8) {
  coeff_cuts <- seq(0, max(unlist(scores)), length.out = 100)
  rHap <- lapply(scores, function(x) {
    sapply(coeff_cuts, function(cutoff) sum(x >= cutoff))})
  df <- data.frame()
  for(hapType in names(rHap)) {
    l <- f * length(unlist(
      scores[!names(rHap) == hapType])) / length(scores[[hapType]])
    df <- dplyr::bind_rows(df,
                           dplyr::data_frame(c = coeff_cuts,
                                             r = unlist(rHap[hapType])) %>%
                             dplyr::mutate(haplotype = hapType,
                                           score = ((r^l)*c)) %>%
                             dplyr::mutate(max.score = score == max(score)))
  }
  dfmax <- dplyr::filter(df, max.score) %>%
    dplyr::select(haplotype, c, r, score)

  list(
    limits = dfmax,
    plt = ggplot(df, aes(r, c, colour = haplotype)) +
      geom_line() +
      geom_point(data = dfmax, size = 3) +
      ylab("Cluster similarity coefficient") +
      xlab("N reads") +
      theme_bw()
  )
}
