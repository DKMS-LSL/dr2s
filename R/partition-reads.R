

# partition_reads -----------------------------------------------------------


#' Partition a SNP matrix into two haplotype groups
#'
#' @param x A \code{[read x position]} SNP matrix.
#' @param shuffle Shuffle positions.
#' @param skip_gap_freq Skip a badly behaved polymorphic position if the
#' frequency of gaps within a haplotype group exceeds \code{skip_gap_freq}.
#'
#' @return A \code{HapPart} object.
#' @export
#' @examples
#' ###
partition_reads <- function(x, shuffle = TRUE, skip_gap_freq = 2/3) {
  x_ <- if (shuffle && NCOL(x) > 1) {
    x[, sample(NCOL(x))]
  } else x

  part_ <- HapPart(read_name = rownames(x_), snp_pos = colnames(x_))
  read_names_  <- as.vector(part_)
  rownames(x_) <- NULL

  slice <- x_[, 1L]
  nm <- top2(slice)

  isA <- slice == nm[1L]
  isB <- slice == nm[2L]

  read_names_0A <- read_names_[isA]
  read_names_0B <- read_names_[isB]

  K(part_) <- K(part_) + 1L

  cls_ <- ifelse(isA, "A", ifelse(isB, "B", "-"))

  Q(part_)[read_names_ %in% read_names_0A] <- Q(part_)[read_names_ %in% read_names_0A] + 1L
  Q(part_)[read_names_ %in% read_names_0B] <- Q(part_)[read_names_ %in% read_names_0B] - 1L

  CLS(part_) <- cls_
  PRT(part_) <- ifelse(cls_ == "-", sample(c("A", "B"), 1), cls_)
  partition_reads_(x_ = x_[, -1, drop = FALSE], part_, read_names_, skip_gap_freq)
}

# Helpers -----------------------------------------------------------------

top2 <- function(slice) {
  lv <- levels(slice)
  ilv <- which(lv != ".")
  tb <- tabulate(slice, nbins = length(lv))[ilv]
  names(tb) <- lv[ilv]
  names(sort(tb, decreasing = TRUE)[1:2])
}

partition_reads_ <- function(x_, part_, read_names_, skip_gap_freq) {
  if (NCOL(x_) == 0) {
    return(part_)
  }

  slice <- x_[, 1L]

  # extract the partitioning scheme from the previous iteration and
  # assign the next polymorphic position according to this scheme
  if (length(nm <- assign_partition(slice, part = PRT(part_), skip_gap_freq)) == 2) {
    isA <- slice == nm[["A"]]
    isB <- slice == nm[["B"]]

    read_names_1A <- read_names_[isA]
    read_names_1B <- read_names_[isB]

    K(part_) <- K(part_) + 1L

    cls_ <- ifelse(isA, "A", ifelse(isB, "B", "-"))

    Q(part_)[read_names_ %in% read_names_1A] <- Q(part_)[read_names_ %in% read_names_1A] + 1L
    Q(part_)[read_names_ %in% read_names_1B] <- Q(part_)[read_names_ %in% read_names_1B] - 1L

    # Classify
    CLS(part_) <- cls_

    # Partition
    PRT(part_) <- ifelse(part_ %in% read_names_[Q(part_) >= 0L], "A", "B")
  } else {
    ## excise offending SNP
    SNP(part_) <- SNP(part_)[-K(part_) -  1L]
    ## excise offending column in classification table
    attr(part_, "classification") <- attr(part_, "classification")[, -K(part_) - 1L]
  }

  #x_ <- x_[, -1, drop = FALSE]
  #(slice <- x_[, 1L])
  Recall(x_[, -1, drop = FALSE], part_, read_names_, skip_gap_freq)
}

assign_partition <- function(slice, part, skip_gap_freq) {
  lvls_  <- levels(slice)
  nbins_ <- length(lvls_)
  sliceA <- slice[part == "A"]
  sliceB <- slice[part == "B"]

  # Tabulate the two parts and calculate the fraction of
  # nucleotides with each part
  A <- array(tabulate(sliceA, nbins = nbins_), dimnames = list(lvls_))
  B <- array(tabulate(sliceB, nbins = nbins_), dimnames = list(lvls_))

  # check frequency of non-nucleotides
  if (A[names(A) == "."]/sum(A) > skip_gap_freq || B[names(B) == "."]/sum(B) > skip_gap_freq) {
    message("Skipping badly behaved polymorphic position")
    return(NULL)
  }

  # remove non-nucleotide
  A <- A[names(A) != "."]
  B <- B[names(B) != "."]

  A <- A/sum(A)
  B <- B/sum(B)

  # An allele will be overrepresented in it's correct partition and under-
  # represented in the opposite partition. Spurious alleles will be randomly
  # distributed across partitions. The two correctly partitioned alleles will
  # show the maximum frequency difference across the partiotions.
  d <- (A - B)^2
  nm <- names(sort(d, decreasing = TRUE)[1:2])

  ## We calculate the determinant of the square matrix formed by
  ## the frequencies of alleles "A" and "B" in parts A and be respectively.
  mx <- rbind(A = A[nm], B = B[nm])

  if (det2(mx) > 0) {
    Anm <- nm[1L]
    Bnm <- nm[2L]
  } else {
    Anm <- nm[2L]
    Bnm <- nm[1L]
  }

  c(A = Anm, B = Bnm)
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
  attr(rs, "classification") <- CLS(x)[i, , drop = FALSE]
  attr(rs, "partition") <- PRT(x)[i]
  class(rs) <- c("HapPart", "character")
  rs
}

#' Cluster membership coefficient q
#'
#' q = Q/K: cluster weight/number of polymorphic positions
#'
#' @export
mcoef <- function(x) UseMethod("mcoef")
#' @export
mcoef.HapPart <- function(x) {
  Q(x)/K(x)
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

## Cluster weight
##
## Number of polymorphic positions that support membership in
## in one over the other haplotype
Q <- function(x) UseMethod("Q")
Q.HapPart <- function(x) {
  attr(x, "q")
}

`Q<-` <- function(x, value) UseMethod("Q<-")
`Q<-.HapPart` <- function(x, value) {
  attr(x, "q") <- value
  x
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

# Summarise and Plot ------------------------------------------------------

#' Plot distribution of haplotype partitioned reads
#'
#' @param x A \code{HapPart} object.
#' @param label Optional plot label.
#'
#' @return A \code{ggplot} object.
#' @export
#' @examples
#' ###
plot_partition_histogram <- function(x, label = "") {
  stopifnot(is(x, "HapPart"))
  ggplot(partition(x)) +
    geom_histogram(aes(x = mcoef, fill = haplotype), binwidth = 0.05) +
    scale_fill_manual(values = PARTCOL()) +
    geom_vline(xintercept = 0L, linetype = "dashed", colour = "grey80") +
    geom_vline(xintercept = c(-0.75, -0.5, -0.25, 0.25, 0.5, 0.75),
               linetype = "dashed", size = 0.5, colour = "grey80") +
    xlim(c(-1.1, 1.1)) +
    xlab("Haplotype membership coefficient") +
    ylab("Number of reads") +
    ggtitle(label = label) +
    theme_bw()
}

#' Tile plot of haplotyped reads.
#'
#' @param x A \code{HapPart} object.
#' @param thin Subsample reads.
#' @param label Optional plot label.
#' @param sort Sort
#' @param name_reads Name reads
#'
#' @return A \code{ggplot} object.
#' @export
#' @examples
#' ###
plot_partition_haplotypes <- function(x, thin = 1, label = "", sort = TRUE, name_reads = FALSE) {
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
    factor(a, levels = 0:3, labels = c("q <= 0.25", "0.25 < q <= 0.50", "0.50 < q <= 0.75", "0.75 < q <= 1"))
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
    labs(x = "Polymorphic positions", y = "Reads", fill = "Haplotype", alpha = "Coefficient") +
    ggtitle(label = label) +
    theme_bw() +
    theme(legend.position = "top")
}


# Optimal partitioning ----------------------------------------------------


optimal_partition_limits <- function(cA, cB, scale = 2, correction = TRUE) {
  coeff_cuts <- seq(0.95, 0.00, by = -0.05)
  scaleA <- scaleB <- scale
  rA <- sapply(coeff_cuts, function(cutoff) sum(cA >= cutoff))
  rB <- sapply(coeff_cuts, function(cutoff) sum(cB >= cutoff))
  if (correction) {
    mAB <- c(max(rA), max(rB))
    i <- which.min(mAB)
    corrfact <- sqrt(2 * (mAB[i]/sum(mAB)))
    scaleA <- if (i == 1) scaleA * corrfact else scaleA
    scaleB <- if (i == 2) scaleB * corrfact else scaleB
  }
  df <- dplyr::bind_rows(
    dplyr::data_frame(c = coeff_cuts, r = rA) %>%
      dplyr::mutate(haplotype = "A", score = (c^scaleA * r), scale = scaleA) %>%
      dplyr::mutate(max.score = score == max(score)),
    dplyr::data_frame(c = coeff_cuts, r = rB) %>%
      dplyr::mutate(haplotype = "B", score = (c^scaleB * r), scale = scaleB) %>%
      dplyr::mutate(max.score = score == max(score))
  )
  dfmax <- dplyr::filter(df, max.score) %>%
    dplyr::select(haplotype, c, r, score, scale)
  list(
    limits = dfmax,
    plt = ggplot(df, aes(r, c, colour = haplotype)) +
      geom_line() +
      geom_point(data = dfmax, size = 3) +
      ylab("Membership coefficient") +
      xlab("N reads") +
      theme_bw()
  )
}
