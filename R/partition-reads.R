

# partition_reads -----------------------------------------------------------


#' Partition a SNP matrix into two haplotype groups
#'
#' @param x A \code{[read x position]} SNP matrix.
#' @param skip_gap_freq Skip a badly behaved polymorphic position if the
#' frequency of gaps within a haplotype group exceeds \code{skip_gap_freq}.
#'
#' @return A \code{HapPart} object.
#' @export
#' @examples
#' ###
partition_reads <- function(x, cl_method="ward.D", min_len = 0.8, skip_gap_freq = 2/3){

  # get SNPs
  ppos <- colnames(x)
  xm <- x[order(rownames(x)),]

  bad_ppos <- apply(xm, 2, function(x) NROW(x[x == "-"])/NROW(x) > skip_gap_freq )
  xm <- xm[,!bad_ppos]
  bad_ppos <- ppos[bad_ppos]

  ## Get the SNPs as sequences
  x_ <- apply(xm, 1, function(t) c(unlist(paste(t, collapse = ""))))

  ## Get only the fraction of reads that contain at least min_len of total SNPs
  # might be evaluated whether really necessary; Seems so. Produces way better clustering!
  # a good cutoff seems around 0.8
  # The remaining reads need to be assigned to clades as well.
  # Maybe build a consensus of the clusters and align each missing read to the best matching cluster?
  x_sub <- x_[nchar(gsub("-", "", x_)) > min_len*NCOL(x)]
  x_sub <- Biostrings::DNAStringSet(x_sub)
  # get hamming distance ## maybe try something else too
  xdist <- Biostrings::stringDist(x_sub, method="hamming")
  xmat <- as.matrix(xdist)

  ## Perform a hierarchical clustering
  hcc <-  hclust(xdist, method = cl_method)
  plot(hcc, labels=FALSE)
  ## do a dynamic cut. Need to be evaluated$
  clusts <- dynamicTreeCut::cutreeHybrid(hcc, distM = xmat, deepSplit = 1)
  # extract the clusters for each sequence
  clades <- as.factor(unlist(lapply(unname(clusts$labels), function(i) rawToChar(as.raw(as.integer(i)+64)))))
  names(clades) <- hcc$labels

  # try cluster quality measure
  # get distance of each read to all clades
  clades <- clust_quality(clades, xmat)
  clades <- arrange(clades, read)

  # Create the partition table
  part_ <- HapPart(read_name = clades$read, snp_pos = colnames(x))
  PRT(part_) <- clades$clade
  HTR(part_) <- hcc
  mcoef(part_) <- clades$coeff
  return(part_)
}

# Helpers -----------------------------------------------------------------

# Get quality measure of all vs all clusters
clust_quality <- function(clades, xmat){
  # Get means for all clusters from distance matrix
  clade_sums <- data.frame(clade = clades)
  for (clade in levels(clades)){
    ownclade <- names(clades[clades == clade])
    sums_ownclade <- unlist(lapply(names(clades), function(t) mean(xmat[t,][names(xmat[t,]) %in% ownclade])))
    clade_sums <- dplyr::bind_cols(clade_sums, !!clade := sums_ownclade)
  }
  clade_sums <- dplyr::bind_cols(read = names(clades), clade_sums)

  ## Calculate coefficients for each vs each cluster for each read
  scored_clade_sums <- data.frame()
  for (cld in levels(clade_sums$clade)){
    clds <- levels(clade_sums$clade)
    tmpcld <- dplyr::filter(clade_sums, clade == cld)

    # The main coefficient is the distance of own cluster to that of all others
    coeff <- data.frame(tmpcld[cld]/rowMeans(tmpcld[clds[!clds == cld]]))
    names(coeff) <- "coeff"
    tmpcld <- dplyr::bind_cols(tmpcld, coeff)
    for (ocld in clds){
      cnames <- c(names(tmpcld), paste0("coeff", ocld))
      ## The finer coefficients are the distance of own cluster to those of all others separately
      indCoef <- tmpcld[cld]/tmpcld[ocld]
      tmpcld <- dplyr::bind_cols(tmpcld, replace(indCoef, indCoef == 1, 0))
      names(tmpcld) <- cnames
    }
    scored_clade_sums <- dplyr::bind_rows(scored_clade_sums, tmpcld)
  }
  scored_clade_sums
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
  class(rs) <- c("HapPart", "character")
  rs
}

## sk: replaced by attr; mcoef is membership coefficient of one over all other clusters
#' Cluster membership coefficient q
#'
#' q = Q/K: cluster weight/number of polymorphic positions
#'
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
