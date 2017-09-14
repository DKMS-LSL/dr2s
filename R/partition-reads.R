

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

## debug ##
# x <- dl3.part$partition$mat
#  cl_method="ward.D"
#  min_len = 0.8
# skip_gap_freq = 2/3
#partn <- partition_reads(x)
#length(partn[mcoef(partn)>1])
#length(partn)
partition_reads <- function(x, cl_method="ward.D", min_len = 0.8, skip_gap_freq = 2/3){
  # get SNPs
  ppos <- colnames(x)
  xm <- x[order(rownames(x)),]
  bad_ppos <- apply(xm, 2, function(x) NROW(x[x == "-"])/NROW(x) > skip_gap_freq )
  xm <- xm[,!bad_ppos]
  bad_ppos <- ppos[bad_ppos]

  ## Get the SNP matrix as sequences
  xseqs <- apply(xm, 1, function(t) c(unlist(paste(t, collapse = ""))))
  xseqs <- Biostrings::DNAStringSet(xseqs)
  # Try hamming dist from decipher without using terminal gaps and also penalize all gaps vs gaps
  xmat <- DECIPHER::DistanceMatrix(xseqs, includeTerminalGaps = FALSE, penalizeGapLetterMatches = TRUE, penalizeGapGapMatches = TRUE, correction = "none", processors = 4)


  ## Get only the fraction of reads that contain at least min_len of total SNPs
  # might be evaluated whether really necessary; Seems so. Produces way better clustering!
  # a good cutoff seems around 0.8
  clustres <- get_clusts(xseqs, xmat, cl_method = cl_method, min_len = min_len)
  subclades <- clustres$clades
  tree <- clustres$tree

  # get distance of each read to all clades and assign the read to the clade with min distance
  clades <- dplyr::data_frame(read = names(xseqs))
  subclades <- dplyr::data_frame(clustclade = subclades, read = names(subclades))
  clades <- dplyr::full_join(clades, subclades, by = "read")
  clades <- assign_partition(clades, xmat)

  # get cluster quality measure; coefficients
  clades <- clust_quality(clades, xmat)

  ## add consenus mapping for read classification
  ## Get consensus sequences for all haplotypes from initial clustering
  #conseqs <- Biostrings::DNAStringSet(lapply(levels(subclades), function(x) get_Hap_consensus(x, xseqs, subclades)))
  #names(conseqs) <- levels(subclades)
  #xconseqs <- c(xseqs, conseqs)

  ## Assign haplotypes to missing reads
  #xmat <- as.matrix(Biostrings::stringDist(xconseqs, method = "hamming"))
  #clades <- unlist(lapply(names(xseqs), function(x) assign_partition(xmat[x, names(conseqs)])))
  #clades<- dplyr::data_frame(read = names(xseqs), clade = clades)
  #subclades <- dplyr::data_frame(clustclade = subclades, read = names(subclades))


  clades <- clades %>%
  dplyr::mutate(correct = dplyr::if_else(is.na(clustclade), NA, dplyr::if_else(clustclade == clade, TRUE, FALSE)))
  ## Correctly classified clades in the initial clustering
  falseClassified <- NROW(dplyr::filter(clades, correct==FALSE))/NROW(dplyr::filter(clades, correct == TRUE))
  message(sprintf( " Corrected classification of %.2f%% reads", 100*falseClassified))

  # Create the partition table
  part_ <- HapPart(read_name = clades$read, snp_pos = colnames(x))
  PRT(part_) <- clades$clade
  HTR(part_) <- tree
  SNP(part_) <- ppos
  K(part_) <- length(ppos)-length(bad_ppos)
  mcoef(part_) <- clades$coeff

  return(part_)
}

# Helpers -----------------------------------------------------------------

get_clusts <- function(xseqs, xmat, min_len = 0.8, cl_method = "ward.D"){
  x_sub <- xseqs[Biostrings::width(gsub("-", "", xseqs)) > min_len*Biostrings::width(xseqs[1])]

  # get distance only for the subset
  submat <- xmat[names(x_sub), names(x_sub)]
  subdist <- dist(submat)
  #subdist <- Biostrings::stringDist(x_sub, method="hamming")
  #submat <- as.matrix(subdist)


  ## Perform a hierarchical clustering
  hcc <-  hclust(subdist, method = cl_method)
  plot(hcc, labels=FALSE)
  ## do a dynamic cut. Need to be evaluated
  clusts <- dynamicTreeCut::cutreeHybrid(hcc, distM = submat, deepSplit = 2)
  # extract the clusters for each sequence
  clades <- as.factor(unlist(lapply(unname(clusts$labels), function(i) rawToChar(as.raw(as.integer(i)+64)))))
  names(clades) <- hcc$labels

  return(list(clades = clades, tree = hcc))
}

# Get quality measure of all vs all clusters
assign_partition <- function(clades, xmat){

  # Get means for all clusters from distance matrix
  clade_sums <- clades
  hptypes <- levels(clades$clustclade)
  for (cld in hptypes){
    ownclade <- dplyr::filter(clades, clustclade == cld)$read
    sums_ownclade <- unlist(lapply(clades$read, function(t) mean(xmat[t,][names(xmat[t,]) %in% ownclade])))
    clade_sums <- dplyr::bind_cols(clade_sums, !!cld := sums_ownclade)
  }
  # get only proper reads, i.e. not NA in the dist matrix. That happens mostly to reads size == 1
  clade_sums <- dplyr::filter(clade_sums, !is.na(A))
  clade <- unlist(lapply(clade_sums$read, function(a) get_haptype(a, clade_sums)))
  dplyr::bind_cols(clade_sums, clade = clade)
}


clust_quality <- function(clades, xmat){
  ## Calculate coefficients for each vs each cluster for each read
  scored_clade_sums <- dplyr::data_frame()
  hptypes <- levels(clades$clustclade)
  for (cld in hptypes){
    tmpcld <- dplyr::filter(clades, clade == cld)
    # The main coefficient is the distance of own cluster to that of all others
    coeff <- tmpcld[cld]/rowMeans(tmpcld[hptypes[!hptypes == cld]])
    names(coeff) <- "coeff"
    tmpcld <- dplyr::bind_cols(tmpcld, coeff)
    for (ocld in hptypes){
      cnames <- c(names(tmpcld), paste0("coeff", ocld))
      ## The finer coefficients are the distance of own cluster to those of all others separately
      indCoef <- tmpcld[cld]/tmpcld[ocld]
      tmpcld <- dplyr::bind_cols(tmpcld, replace(indCoef, indCoef == 1, 0))
      names(tmpcld) <- cnames
    }
    scored_clade_sums <- dplyr::bind_rows(scored_clade_sums, tmpcld)
  }
  dplyr::arrange(scored_clade_sums, read)
}


get_haptype <- function(name, clade_sums){
  entry <- dplyr::filter(clade_sums, read == name)
  value <- min(entry[hptypes])
  names(entry[which(entry == value)])
}

# get_Hap_consensus <- function(cld, x_, clades){
#   seqsCld <- x_[names(x_) %in% names(clades[clades == cld])]
#   seq <- Biostrings::DNAString(Biostrings::consensusString(seqsCld, ambiguityMap="N", threshold = 0.5))
#   seq
# }

#assign_partition <- function(mat){
#  haplotypes <- names(mat[mat == min(mat)])
#  if (length(haplotypes) > 1){sample(haplotypes, 1)}else{haplotypes}
#}

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
    geom_histogram(aes(x = mcoef, fill = haplotype), binwmeidth = 0.05) +
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


optimal_partition_limits <- function(scores) {
  coeff_cuts <- seq(0.05, 1.10, by = 0.05)
  rHap <- lapply(scores, function(x) sapply(coeff_cuts, function(cutoff) sum(x <= cutoff)))

  df <- dplyr::data_frame()
  for(hapType in names(rHap)){
    df <- dplyr::bind_rows(df,
    dplyr::data_frame(c = coeff_cuts, r = unlist(rHap[hapType])) %>%
      dplyr::mutate(haplotype = hapType, score = (r/c)) %>%
      dplyr::mutate(max.score = score == max(score)))
  }

  dfmax <- dplyr::filter(df, max.score) %>%
    dplyr::select(haplotype, c, r, score)

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
