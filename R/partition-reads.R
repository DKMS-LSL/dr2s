

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

# debug
 # x <- mat
 # cl_method="ward.D"
 # min_len = 0.8
 # skip_gap_freq = 2/3
 # deepSplit = 1

#   x <- dedk$partition
#x <-   dedk.part$partition$mat
#x <- matLR
#nrow(matLR)
#nrow(matSR)
#x <- dedk.split$partition$mat
#x <- dl1_2$partition$mat


partition_reads <- function(x, cl_method="ward.D", min_len = 0.8, skip_gap_freq = 2/3, deepSplit = 1, outdir = "./outdir"){
  # get SNPs
  ppos <- colnames(x)
  xm <- x[order(rownames(x)),]
  bad_ppos <- apply(xm, 2, function(x) NROW(x[x == "-"])/NROW(x) > skip_gap_freq )
  xm <- xm[,!bad_ppos]
  bad_ppos <- ppos[bad_ppos]

  ## Get the SNP matrix as sequences
  xseqs <- get_seqs_from_mat(xm)

  #HMM
    # # Write xseqs as an aligned FASTA file
    # seqfile = file.path(outdir, "partSeqs.fa")
    # Biostrings::writeXStringSet(xseqs, seqfile)

  # Try hamming dist from decipher without using terminal gaps and also penalize all gaps vs gaps
  # Swith to a custom Position Specific Distance Matrix
  # xmat <- DECIPHER::DistanceMatrix(xseqs,
  #                                  includeTerminalGaps = FALSE,
  #                                  penalizeGapLetterMatches = FALSE,
  #                                  penalizeGapGapMatches = FALSE,
  #                                  correction = "none",
  #                                  processors = 4,
  #                                  verbose = FALSE)

  ## Get only the fraction of reads that contain at least min_len of total SNPs
  clustres <- get_clusts(xseqs, cl_method = cl_method, min_len = min_len, deepSplit = deepSplit)
  subclades <- clustres$clades
  tree <- clustres$tree
  # debug
  # plot(tree, labels = F)
  hptypes <- levels(subclades)

  # Get scores and assign clades by freq mat
  ## Position Weight Matrix: Use frequency plus pseudocount/ basefrequency (here 0.25 for each).
  msa <- lapply(levels(subclades), function(x) xseqs[names(subclades[subclades == x])])
  names(msa) <- hptypes
  mats <- lapply(msa, function(x) create_PWM(x))

  # ToDo refine chimera finding
  seqs <- Biostrings::DNAStringSet(sapply(msa, function(x) unlist(simple_consensus(t(Biostrings::consensusMatrix(x)[c(VALID_DNA(), "-", "+"),])))))
  rC <- sort(find_chimeric(seqs))
  mats <- mats[rC]

  scores <- dplyr::bind_rows(lapply(1:length(xseqs),function(s) get_scores(s, xseqs, mats)))

  # ToDo: implement lenient clade assignment
  # srpartition %>%
  #   dplyr::group_by(read, haplotype) %>%
  #   dplyr::mutate(clade = prod(prob)) %>%
  #   dplyr::ungroup() %>%
  #   dplyr::group_by(read) %>%
  #   dplyr::mutate(max = max(clade)) %>%
  #   dplyr::group_by(read, haplotype) %>%
  #   dplyr::select(read, haplotype, clade, max) %>%
  #   dplyr::distinct() %>%
  #   dplyr::group_by(read) %>%
  #   dplyr::filter(abs(1-(clade/max)) < 0.1) %>%
  #   dplyr::mutate(exclusive = n() > 1) %>%
  #   dplyr::ungroup()

  # debug
  #hist(dplyr::filter(clades, clade == "A")$score)
  clades <- scores %>%
    dplyr::group_by(read) %>%
    dplyr::slice(which.max(score)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(clustclade = ifelse(read %in% names(subclades),
                                      as.character(subclades[read]),
                                      NA)) %>%
    dplyr::mutate(correct = dplyr::if_else(clustclade == clade, TRUE, FALSE))

  ######## Version using HMMER #######################
  ##
  # Write seqs as fasta for each cluster
  #HMM
  # partFiles <- lapply(levels(subclades), function(x) file.path(outdir, paste(c(x, "fa"), collapse = ".")))
  # names(partFiles) <- levels(subclades)
  # lapply(levels(subclades), function(x) Biostrings::writeXStringSet(xseqs[names(subclades[subclades == x])], partFiles[[x]]))

  # Make HMM for each cluster and map/search the long reads
  # extClusters <- extend_clusters(partFiles, seqfile, outdir, force)

  # cladeshmm <- extClusters %>%
  #   dplyr::group_by(read) %>%
  #   dplyr::slice(which.max(score)) %>%
  #   dplyr::ungroup() %>%
  #   dplyr::mutate(clustclade = ifelse(read %in% names(subclades),
  #                                     as.character(subclades[read]),
  #                                     NA)) %>%
  #   dplyr::mutate(correct = dplyr::if_else(clustclade == clade, TRUE, FALSE))
  # falseClassifiedhmm <- NROW(dplyr::filter(cladeshmm, correct==FALSE))/NROW(dplyr::filter(cladeshmm, correct == TRUE))
  # message(sprintf( " Corrected classification of %.2f%% reads", 100*falseClassifiedhmm))
  #####################################################


  ## Correctly classified clades in the initial clustering
  falseClassified <- NROW(dplyr::filter(clades, correct==FALSE))/NROW(dplyr::filter(clades, correct == TRUE))
  message(sprintf( " Corrected classification of %.2f%% reads", 100*falseClassified))

  # Create the partition table
  part_ <- HapPart(read_name = clades$read, snp_pos = colnames(x))
  PRT(part_) <- clades$clade
  HTR(part_) <- tree
  SNP(part_) <- ppos
  K(part_) <- length(ppos)-length(bad_ppos)
  oc(part_) <-as.character(unique(subclades))
  mcoef(part_) <- clades$score
  scores(part_) <- scores
#  mats(part_) <- mats


  return(part_)
}

## debug ##
# part_7 <-partition_reads(x, cl_method="ward.D", min_len = 0.7, skip_gap_freq = 2/3, deepSplit = 1)
#  part_8 <-partition_reads(x, cl_method="ward.D", min_len = 0.8, skip_gap_freq = 2/3, deepSplit = 1)
#  part_9 <-partition_reads(x, cl_method="ward.D", min_len = 0.91, skip_gap_freq = 2/3, deepSplit = 1)
#  levels(partition(part_7)$haplotype)
#  levels(partition(part_8)$haplotype)
#  levels(partition(part_9)$haplotype)
#  HTR(part_7)
#  HTR(part_8)
#  b <- HTR(part_9)
#  length(b$labels)/ nrow(x)
#
#  plot_partition_tree(part_7)
#  plot_partition_tree(part_8)
#  plot_partition_tree(part_9)
#
#  part_
#  part_9 <-partition_reads(x, cl_method="ward.D", min_len = 0.93, skip_gap_freq = 2/3, deepSplit = 1)
#  part_a <-partition_reads(x, cl_method="ward.D", min_len = 0.85, skip_gap_freq = 2/3, deepSplit = 1)
#  plot_partition_tree(part_a)
#  HTR(part_a)
# # #self <- dedk.part
# # #self$getPartition()
#  x <- run$partition$mat
#  x <- dedk$partition$mat
#x <- hla.part$partition$mat

#x
#class(x)
#  cl_method="ward.D"
#  min_len = 0.8
# skip_gap_freq = 2/3
# min_reads_frac = 1/3
#deepSplit = 1
# partn <- partition_reads(x)
#length(partn[mcoef(partn)>1])
#length(partn)
#min_reads_frac <- 1/3

# Helpers -----------------------------------------------------------------

get_clusts <- function(xseqs, xmat, min_len = 0.80, cl_method = "ward.D", deepSplit = 1, min_reads_frac = 1/3){

  # Need more sequences for more snps
  adjusted_min_reads_frac <- min_reads_frac + (0.03*Biostrings::width(xseqs[1])/50)
  # get long reads
  min_lens <- seq(1, .7, -0.01)
  len_counts <- sapply(min_lens, function(minl) length(xseqs[Biostrings::width(gsub("-", "", xseqs)) > minl*Biostrings::width(xseqs[1])])/length(xseqs))
  names(len_counts) <- min_lens
  min_len <- max(as.numeric(names(len_counts[which(len_counts>adjusted_min_reads_frac)][1])), min_len)

  message("Using only long reads containing ", min_len*100, "% of all SNPs")
  x_sub <- xseqs[Biostrings::width(gsub("-", "", xseqs)) > min_len*Biostrings::width(xseqs[1])]

  ## Consensus matrix with pseudocount
  message("  Constructing a Position Specific Distance Matrix of ", length(xseqs), " sequences ...")
  consmat  <- Biostrings::consensusMatrix(x_sub, as.prob = TRUE)[VALID_DNA(),] + 1/length(xseqs)
  ## Create seq matrix as input for cpp_PSDM
  x_sub_tmp <- as.matrix(x_sub)
  x_sub_mat<- plyr::revalue(x_sub_tmp, c("G" = 1, "A" = 2, "T" = 3, "C" = 4, "-" = 5))
  x_sub_mat<- sapply(x_sub_mat, as.numeric)
  dim(x_sub_mat) <- dim(x_sub_tmp)
  rm(x_sub_tmp)

  ## Get Position Specific Distance Matrix
  dist <- cpp_PSDM(consmat, x_sub_mat)
  colnames(dist) <- names(x_sub)
  rownames(dist) <- names(x_sub)
  dist <- as.dist(dist)

  ## Perform a hierarchical clustering
  hcc <-  hclust(dist, method = cl_method)
  ## do a dynamic cut. Need to be evaluated
  clusts <- dynamicTreeCut::cutreeHybrid(hcc, distM = as.matrix(dist), deepSplit = deepSplit)
  # extract the clusters for each sequence
  clades <- as.factor(sapply(unname(clusts$labels), function(i) rawToChar(as.raw(as.integer(i)+64))))
  names(clades) <- hcc$labels
  return(list(clades = clades, tree = hcc))
}

get_scores <- function(s, xseqs, mats){
  seq <- as.character(xseqs[[s]])
  seq <- unlist(strsplit(seq, split = ""))
  read <- names(xseqs[s])

  b <- sapply(mats, function(t) sum( sapply(1:length(seq), function(x) t[seq[x],x] ) ) )
  t <- data.frame(read, b, names(b))#, max(b), names(which(b == max(b)))))
  t$read <- as.character(t$read)
  names(t) <- c("read", "score", "clade")
  t
}

find_chimeric <- function(seqs) {
  forward <- lapply(1:(Biostrings::width(seqs[1])/2), function(x) DECIPHER::DistanceMatrix(Biostrings::DNAStringSet(seqs, start = 1, width = x), verbose = FALSE))
  reverse <- lapply(1:(Biostrings::width(seqs[1])/2), function(x) DECIPHER::DistanceMatrix(Biostrings::DNAStringSet(Biostrings::reverse(seqs), start = 1, width = x), verbose = FALSE))
  f <- DECIPHER::DistanceMatrix(Biostrings::DNAStringSet(seqs, start = 1, width = Biostrings::width(seqs)/2), verbose = FALSE)
  r <- DECIPHER::DistanceMatrix(Biostrings::DNAStringSet(Biostrings::reverse(seqs), start = 1, width = Biostrings::width(seqs)/2), verbose = FALSE)
  seqs
  dm <- foreach(h = names(seqs), .combine = 'cbind') %:%
   foreach(hp = names(seqs)) %do% {
    f <- sapply(forward, function(x) x[h,hp])
    r <- sapply(reverse, function(x) x[h,hp])
    mean(f+r)
   }
  rownames(dm) <- names(seqs)
  colnames(dm) <- names(seqs)
  addC <- names(dm[1,][dm[1,]>1])
  c("A", addC)

  #rC <-arrayInd(which(dm == unlist(dm[dm > 1])), dim(dm))
  #rC <- rownames(dm)[rC[,1]]
}

## TODO: RM ????
extend_clusters <- function(partFiles, seqfile, outdir, force = FALSE) {

  # Construct hmm for each haplotype
  hmmFiles <- lapply(partFiles, function(x) run_hmmer("hmmbuild", x, outdir, force = TRUE))
  # Merge and prepare hmm database
  hmmDB <- run_hmmer("hmmpress", lapply(hmmFiles, function(x) x$outfile), outdir, force = TRUE)
  # Search and score all seqs
  hmmerResFile <- run_hmmer("hmmscan", list(db = hmmDB$outfile, seqs = seqfile), outdir, force = TRUE)
  hmmerRes <- readr::read_table2(hmmerResFile$outfile, comment = "#", col_names = FALSE)
  hmmerRes <- hmmerRes %>%
    dplyr::select(X3, X1, X9)
  names(hmmerRes) <- c("read", "clade", "score")
  hmmerRes
}
#
# # ToDo: RM
# # Get quality measure of all vs all clusters
# assign_partition <- function(clades, xmat){
#   # Get means for all clusters from distance matrix
#   clade_sums <- clades
#   hptypes <- levels(clades$clustclade)
#
#   for (cld in hptypes){
#     ownclade <- dplyr::filter(clades, clustclade == cld)$read
#     sums_ownclade <- sapply(clades$read, function(t) mean(xmat[t,][names(xmat[t,]) %in% ownclade]))
#     clade_sums <- dplyr::bind_cols(clade_sums, !!cld := sums_ownclade)
#   }
#   # get only proper reads, i.e. not NA in the dist matrix. That happens mostly to reads size == 1
#   proper_reads <- clade_sums %>%
#     dplyr::select(read,hptypes) %>%
#     dplyr::filter(complete.cases(.)) %>%
#     .$read
#   clade_sums <- dplyr::filter(clade_sums, read %in% proper_reads)
#
#   # get haplotype
#   clade <- sapply(clade_sums$read, function(a) get_haptype(a, clade_sums))
#
#   dplyr::bind_cols(clade_sums, clade = clade)
# }
#
# # Todo: RM
# clust_quality <- function(clades, xmat){
#   ## Calculate coefficients for each vs each cluster for each read
#   scored_clade_sums <- dplyr::data_frame()
#   hptypes <- levels(clades$clustclade)
#   for (cld in hptypes){
#     tmpcld <- dplyr::filter(clades, clade == cld)
#     # The main coefficient is the distance of own cluster to that of all others
#     coeff <- 1-tmpcld[cld]/rowMeans(tmpcld[hptypes[!hptypes == cld]])
#     names(coeff) <- "coeff"
#     tmpcld <- dplyr::bind_cols(tmpcld, coeff)
#     for (ocld in hptypes){
#       cnames <- c(names(tmpcld), paste0("coeff", ocld))
#       ## The finer coefficients are the distance of own cluster to those of all others separately
#       indCoef <- 1-tmpcld[cld]/tmpcld[ocld]
#       tmpcld <- dplyr::bind_cols(tmpcld, replace(indCoef, indCoef == 1, 0))
#       names(tmpcld) <- cnames
#     }
#     scored_clade_sums <- dplyr::bind_rows(scored_clade_sums, tmpcld)
#   }
#   dplyr::arrange(scored_clade_sums, read)
# }
#
# #Todo RM
# get_haptype_no_reassign <- function(name, clade_sums){
#   entry <- dplyr::filter(clade_sums, read == name)
#   if (!is.na(entry$clustclade)){
#    return(as.character(entry$clustclade))
#   }else{
#     value <- min(entry[levels(clade_sums$clustclade)])
#     return(names(entry[which(entry == value)]))
#   }
# }
#
# # Todo: RM
# get_haptype <- function(name, clade_sums){
#   entry <- dplyr::filter(clade_sums, read == name)
#   value <- min(entry[levels(clade_sums$clustclade)])
#   return(names(entry[which(entry == value)]))
# }

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
    geom_vline(xintercept = c(limits["A"], -limits["B"]), linetype = "dashed", colour = "grey80") +
    # geom_vline(xintercept = c(-0.75, -0.5, -0.25, 0.25, 0.5, 0.75),
    #            linetype = "dashed", size = 0.5, colour = "grey80") +
    # xlim(c(-1.1, 1.1)) +
    xlab("Haplotype membership coefficient") +
    ylab("Number of reads") +
    ggtitle(label = label) +
    theme_bw()
}

#' Plot distribution of haplotype partitioned reads for more than two haplotypes
#'
#' @param x A \code{HapPart} object.
#' @param label Optional plot label.
#'
#' @return A \code{ggplot} object.
#' @export
#' @examples
#' ###
#'
#x <- dedk.part$getPartition()
#limits <- dedk.part$getLimits()
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
                     linetype = "dashed", size = 1)+
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


plot_radar_partition <- function(x){
  stopifnot(is(x, "HapPart"))

  df <- dplyr::full_join(scores(x), partition(x), by = "read")

  # use bigger size for 2d radarplot.
  size <- ifelse(length(levels(df$haplotype)) == 2, 4, 0.05)
  ggplot(df, aes(x = clade, y = score)) +
    geom_polygon(aes(group = read, color = haplotype), fill = NA, size = size, show.legend = FALSE, alpha = 1) +
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

#' Plot tree of first clustering
#'
#' @param x A \code{HapPart} object.
#'
#' @return A \code{ggplot} object.
#' @export
#' @examples
#' ###
#x <- self$getPartition()
plot_partition_tree <- function(x){
  stopifnot(is(x, "HapPart"))
  # ToDo: test if packages are loaded
  tree <- HTR(x)
  k <- length(oc(x))

  dendr <- ggdendro::dendro_data(tree, type = "rectangle")
  dendr$labels <- dendr$labels %>%
    dplyr::mutate(label = as.character(label))

  clust <- cutree(tree, k)
  clust <- dplyr::data_frame(label = names(clust), haplotype = clust)
  # clust <- partition(x) %>%
  #   dplyr::rename(label = read) %>%
  #   dplyr::select(-mcoef)
  #
  # clust <- scores(x) %>%
  #   dplyr::filter(!is.na(clustclade)) %>%
  #   dplyr::select(read, clustclade) %>%
  #   dplyr::rename( haplotype = clustclade, label = read)
  #
  dendr$labels <- dplyr::left_join(dendr$labels, clust, by="label")

  height <- unique(dendr$segments$y)[order(unique(dendr$segments$y), decreasing = TRUE)]
  cut.height <- mean(c(height[k], height[k-1]))
  dendr$segments <- dendr$segments %>%
    dplyr::mutate(line =  dplyr::if_else(
      y == yend & y > cut.height, 1, dplyr::if_else(
        yend > cut.height,1, 2)))

  dendr$segments$cluster <- c(-1, diff(dendr$segments$line))
  change <- which(dendr$segments$cluster == 1)
  for (i in 1:k) dendr$segments$cluster[change[i]] = i + 1
  dendr$segments <- dendr$segments %>%
    dplyr::mutate(cluster = dplyr::if_else( line == 1, 1, ifelse( cluster == 0, NA, cluster)))
  dendr$segments$cluster <- sapply(1:NROW(dendr$segments$cluster), function(x) getCl(x, dendr$segments$cluster, change))

  # Correct order

  labs <- c("N", oc(x))#LETTERS[as.numeric(names(sort(table(clust$haplotype), decreasing = TRUE)))])
  dendr$segments$cluster <- labs[dendr$segments$cluster]
  # Make plot
  p <- ggplot() +
    geom_segment(data=ggdendro::segment(dendr), aes(x = x, y = y, xend = xend, yend = yend, color = factor(cluster)), size = 1.25) +
    scale_color_manual(values = PARTCOL(),
                          name = "Haplotype",
                          breaks = LETTERS[1:k],
                          labels =  LETTERS[1:k]) +
    geom_hline(yintercept = cut.height, color = "blue") +
    ggtitle("Initial clustering of haplotypes") +
    # ggdendro::theme_dendro() +
    theme_bw() +
    theme(legend.position = "bottom",
          axis.title = element_blank(),
          axis.text  = element_blank(),
          axis.ticks  = element_blank()
          )
  p
}


# Plot helper
getCl <- function(n, cluster, change) {
  ifelse(!is.na(cluster[n]),
         cluster[n],
         cluster[change[max(which(change < n))]]
  )
}
# Helper function for the plotting of radar charts
coord_radar <- function (theta = "x", start = 0, direction = 1)  {
  theta <- match.arg(theta, c("x", "y"))
  r <- if (theta == "x")
    "y"
  else "x"
  ggproto("CoordRadar", CoordPolar, theta = theta, r = r, start = start,
          direction = sign(direction),
          is_linear = function(coord) TRUE)
}

# Optimal partitioning ----------------------------------------------------


optimal_partition_limits <- function(scores, f = 0.3) {
  coeff_cuts <- seq(0, max(unlist(scores)), length.out = 100)
  rHap <- lapply(scores, function(x) sapply(coeff_cuts, function(cutoff) sum(x >= cutoff)))
  df <- dplyr::data_frame()

  f <- 0.3
  df <- data.frame()
  for(hapType in names(rHap)){
    l <- f * length(unlist(scores[!names(rHap) == hapType]))/ length(scores[[hapType]])
    df <- dplyr::bind_rows(df,
    dplyr::data_frame(c = coeff_cuts, r = unlist(rHap[hapType])) %>%
      dplyr::mutate(haplotype = hapType, score = ((r^l)*c)) %>%
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
