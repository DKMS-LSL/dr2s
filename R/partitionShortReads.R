
# Partition ShortReads ------------------------------------------------------------------

#' Get a score for each shortread in a bam file on specific positions. Scores based on a consensus matrix of longreads.
#'
#' @param ppos List of polymorphic positions used for initial clustering.
#' @param refname Name of the reference.
#' @param bamfile BAM file path.
#' @param mats Consensus matrix from longread clustering at positions in ppos.
#' @param cores Number of cores or "auto"
#'
#' @details
#' Returns a \code{srpartition} object:
#' A \code{list} with slots:
#' \describe{
#'   \item{bamfile}{<character>; Path to the bam file used to construct the pileup}
#'   \item{ppos}{<numeric>; Positions used for scoring}
#'   \item{srpartition}{A \code{data.frame} with colums:
#'     \describe{
#'       \item{read}{<character>; The read name}
#'       \item{pos}{<integer>; Genomic position of base}
#'       \item{Haplotype}{<character>; Haplotype}
#'       \item{prob}{<integer>; probability for a base at a position}
#'     }
#'   }
#'   \item{consmat}{Consensus matrix}
#' }
#'
#' @return A \code{srpatition} object. See \code{Details}.
#' @export
#' @examples
#' ###
get_SR_partition_scores <- function(ppos, refname, bamfile, mats, cores = "auto"){
  stopifnot(
    requireNamespace("parallel", quietly = TRUE),
    requireNamespace("doParallel", quietly = TRUE),
    requireNamespace("GenomicAlignments", quietly = TRUE)
  )

  # Register parallel worker
  if (cores == "auto") {
    cores <- parallel::detectCores()
    cores <- ifelse(is.na(cores), 1, cores/2)
  }
  stopifnot(is.numeric(cores))
  doParallel::registerDoParallel(cores = cores)

  flog.info("  Partition shortreads on %s positions ...", length(ppos),  name = "info")
  res <- foreach (pos = ppos, .combine = "rbind") %dopar% {
    message("Get reads on position ", pos)
    param <- paste(paste(refname, pos, sep = ":"), pos, sep = "-")
    stack <- GenomicAlignments::stackStringsFromBam(bamfile, param = param, use.names = TRUE)
    dplyr::bind_cols(read = rep(names(stack), each = length(mats)), dplyr::bind_rows(lapply(stack, function(x) part_read(x, mats, pos))))
  }
  structure(
    list(
      bamfile     = bamfile,
      ppos        = ppos,
      srpartition = res,
      consmat     = mats
    ),
    class = "srpartition"
  )
}

##' Get the highest scoring haplotype for each read
#'
#' @param srpartition srpartition dataframe with scores per read and haplotype.
#' @param diffThreshold report reads which scores differ only below threshold to both haplotypes
#' @details
#' Returns a \code{data.frame} object with columns:
#' \describe{
#'   \item{read}{<character>; The read name}
#'   \item{Haplotype}{<character>; Haplotype}
#' }
#'
#' @return A \code{data.frame} object. See \code{Details}.
#' @export
#' @examples
#' ### Score
#
score_highest_SR <- function(srpartition, diffThreshold = 0.001) {
  # debug
# srpartition <- srpartitionbc$srpartition
# srpartitionbc <- srpartition

# microbenchmark::microbenchmark(
#   sr <- srpartition %>%
#     dplyr::group_by(read, haplotype) %>%
#     dplyr::mutate(clade = sum(prob)) %>%
#     dplyr::ungroup() %>%
#     dplyr::group_by(read) %>%
#     dplyr::mutate(max = max(clade)) %>%
#     dplyr::group_by(read, haplotype) %>%
#     dplyr::select(read, haplotype, clade, max) %>%
#     dplyr::distinct(),
#   a <- unique(s[, clade := sum(prob), by = list(read, haplotype)] # Get the sum of each read and hptype
#               [,max := max(clade), by = read] # get the max of the sums of each
#               [, !c("prob", "pos")]) # dismiss the prob and pos which we dont need anymore
# )

# microbenchmark::microbenchmark(
#   sr2 <- sr %>%
#     dplyr::group_by(read) %>%
#     dplyr::filter(abs(1-(clade/max)) < diffThreshold) %>%
#     dplyr::mutate(exclusive = n() == 1) %>%
#     dplyr::ungroup(),
#   b <- data.table::as.data.table(sr)[abs(1-(clade/max)) < diffThreshold][, exclusive:= .N, by = read], times = 2
# )

  ## Using data.table is WAY faster than dplyr! 100x
  sr <- unique(data.table::as.data.table(srpartition)[, clade := sum(prob), by = list(read, haplotype)] # Get the sum of each read and hptype
              [,max := max(clade), by = read] # get the max of the sums of each
              [, !c("prob")]) # dismiss the prob and pos which we dont need anymore
  sr2 <- sr[abs(1-(clade/max)) < diffThreshold][, exclusive:= .N, by = read]


  flog.info(" Calculate shortread scoring with cutoff: %s ...", diffThreshold, name = "info")
  correctScoring <- sapply(unique(sr2$pos), function(position) check_SR_scoring(position, sr2))

  if (!all(correctScoring))
    flog.info(" Scoring not meeting proper clustering, refine diffThreshold ...", diffThreshold, name = "info")

  while(!all(correctScoring)){
    diffThreshold <- diffThreshold + 0.01
    flog.info(" Calculate shortread scoring with cutoff: ", diffThreshold, name = "info")
    # Check again with new threshold if there are enough reads at each position
    sr2 <- sr[abs(1-(clade/max)) < diffThreshold][, exclusive:= .N, by = read]
    correctScoring <- sapply(unique(sr2$pos), function(position) check_SR_scoring(position, sr2))
  }
  flog.info(" Use overlap cutoff of %s for shortread scoring ...", diffThreshold, name = "info")
  return(sr2)
}

write_part_fq <- function(fq, srFastqHap, dontUseReads = NULL, useReads = NULL) {
  fqstream = ShortRead::FastqStreamer(fq)
  file_delete_if_exists(srFastqHap)
  repeat {
    sr = ShortRead::yield(fqstream)
    if (length(sr) == 0) break
    fqnames <- as.character(ShortRead::id(sr))
    fqnames <- sub(" .*$", "", fqnames)
    if (is.null(dontUseReads)) {
      useReads <- which(fqnames %in% useReads)
    } else {
      # useReads = qnames
      useReads <- which(!fqnames %in% dontUseReads)
    }
    flog.info("  Using %s of %s reads", length(useReads), length(fqnames), name = "info")
    sr <- sr[useReads]
    ShortRead::writeFastq(sr, srFastqHap, mode="a", compress = TRUE)
  }
  close(fqstream)
}
# --- Helper ---

part_read <- function(a, mats, pos){
  l <- sapply(mats, function(x) x[,pos][as.character(a)])
  list(prob = l,  haplotype = names(mats), pos = rep(pos, length(l)))
}

check_SR_scoring <- function(position, sr2){
  reads_at_pos <- unique(sr2[pos == position]$read)
  score_ok <- sapply(unique(sr2$haplotype), function(hp) sum(reads_at_pos %in% sr2[haplotype == hp]$read)/length(reads_at_pos) > 0.001)
  return(score_ok)
}
