
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

  message("Partition shortreads on ", length(ppos), " positions")
  res <- foreach (pos = ppos, .combine = "rbind") %dopar% {
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
score_highest_SR <- function(srpartition, diffThreshold = 0.2) {
  srpartition %>%
    dplyr::group_by(read, haplotype) %>%
    dplyr::mutate(clade = prod(prob)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(read) %>%
    dplyr::mutate(max = max(clade)) %>%
    dplyr::group_by(read, haplotype) %>%
    dplyr::select(read, haplotype, clade, max) %>%
    dplyr::distinct() %>%
    dplyr::group_by(read) %>%
    dplyr::filter(abs(1-(clade/max)) < diffThreshold) %>%
    dplyr::mutate(exclusive = n() == 1) %>%
    dplyr::ungroup()
}

write_part_fq <- function(fq, srFastqHap, dontUseReads = dontUseReads) {
  file_delete_if_exists(srFastqHap)
  fqstream = ShortRead::FastqStreamer(fq)
  repeat {
    sr = ShortRead::yield(fqstream)
    if (length(sr) == 0) break
    fqnames <- as.character(ShortRead::id(sr))
    fqnames <- sub(" .*$", "", fqnames)
    sr <- sr[which(!fqnames %in% dontUseReads)]
    ShortRead::writeFastq(sr, srFastqHap, mode="a", compress = TRUE)
  }
  close(fqstream)
}
# --- Helper ---
part_read <- function(a, mats, pos){
  l <- sapply(mats, function(x) x[,pos][as.character(a)])
  names(l) <- names(mats)
  list(prob = l,  haplotype = names(l), pos = rep(pos, length(l)))
}

