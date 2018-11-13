
# Partition Shortreads ---------------------------------------------------------

#' Get a score for each shortread in a bam file on specific positions.
#' Scores based on a consensus matrix of longreads.
#'
#' @param bamfile BAM file path.
#' @param mats Consensus matrix from longread clustering at positions in ppos.
#' @details
#' Returns a \code{srpartition} object:
#' A \code{list} with slots:
#' \describe{
#'   \item{bamfile}{<character>; Path to the bam file used to construct the
#'   pileup}
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
getSRPartitionScores <- function(bamfile, mats, ...) {
  assert_that(
    requireNamespace("parallel", quietly = TRUE),
    requireNamespace("doParallel", quietly = TRUE),
    requireNamespace("GenomicAlignments", quietly = TRUE))

  indent <- list(...)$indent %||% indentation()

  ## make sure that all available cores get utilised
  op <- options("dr2s.max.cores")
  if (!is.null(op$dr2s.max.cores)) {
    options(dr2s.max.cores = NULL)
    on.exit(options(op))
  }

  ## Get polymorphic positions
  ppos <- colnames(mats[[1]])
  ## Register as many workers as necessary or available
  workers <- min(length(ppos), .getIdleCores())
  bpparam <- BiocParallel::MulticoreParam(workers = workers)
  bam <- Rsamtools::BamFile(bamfile)
  ## Get reference from bamfile
  refname <- GenomeInfoDb::seqnames(Rsamtools::seqinfo(bam))
  if (workers > 1) {
    Rsamtools::open.BamFile(bam)
    on.exit(Rsamtools::close.BamFile(bam))
  }
  flog.info("%sUsing %s workers to partition shortreads on %s positions",
            indent(), workers, length(ppos), name = "info")
  res <- do.call(dplyr::bind_rows,
    suppressWarnings(BiocParallel::bplapply(as.integer(ppos),
    function(pos, refname, bamfile, mats) {
      #message("Get reads on position ", pos)
      param <- GenomicRanges::GRanges(
        seqnames = refname, ranges = IRanges::IRanges(start = pos, end = pos))
      stack <- GenomicAlignments::stackStringsFromBam(bamfile, param = param, use.names = TRUE)
      dplyr::bind_cols(read = rep(names(stack), each = length(mats)),
                       dplyr::bind_rows(lapply(stack, function(x)
                         .partRead(x, mats, pos))))
  }, refname = refname, bamfile = bam, mats = mats, BPPARAM = bpparam)))
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
#' @param diffThreshold report reads which scores differ only
#' below threshold to both haplotypes. Iteratively increases diffThreshold until
#' all positions differ below threshold.
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
scoreHighestSR <- function(srpartition, diffThreshold = 0.001, ...) {

  indent <- list(...)$indent %||% indentation()

  sr <- unique(data.table::as.data.table(srpartition)
               # Get the sum of each read and hptype
               [, clade := sum(prob), by = list(read, haplotype)]
               # get the max of the sums of each
               [, max := max(clade), by = read]
               # dismiss the prob and pos which we dont need anymore
               [, !c("prob")])

  srtmp <- NULL
  sr2 <- NULL
  while (TRUE) {
    flog.info("%sCalculate shortread scoring with cutoff = %s", indent(), diffThreshold, name = "info")
    if (NROW(srtmp) > 0) {
      srtmp <- sr[pos %in% names(which(!correctScoring))]
      srtmp <- srtmp[abs(1 - (clade/max)) < diffThreshold]
    } else {
      srtmp <- sr[abs(1 - (clade/max)) < diffThreshold]
    }
    correctScoring <- vapply(unique(srtmp$pos),
                             function(position, srtmp)
                               .checkSRScoring(position,
                                               srtmp[pos == position]),
                             srtmp = srtmp, FUN.VALUE = logical(1))
    if (all(correctScoring)) {
      if (NROW(sr2) == 0)
        sr2 <- srtmp
      break
    }
    sr2 <- dplyr::bind_rows(sr2, srtmp[pos %in% names(which(correctScoring))])
    diffThreshold <- 1.2*diffThreshold
  }
  dtplyr::tbl_dt(sr2)
}

.writePartFq <- function(fq, fqPart, dontUse = NULL, doUse = NULL, ...) {
  indent <- list(...)$indent %||% indentation()
  fqstream <- ShortRead::FastqStreamer(fq)
  .fileDeleteIfExists(fqPart)
  repeat {
    sr <- ShortRead::yield(fqstream)
    if (length(sr) == 0) break
    fqnames <- as.character(ShortRead::id(sr))
    fqnames <- sub(" .*$", "", fqnames)
    if (is.null(dontUse)) {
      doUse <- which(fqnames %in% doUse)
    } else {
      # useReads = qnames
      doUse <- which(!fqnames %in% dontUse)
    }
    flog.info("%sUsing %s of %s reads", indent(), length(doUse), length(fqnames), name = "info")
    ShortRead::writeFastq(sr[doUse], fqPart, mode = "a", compress = TRUE)
  }
  close(fqstream)
}

####### Helper #############
.partRead <- function(a, mats, pos) {
  pos <- as.character(pos)
  l <- vapply(mats, function(x) x[, pos][as.character(a)], FUN.VALUE = double(1))
  list(prob = l,  haplotype = names(mats), pos = rep(pos, length(l)))
}

.checkSRScoring <- function(position, dfpos){
  scoreOk <- all(vapply(unique(dfpos$haplotype), function(hp, dfpos)
    sum( dfpos$read %in% dfpos[dfpos$haplotype == hp]$read)/NROW(dfpos) > 0.2,
    dfpos = dfpos, FUN.VALUE = logical(1)))
  return(scoreOk)
}
