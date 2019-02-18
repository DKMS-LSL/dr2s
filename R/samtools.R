.bamSortIndex <- function(samfile,
                          reffile,
                          minMapq = 0,
                          threads = "auto",
                          threadmem = "1G",
                          clean = TRUE) {
  if (threads == "auto") {
    threads <- .getIdleCores()
  }
  assert_that(is.numeric(threads))
  samfile <- normalizePath(samfile, mustWork = TRUE)
  reffile <- normalizePath(reffile, mustWork = TRUE)
  ext <- sprintf("%s.sorted",
                 if (minMapq > 0)
                   ".MAPQ" %<<% minMapq
                 else
                   "")
  #ext <- ".sorted"
  samtoolsPath <- Sys.which("samtools")
  sorted <- sub("\\.sam(\\.gz)?", dot(c(ext, "bam")), samfile)
  if (nzchar(samtoolsPath)) {
    ## -F260 exclude 'read unmapped', 'not primary alignment'
    ## Better use -F2308 to also exclude chimeric reads!
    tmp <- tempfile()
    view <- sprintf("view -@%s -F2308 -q%s -bT '%s' '%s' > '%s'",
                    threads, minMapq, reffile, samfile, tmp)
    sort <- sprintf("sort -m%s -@%s -o '%s' '%s'",
                    threadmem, threads, sorted, tmp)
    index <- sprintf("index %s", sorted)
    ## Don't execute if file exists
    if (!file.exists(sorted)) {
      system2(samtoolsPath, view)
      system2(samtoolsPath, sort)
      system2(samtoolsPath, index)
    } else {
      warning(sprintf("file <%s> already exists", sorted))
    }
  } else {
    ## Sam to bam
    bamfile <- Rsamtools::asBam(samfile, tempfile(), overwrite = TRUE,
                                indexDestination = TRUE)
    ## Filter the bam by flag
    flag <- Rsamtools::scanBamFlag(isUnmappedQuery = FALSE, isSecondaryAlignment = FALSE)
    ## remove also the in Rsamtools unsupported supplmentary alignment flags
    filter <- S4Vectors::FilterRules(list(flag = function(x) x$flag <= flag[2]))
    param <- Rsamtools::ScanBamParam(
      flag = flag,
      what = Rsamtools::scanBamWhat())
    sorted <- Rsamtools::filterBam(file = bamfile, dest = sorted,
                                   param = param, filter = filter)
    ## Index the bam
    Rsamtools::indexBam(sorted)
  }
  ## Clean up only if the bamfile now exists
  if (clean && file.exists(sorted)) {
    unlink(samfile)
  }
  sorted
}

#' Subsample a bam file according to coverage
#'
#' @param bamfile The bamfile that will be subsampled
#' @param windowSize The window size for subsampling. The number of reads
#'   precised by sampleSize argument will be sampled in this window.
#' @param sampleSize The desired new coverage of the bamfile.
#' @param what What to extract from the bam file to write to the new one.
#'   If NULL defaults to the complete output of
#'   \code{\link[Rsamtools]{scanBamWhat}}.
#' @return A list containing the filename of the original bamfile, the name of
#' the new bamfile, the reference length and the sampleSize
#'
#' @export
subSampleBam <- function(bamfile, windowSize = NULL, sampleSize = 100,
                         fragmentReads = FALSE, fragmentWidth = 1000,
                         what = NULL, clusteredReads = NULL) {
  assert_that(
    file.exists(bamfile),
    endsWith(bamfile, ".bam"),
    is.numeric(sampleSize))

  bam <- Rsamtools::BamFile(bamfile)
  Rsamtools::open.BamFile(bam)
  on.exit(Rsamtools::close.BamFile(bam))

  ## Get everything from the bamfile if nothing else is specified
  if (is.null(what))
    what <- Rsamtools::scanBamWhat()
  assert_that(is.character(what))

  alignmentBam <-  GenomicAlignments::readGAlignments(
    bam, param = Rsamtools::ScanBamParam(what = what), use.names = TRUE)
  ## Use the median read length if no window size is given
  if (is.null(windowSize))
    windowSize <- IRanges::median(GenomicAlignments::qwidth(alignmentBam))
  assert_that(is.numeric(windowSize))
  geneLength <- GenomeInfoDb::seqlengths(GenomeInfoDb::seqinfo(bam))

  ## Sample only for reads < 0.5 * geneLength
  if (0.5*geneLength > windowSize) {
    windows <- seq(from = windowSize, to = geneLength, by = windowSize)
    sampledAlignmentBam <- do.call(c, lapply(windows, function(i, m, maxCov, windowSize) {
      readMid <- GenomicAlignments::start(m) + floor(windowSize/2)
      m <- m[readMid > i - windowSize & readMid < i]
      if (maxCov <= length(m)) {
        return(sample(m, maxCov))
      } else {
        return(m)
      }
    }, m = alignmentBam, maxCov = sampleSize, windowSize = windowSize))
  } else if (length(names(alignmentBam)) > sampleSize) {
    ## Use only longreads of desired lengths, i.e. between .9 and 1.1 of reference length if there are enough
    lens <- GenomicAlignments::qwidth(alignmentBam)
    sampledAlignmentBam <- sample(alignmentBam[lens > 0.9 * geneLength & lens < 1.1 * geneLength], size = sampleSize)
    missingReads <- max(c(sampleSize - length(sampledAlignmentBam), 0))
    sampledAlignmentBam <- c(sampledAlignmentBam,
                             sample(alignmentBam[!names(alignmentBam) %in% names(sampledAlignmentBam)], missingReads))
  } else {
    sampledAlignmentBam <- alignmentBam
  }
  ## Split the reads if specified.
  if (fragmentReads) {
    sampledAlignmentBam <- .fragmentReads(alignment = sampledAlignmentBam, fragmentLength = fragmentWidth)
  }

  # smpld <- mapper$srpartition$A$srpartition$haplotypes
  if (!is.null(clusteredReads)) {
    clustered = ifelse(names(sampledAlignmentBam) %in% clusteredReads, 1, 2)
    S4Vectors::mcols(sampledAlignmentBam)$pt <- clustered
  }
  newBamfile <- gsub(".bam", ".sampled.bam", bamfile)
  rtracklayer::export(sampledAlignmentBam, newBamfile)
  list(
     original        = bamfile,
     sampled         = newBamfile,
     referenceLength = geneLength,
     sampleSize      = sampleSize,
     referenceName   = GenomeInfoDb::seqnames(GenomeInfoDb::seqinfo(bam))
  )
}

.fragmentReads <- function(alignment, fragmentLength = 1000) {
  assert_that(is(alignment, "GAlignments"), is.count(fragmentLength))
  fragmentAlignment <- GenomicAlignments::GAlignmentsList(
    lapply(seq_along(alignment), function(ii, fragmentLength) {
      # ii <- 55
      a <- alignment[ii]
      readWidth <- GenomicAlignments::qwidth(a)
      windowLen <- ceiling(readWidth/fragmentLength)
      wi <- floor(seq(from = 1, readWidth, length.out = windowLen))
      if (length(wi) > 2) {
        wRanges <- IRanges::IRanges(start = c(1, wi[2:(windowLen - 1)] + 1), end = wi[2:windowLen])
        aa <- rep(a, windowLen - 1)
        a <- GenomicAlignments::qnarrow(aa, wRanges)
        S4Vectors::mcols(a)$seq <- GenomicAlignments::narrow(S4Vectors::mcols(a)$seq, wRanges)
        S4Vectors::mcols(a)$qual <- GenomicAlignments::narrow(S4Vectors::mcols(a)$qual, wRanges)
      }
      a
    }, fragmentLength = fragmentLength))
  fragmentAlignment <- unlist(fragmentAlignment, recursive = TRUE, use.names = TRUE)
  fragmentAlignment[order(GenomicAlignments::start(fragmentAlignment))]
}


