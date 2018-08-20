.bamSortIndex <- function(samfile,
                          reffile,
                          minMapq = 0,
                          threads = 12,
                          threadmem = "1G",
                          clean = TRUE,
                          force = FALSE) {
  
  samfile <- normalizePath(samfile, mustWork = TRUE)
  reffile <- normalizePath(reffile, mustWork = TRUE)
  ext <- sprintf("%s.sorted",
                 if (minMapq > 0)
                   minMapq %<<% "MAPQ"
                 else
                   "")
  
  samtoolsPath <- Sys.which("samtools")
  sorted <- sub("\\.sam(\\.gz)?", paste(ext, "bam", sep = "."), samfile)
  
  if (nzchar(samtoolsPath)) {
    ## -F260 exclude 'read unmapped', 'not primary alignment'
    ## Better use -F2308 to also exclude chimeric reads!
    tmp <- tempfile()
    view <- sprintf("view -@%s -F2308 -q%s -bT '%s' '%s' > '%s'",
                    threads, minMapq, reffile, samfile, tmp)
    sort <- sprintf("sort -m%s -@%s -o '%s' '%s'",
                    threadmem, threads, sorted, tmp)
    index <- sprintf("index %s", sorted)
    ## Don't execute if file exists and force is false
    if (force || !file.exists(sorted)) {
      system2(samtoolsPath, view)
      system2(samtoolsPath, sort)
      system2(samtoolsPath, index)
    }
  } else {
    ## Sam to bam
    bamfile <- asBam(samfile, tempfile(), indexDestination = TRUE, overwrite = TRUE)
    
    ## Filter the bam by flag
    flag <- scanBamFlag(isUnmappedQuery = FALSE,
                         isSecondaryAlignment = FALSE)
    ## remove also the in Rsamtools unsupported supplmentary alignment flags
    filter <- S4Vectors::FilterRules(list(flag = function(x) x$flag <= flag[2]))
    param <- ScanBamParam(
      flag = flag,
      what = scanBamWhat())
    sorted <- filterBam(file = bamfile, dest = sorted,
                      param = param, filter = filter)
    ## Index the bam
    indexBam(sorted)
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
                         what = NULL) {
  assert_that(
    file.exists(bamfile),
    endsWith(bamfile, ".bam"),
    is.numeric(sampleSize))
  bam <- BamFile(bamfile)
  
  ## Get everything from the bamfile if nothing else is specified
  if(is.null(what))
    what <- scanBamWhat()
  assert_that(is.character(what))
  
  alignmentBam <-  GenomicAlignments::readGAlignments(bam,
                                                      param = ScanBamParam(
                                                        what=scanBamWhat()), 
                                                      use.names = TRUE)
  ## Split the reads if specified. For longreads.
  if (fragmentReads) {
    alignmentBam <- .fragmentReads(alignmentBam, fragmentLength = fragmentWidth)
  }
  ## Use the median read length if no window size is given
  if (is.null(windowSize))
    windowSize <- IRanges::median(GenomicAlignments::qwidth(alignmentBam))
  assert_that(is.numeric(windowSize))
  geneLength <- GenomeInfoDb::seqlengths(GenomeInfoDb::seqinfo(bam))
  
  ## Sample only for reads < 0.5 * geneLength 
  if (0.5*geneLength > windowSize) {
    windows <- seq(from = windowSize, to = geneLength, by = windowSize)
    
    sampledAlignmentBam <- do.call(c, 
                                   lapply(windows, function(i, m, maxCov, 
                                                               windowSize) {
      readMid <- start(m)+floor(windowSize/2)
      m <- m[readMid > i - windowSize & readMid < i]
      if (maxCov <= length(m)) {
        return(sample(m, maxCov))
      } else {
        return(m)
      }
    }, m = alignmentBam, maxCov = sampleSize, windowSize = windowSize))
  } else {
    if (sampleSize <= length(alignmentBam)) {
      sampledAlignmentBam <- sample(alignmentBam, sampleSize)
    } else {
      sampledAlignmentBam <- alignmentBam
    }
  }
  
  newBamfile <- gsub(".bam", ".sampled.bam", bamfile)
  export(sampledAlignmentBam, newBamfile)
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
    lapply(seq_along(alignment),  function(ii, fragmentLength) {
      a <- alignment[ii]
      readWidth <- GenomicAlignments::qwidth(a)
      windowLen <- ceiling(readWidth/fragmentLength)
      wi <- floor(seq(from = 1, readWidth, length.out = windowLen))[2:windowLen]
      b <- do.call(c, lapply(seq_along(wi), function(w, a) {
        start <- ifelse(w == 1, 1, wi[w-1]+1)
        b <- GenomicAlignments::qnarrow(a, start = start, end = wi[w])
        seq <- b@elementMetadata$seq
        seq <- Biostrings::subseq(seq, start = start, end = wi[w])
        b@elementMetadata$seq <- seq
        qual <- b@elementMetadata$qual
        qual <- Biostrings::subseq(qual, start = start, end = wi[w])
        b@elementMetadata$qual <- qual
        b
      }, a = a))
    }, fragmentLength = fragmentLength))
  fragmentAlignment <-unlist(fragmentAlignment, recursive = TRUE, 
                             use.names = TRUE)
  fragmentAlignment[order(GenomicAlignments::start(fragmentAlignment))]
}


