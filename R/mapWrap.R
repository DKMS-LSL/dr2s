## mapping wrapper functions
##
mapReads <- function(
  mapfun, maplabel, reffile, refname, readfile, readtype, opts = NULL, outdir,
  includeDeletions, includeInsertions, callInsertions = FALSE,
  callInsertionThreshold = 0.15, clip = FALSE, distributeGaps = FALSE,
  removeError = TRUE, updateBackgroundModel = FALSE, topx = 0, clean, force = FALSE, ...) {

  dots   <- list(...)
  indent <- dots$indent %||% indentation()

  ## Run mapper
  flog.info("%sMap <%s> reads <%s> to reference <%s>", indent(),
            readtype, comma(names(readfile)), names(reffile), name = "info")
  samfile <- mapfun(reffile, refname, readfile, readtype, outdir, maplabel, opts)

  ## collect minMapq for use in .bamSortIndex
  # minMapq = 0
  minMapq <- dots$minMapq %||% dots$min_mapq %||% 0
  if (clip) {
    flog.info("%sTrim softclips and polymorphic ends", indent(), name = "info")
    ## Run bam - sort - index pipeline
    bamfile <- .bamSortIndex(samfile = samfile, reffile = reffile,
                             minMapq = minMapq, clean = clean)
    ## Trim softclips
    fq <- .trimSoftclippedEnds(bam = Rsamtools::scanBam(bamfile)[[1]], preserveRefEnds = TRUE)
    ## Trim polymorphic ends
    fq <- .trimPolymorphicEnds(fq)
    ## Write new shortread file to disc
    fqdir  <- .dirCreateIfNotExists(file.path(outdir, "trimmed"))
    fqfile <- gsub(".fastq(.gz)?", ".trimmed.fastq.gz", basename(readfile[1]))
    readfile <- .fileDeleteIfExists(file.path(fqdir, fqfile))
    ShortRead::writeFastq(fq, readfile, compress = TRUE)
    .fileDeleteIfExists(bamfile)
    ## Rerun mapper
    flog.info("%sRemap trimmed reads <%s> to reference <%s>", indent(), fqfile, names(reffile), name = "info")
    samfile <- mapfun(reffile, refname, readfile, readtype, outdir, maplabel,
                      opts = list(A = 1, B = 4, O = 2))
    # cleanup
    .fileDeleteIfExists(readfile)
    .fileDeleteIfExists(fqdir)
  }

  ## Run bam - sort - index pipeline
  #flog.info("%sSort and index", indent(), name = "info")
  bamfile <- .bamSortIndex(samfile = samfile, reffile = reffile, minMapq = minMapq, force = force)

  ## NOTE: "auto" > 0 evaluates to TRUE, FALSE > 0 evaluates to FALSE
  ## dots <- list()
  ## indelRate <- NULL
  if (topx > 0 && readtype != "illumina") {
    pickiness  <- dots$pickiness  %||% 1
    lowerLimit <- dots$lowerLimit %||% 20
    indelRate  <- dots$indelRate
    reads <- .pickTopXReads(bamfile, topx, pickiness, lowerLimit, indelRate, indent = incr(indent))
    alignmentBam <- GenomicAlignments::readGAlignments(
      file = Rsamtools::BamFile(bamfile),
      param = Rsamtools::ScanBamParam(
        what = Rsamtools::scanBamWhat(),
        which = IRanges::IRangesList()),
      use.names = TRUE)
    if (length(reads) < length(alignmentBam))
      rtracklayer::export(alignmentBam[reads], bamfile)
  }

  ## Calculate pileup from graphmap produced SAM file
  ## pParam = .collectPileupParams(includeDeletion = includeDeletions, includeInsertions = includeInsertions)
  ## pileup <- pileup(bamfile, reffile, readtype, pParam = pParam)
  #flog.info("%sPile up", indent(), name = "info")
  pileup <- pileup(bamfile, reffile, readtype, indent = incr(indent),
                   pParam = .collectPileupParams(
                     includeDeletion = includeDeletions,
                     includeInsertions = includeInsertions,
                     ...))

  if (distributeGaps) {
    #flog.info("%sDistribute gaps", indent(), name = "info")
    consmat(pileup) <- .distributeGaps(mat = consmat(pileup),
                                       bamfile = bampath(pileup),
                                       removeError = removeError,
                                       indent = incr(indent))
  }

  if (callInsertions && is.null(ins(consmat(pileup)))) {
    #flog.info("%sCall insertions", indent(), name = "info")
    ## TODO check threshold
    pileup <- .pileupIncludeInsertions(x = pileup,
                                       threshold = callInsertionThreshold,
                                       indent = incr(indent))
  }

  if (updateBackgroundModel) {
    noise <- .getGapErrorBackground(consmat(pileup), n = 5)
    flog.info("%sEstimate indel noise <%0.3g> to update background model",
              incr(indent)(), noise, name = "info")
    indelRate(pileup) <- noise
  }

  if (topx > 0 ) {
    reads(pileup) <- reads
  }

  pileup
}
