## mapping wrapper functions
##
mapReads <- function(
  mapfun, maplabel, reffile, refname, readfile, readtype, opts = NULL, outdir,
  includeDeletions, includeInsertions, callInsertions = FALSE,
  callInsertionThreshold = 0.15, clip = FALSE, distributeGaps = FALSE,
  removeError = TRUE, update_background_model = FALSE, topx = 0, force, clean,
  ...) {

  dots   <- list(...)
  indent <- dots$indent %||% indentation()

  ## Run mapper
  flog.info("%sMap <%s> reads <%s> to reference <%s>", indent(),
            readtype, comma(names(readfile)), names(reffile), name = "info")
  samfile <- mapfun(reffile, refname, readfile, readtype, outdir, maplabel,
                    opts, force)

  ## collect minMapq for use in .bamSortIndex
  # minMapq = 25
  minMapq <- dots$min_mapq %||% dots$minMapq %||% 0

  if (clip) {
    flog.info("%sTrim softclips and polymorphic ends", indent(), name = "info")
    ## Run bam - sort - index pipeline
    bamfile <- .bamSortIndex(samfile = samfile, reffile = reffile,
                             minMapq = minMapq, force = force, clean = clean)
    ## Trim softclips
    fq <- .trimSoftclippedEnds(bam = scanBam(bamfile)[[1]], preserveRefEnds = TRUE)
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
                      opts = list(A = 1, B = 4, O = 2), force)
    # cleanup
    .fileDeleteIfExists(fqout)
    .fileDeleteIfExists(fqdir)
  }

  ## Run bam - sort - index pipeline
  #flog.info("%sSort and index", indent(), name = "info")
  bamfile <- .bamSortIndex(samfile = samfile, reffile = reffile,
                           minMapq = minMapq, force = force)

  ## NOTE: "auto" > 0 evaluates to TRUE
  ## dots <- list()
  if (topx > 0 && readtype != "illumina") {
    pickiness      <- dots$pickiness %||% 1
    lower_limit    <- dots$lower_limit %||% 1
    del_error_rate <- dots$del_error_rate
    reads <- .pickTopXReads(bamfile, topx, pickiness, lower_limit,
                            del_error_rate, indent = incr(indent))
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

  if (update_background_model) {
    delError <- .getGapErrorBackground(consmat(pileup), n = 5)
    flog.info("%sEstimate indel noise <%0.3g> to update background model",
              incr(indent)(), delError, name = "info")
    pileup$meta$del_error_rate <- delError
  }

  if (topx > 0 ) {
    reads(pileup) <- reads
  }

  pileup
}
