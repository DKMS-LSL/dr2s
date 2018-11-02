## mapping wrapper functions
##
mapReads <- function(
  mapFun, maptag, reffile, allele, readfile, readtype, opts = NULL,
  refname = "", includeDeletions, includeInsertions, callInsertions = FALSE,
  clip = FALSE, distributeGaps = FALSE, removeError = TRUE, topx = 0,
  outdir, force, clean, ...) {

  indent <- list(...)$indent %||% indentation()
  flog.info("%sMap <%s> reads <%s> to reference <%s>", indent(),
            readtype, comma(names(readfile)), names(reffile), name = "info")

  ## Run mapper
  samfile <- mapFun(reffile, readfile, readtype, allele, refname, force, outdir, opts)
  ## collect minMapq for use in .bamSortIndex
  # minMapq = 0
  minMapq <- list(...)$min_mapq %||% list(...)$minMapq %||% 0

  if (clip) {
    flog.info("%sTrim softclips and polymorphic ends", indent(), name = "info")
    ## Run bam - sort - index pipeline
    bamfile <- .bamSortIndex(samfile = samfile, reffile = reffile,
                             minMapq = minMapq, force = force, clean = TRUE)
    ## Trim softclips
    fq <- .trimSoftclippedEnds(bam = scanBam(bamfile)[[1]], preserveRefEnds = TRUE)
    ## Trim polymorphic ends
    fq <- .trimPolymorphicEnds(fq)
    ## Write new shortread file to disc
    fqdir  <- .dirCreateIfNotExists(file.path(outdir, "trimmed"))
    fqfile <- gsub(".fastq(.gz)?", ".trimmed.fastq", basename(readfile[1]))
    fqout  <- .fileDeleteIfExists(file.path(fqdir, fqfile))
    ShortRead::writeFastq(fq, fqout, compress = TRUE)
    .fileDeleteIfExists(bamfile)
    ## Rerun mapper
    flog.info("%sRemap trimmed reads <%s> to reference <%s>", indent(), fqfile, names(reffile), name = "info")
    samfile <- mapFun(reffile, fqout, readtype, allele, refname, force, outdir, opts = list(A = 1, B = 4, O = 2))
    # cleanup
    .fileDeleteIfExists(fqout)
    .fileDeleteIfExists(fqdir)
  }

  ## Run bam - sort - index pipeline
  flog.info("%sSort and index", indent(), name = "info")
  bamfile <- .bamSortIndex(samfile = samfile, reffile = reffile,
                           minMapq = minMapq, force = force)

  if (topx > 0 && readtype != "illumina") {
    flog.info("%sExtract the top-scoring %s reads", indent(), topx, name = "info")
    reads <- .topXReads(bamfile, n = topx, indent = incr(indent))
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
  flog.info("%sPile up", indent(), name = "info")
  pileup <- pileup(bamfile, reffile, readtype, indent = incr(indent),
                   pParam = .collectPileupParams(
                     includeDeletion = includeDeletions,
                     includeInsertions = includeInsertions,
                     ...))
  if (distributeGaps) {
    flog.info("%sDistribute gaps", indent(), name = "info")
    consmat(pileup) <- .distributeGaps(mat = consmat(pileup),
                                       bamfile = path(pileup),
                                       removeError = removeError,
                                       indent = incr(indent))
  }

  if (callInsertions && is.null(ins(consmat(pileup)))) {
    flog.info("%sCall insertions", indent(), name = "info")
    ## TODO check threshold
    pileup <- .pileupIncludeInsertions(x = pileup,
                                       threshold = 0.15,
                                       indent = incr(indent))
  }

  if (topx > 0) {
    reads(pileup) <- reads
  }

  pileup
}
