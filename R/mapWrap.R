## mapping wrapper functions
##
mapReads <- function(
  mapFun, maptag, reffile, allele, readfile, readtype, opts = NULL,
  refname = "", includeDeletions, includeInsertions, callInsertions = FALSE,
  clip = FALSE, distributeGaps = FALSE, removeError = TRUE, topx = 0,
  outdir, force, clean, ...) {

  ## Run mapper
  flog.info("  Mapping to <%s> ...", names(reffile), name = "info")
  samfile <- mapFun(reffile, readfile, readtype, allele, refname, force, outdir, opts)

  ## collect minMapq for use in .bamSortIndex
  # dots <- list(minMapq = 0)
  dots <- list(...)
  minMapq <- dots$min_mapq %||% dots$minMapq %||% 0

  if (clip) {
    flog.info("  Trimming softclips and polymorphic ends ...", name = "info")
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
    flog.info("   Mapping trimmed short reads against latest consensus ...", name = "info")
    samfile <- mapFun(reffile, readfile, readtype, allele, refname, force,
                      outdir, opts = list(A = 1, B = 4, O = 2))
    # cleanup
    .fileDeleteIfExists(fqout)
    .fileDeleteIfExists(fqdir)
  }

  ## Run bam - sort - index pipeline
  flog.info("  Indexing ...", name = "info")
  bamfile <- .bamSortIndex(samfile = samfile, reffile = reffile,
                           minMapq = minMapq, force = force)

  if (topx > 0 && readtype != "illumina") {
    flog.info("  Extracting the top-scoring %s longreads ...", topx, name = "info")
    reads <- .topXReads(bamfile, n = topx)
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
  flog.info("  Piling up ...", name = "info")
  pileup <- pileup(bamfile, reffile, readtype,
                   pParam = .collectPileupParams(
                     includeDeletion = includeDeletions,
                     includeInsertions = includeInsertions,
                     ...))
  if (distributeGaps) {
    flog.info("  Distributing gaps ...", name = "info")
    consmat(pileup) <- .distributeGaps(mat = consmat(pileup),
                                       bamfile = path(pileup),
                                       removeError = removeError)
  }

  if (callInsertions && is.null(ins(consmat(pileup)))) {
    flog.info("  Calling insertions ...", name = "info")
    ## TODO check threshold
    pileup <- .pileupIncludeInsertions(x = pileup, threshold = 0.15)
  }

  if (topx > 0) {
    reads(pileup) <- reads
  }

  pileup
}
