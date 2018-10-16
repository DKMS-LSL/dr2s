## mapping wrapper functions
##
mapReads <- function(
  maptag, reffile, readfile, allele, readtype, opts = NULL,  refname = "",
  optsname = "", force, outdir, minMapq = 0, clean,  threshold, maxDepth = 1e4,
  minBaseQuality = 3,  minNucleotideDepth = 3, refseq = NULL, includeInsertions,
  mapFun,   clip = FALSE, distributeGaps = FALSE, includeDeletions, 
  callInsertions = FALSE, removeError = TRUE, topx = 0) {
  
  ## Run mapper
  flog.info("  Mapping ...", name = "info")
  samfile <- mapFun(reffile = reffile, readfile = readfile, allele = allele,
    readtype = readtype, opts = opts, refname  = refname, optsname = optsname,
    force = force, outdir = outdir)

  if (clip) {
    flog.info("  Trimming softclips and polymorphic ends ...",
              name = "info")
    ## Run bam - sort - index pipeline
    bamfile <- .bamSortIndex(samfile = samfile, reffile = reffile,
                             minMapq = minMapq, force = force, clean = TRUE)
    ## Trim softclips
    fq <- .trimSoftclippedEnds(bam = scanBam(bamfile)[[1]],
                               preserveRefEnds = TRUE)
    ## Trim polymorphic ends
    fq <- .trimPolymorphicEnds(fq)
    ## Write new shortread file to disc
    fqdir  <- .dirCreateIfNotExists(file.path(outdir,
                                              "mapFinal"))
    fqfile <- gsub(".fastq(.gz)?", ".trimmed.fastq", basename(readfile[1]))
    fqout  <- .fileDeleteIfExists(file.path(fqdir, fqfile))
    ShortRead::writeFastq(fq, fqout, compress = TRUE)
    .fileDeleteIfExists(bamfile)
    ## Rerun mapper
    flog.info("   Mapping trimmed short reads against latest consensus ...",
              name = "info")
    samfile <- mapFun(
      reffile = reffile, readfile = fqout, allele = maptag, readtype = readtype,
      opts = list(A = 1, B = 4, O = 2), refname = refname,
      optsname = optstring(opts), force = force, outdir = outdir)
    # cleanup
    .fileDeleteIfExists(fqout)
  }

  ## Run bam - sort - index pipeline
  flog.info("  Indexing ...", name = "info")
  bamfile <- .bamSortIndex(samfile = samfile, reffile = reffile,
                           minMapq = minMapq, force = force)

  if (topx > 0) {
    flog.info("  Extracting the top %s reads ...", topx, name = "info")
    reads <- .topXReads(bamfile, refseq, n = topx)
    bam <- BamFile(bamfile)
    alignmentBam <-  GenomicAlignments::readGAlignments(bam,
                                                        param = ScanBamParam(
                                                          what=scanBamWhat()), 
                                                        use.names = TRUE)
    rtracklayer::export(alignmentBam[reads], bamfile)
  }
  ## Calculate pileup from graphmap produced SAM file
  flog.info("  Piling up ...", name = "info")
  pileup <- Pileup(bamfile, threshold, max_depth = maxDepth,
    min_base_quality = minBaseQuality, min_mapq = minMapq,
    min_nucleotide_depth = minNucleotideDepth,
    include_deletions = includeDeletions,
    include_insertions = includeInsertions)
  

  if (distributeGaps) {
    flog.info("  Distributing gaps ...", name = "info")
    pileup$consmat <- .distributeGaps(mat = pileup$consmat,
                                      bamfile = bamfile,
                                      reference = refseq,
                                      removeError = removeError)
  }

  if (callInsertions && is.null(ins(pileup$consmat))) {
    flog.info("  Calling insertions ...", name = "info")
    ## TODO check threshold
    pileup <- .pileupIncludeInsertions(x = pileup, threshold = 0.15, readtype = readtype)
  }
  if (topx > 0) 
    pileup$reads <- reads
  pileup
}
