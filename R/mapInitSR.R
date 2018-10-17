mapInitSR <- function(self, threshold = 0.2, opts = list(), optsname = "",
                      minBaseQuality = 3, minMapq = 50, maxDepth = 1e4,
                      minNucleotideDepth = 3, includeDeletions = TRUE,
                      includeInsertions = TRUE, callInsertions = TRUE,
                      clip = FALSE, distributeGaps = FALSE, removeError = TRUE,
                      topx = 0, outdir, force, clean) {

  ## get flags
  forceMapping <- self$getForceMapping()
  microsat     <- self$getMicrosatellite()

  ## Primary mapping
  mapLabel     <- "mapInit1"
  reffile      <- self$getRefPath()
  allele       <- self$getReference()
  readtype     <- self$getSrdType()
  maptag <- paste(mapLabel, paste(litArrows(c(allele, readtype,
                                              self$getSrMapper(),
                                              optstring(opts, optsname))),
                                  collapse = " "))

  flog.info(" Map shortreads to provided reference", name = "info")
  pileup <- mapReads(
    mapFun = self$getSrMapFun(), maptag = maptag, reffile = reffile,
    refseq = NULL, allele = allele, readfile = self$getShortreads(),
    readtype = readtype, threshold = threshold, opts = opts,
    optsname = optsname, refname = "", minBaseQuality = minBaseQuality,
    minMapq = minMapq, maxDepth = maxDepth, minNucleotideDepth = minNucleotideDepth,
    includeDeletions = includeDeletions, includeInsertions = includeInsertions,
    callInsertions = TRUE, clip = FALSE, distributeGaps = FALSE,
    removeError = TRUE, topx = 0, outdir = outdir, force = force, clean = clean)

  # ## TODO: maybe bum this?
  # if (filterScores) {
  #   flog.info(" Filter reads with low alignment scores", name = "info")
  #   ## Run bam - sort - index pipeline
  #   bamfile <- .bamSortIndex(samfile, self$getRefPath(),
  #                             minMapq, force = force, clean = TRUE)
  #   ## Filter Reads
  #   bam <- scanBam(bamfile,
  #                  param = ScanBamParam(tag = "AS",
  #                  what = c("qname", "pos", "cigar")))[[1]]
  #   readfilter <- .filterReads(bam = bam, preserveRefEnds = TRUE)
  #   .fileDeleteIfExists(bamfile)
  #
  #   flog.info(" Write new shortread fastqs to file", name = "info")
  #   fqs <- self$getShortreads()
  #   fqdir <- .dirCreateIfNotExists(file.path(outdir,self$getSrdType()))
  #   # write fastq's
  #   readfile <- c()
  #   readfile <- foreach(fq = fqs, .combine = c) %do% {
  #     srFastqHap = file.path(fqdir, basename(fq))
  #     .writePartFq(fq = fq, srFastqHap = srFastqHap,
  #                   dontUseReads = readfilter)
  #     srFastqHap
  #   }
  #   # set new shortread directory
  #   self$setConfig("filteredShortreads", self$relPath(fqdir))
  #
  #   flog.info(" Map filtered shortreads to provided reference", name = "info")
  #   ## Rerun mapper
  #   flog.info("  Mapping ...", name = "info")
  #   samfile <- mapFun(
  #     reffile  = self$getRefPath(),
  #     readfile = readfile,
  #     allele   = self$getReference(),
  #     readtype = self$getSrdType(),
  #     opts     = opts,
  #     refname  = "",
  #     optsname = optsname,
  #     force    = force,
  #     outdir   = outdir
  #   )
  # }

  ## Check if the coverage is somewhat equally distributed
  if (max(rowSums(pileup$consmat)) /
      quantile(rowSums(pileup$consmat), 0.75) > 5) {
    plotFile <- self$absPath("plot.MapInit.SR.problem.pdf")
    .checkCoverage(pileup, forceMapping, plotFile, maptag)
  }

  ## calc initial consensus
  flog.info(" Construct initial consensus from shortreads", name = "info")
  ## Get conseq
  conseqName <- "Init.consensus." %<<% sub(".bam", "", basename(pileup$bamfile))
  conseqPath <- file.path(outdir, conseqName %<<% ".fa")
  conseq <- .getWriteConseq(pileup = pileup, name = "mapInit1",
                            type = "prob",  threshold = threshold,
                            forceExcludeGaps = TRUE, conseqPath = conseqPath)

  if (microsat) {
    mapLabel <- "mapInit1.2"
    reffile  <- conseqPath
    allele   <- conseqName
    readtype <- self$getSrdType()
    maptag   <- paste(mapLabel, paste0(litArrows(c(conseqName, readtype,
                                                   self$getSrMapper(),
                                                   optstring(opts, optsname))),
                                       collapse = " "))
    flog.info(" Refine microsatellites or repeats by extending the reference",
              name = "info")
    flog.info(" Remap shortreads to initial consensus from shortreads",
              name = "info")
    pileup <- mapReads(
      mapFun = self$getSrMapFun(), maptag = maptag, reffile = reffile,
      refseq = NULL, allele = allele, readfile = self$getShortreads(),
      readtype = readtype, threshold = threshold, opts = opts,
      optsname = optsname, refname = "", minBaseQuality = minBaseQuality,
      minMapq = minMapq, maxDepth = maxDepth, minNucleotideDepth = minNucleotideDepth,
      includeDeletions = includeDeletions, includeInsertions = includeInsertions,
      callInsertions = TRUE, clip = FALSE, distributeGaps = FALSE,
      removeError = TRUE, topx = 0, outdir = outdir, force = force, clean = clean)

    # Infer initial consensus
    flog.info(" Construct second consensus from shortreads " %<<%
                "with refined repeats", name = "info")
    conseqName <- "Init.consensus.2" %<<%
      sub(".bam", "", basename(pileup$bamfile))
    conseqPath <- file.path(outdir, conseqName %<<% ".fa")
    conseq <- .getWriteConseq(pileup, name = "mapInit1.2",
                              type = "prob",  threshold = threshold,
                              forceExcludeGaps = TRUE, conseqPath = conseqPath)
  }

  mapInitSR1 = structure(
    list(
      reads   = self$relPath(self$getShortreads()),
      bamfile = self$relPath(pileup$bamfile),
      pileup  = pileup,
      tag     = maptag,
      conseq  = conseq,
      seqpath = self$relPath(conseqPath),
      ref     = conseqName
    ),
    class  = c("mapInit", "list")
  )

  ## Second mapping to infer polymorphic positions
  ## from same reference as longreads
  mapLabel <- "mapInit2"
  reffile  <- self$absPath(mapInitSR1$seqpath)
  allele   <- mapInitSR1$ref
  readtype <- self$getSrdType()
  maptag   <- paste(mapLabel, paste0(litArrows(c(allele, readtype,
                                                 self$getSrMapper(),
                                                 optstring(opts, optsname))),
                                     collapse = " "))

  flog.info(" Remap shortreads to consensus for SNP calling", name = "info")
  pileup <- mapReads(
    mapFun = self$getSrMapFun(), maptag = maptag, reffile = reffile,
    refseq = NULL, allele = allele, readfile = self$getShortreads(),
    readtype = readtype, threshold = threshold, opts = opts,
    optsname = optsname, refname = "", minBaseQuality = minBaseQuality,
    minMapq = minMapq, maxDepth = maxDepth, minNucleotideDepth = minNucleotideDepth,
    includeDeletions = TRUE, includeInsertions = FALSE, callInsertions = FALSE,
    clip = FALSE, distributeGaps = FALSE, removeError = TRUE, topx = 0,
    outdir = outdir, force = force, clean = clean)

  mapInitSR2 = structure(
    list(
      reads   = self$relPath(self$getShortreads()),
      bamfile = self$relPath(pileup$bamfile),
      pileup  = pileup,
      tag     = maptag,
      conseq  = conseq,
      seqpath = self$relPath(conseqPath),
      ref     = conseqName
    ),
    class = c("mapInit", "list")
  )

  list(
    mapInitSR1 = mapInitSR1,
    mapInitSR2 = mapInitSR2
  )
}
