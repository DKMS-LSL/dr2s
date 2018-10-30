mapInitSR <- function(self, threshold = 0.2, opts = list(), includeDeletions = TRUE,
                      includeInsertions = TRUE, callInsertions = TRUE, clip = FALSE,
                      distributeGaps = FALSE, removeError = TRUE, topx = 0,
                      outdir, force, clean, ...) {

  ## get flags
  forceMapping <- self$getForceMapping()
  microsat     <- self$getMicrosatellite()

  ## Primary mapping
  mapLabel <- "mapInit1"
  reffile  <- self$getRefPath()
  allele   <- self$getReference()
  readtype <- self$getSrdType()
  maptag   <- paste(mapLabel, paste(litArrows(c(allele, readtype,
                                                self$getSrdMapper(),
                                                optstring(opts))),
                                    collapse = " "))

  flog.info(" Map shortreads to provided reference <%s>", self$getReference(), name = "info")
  pileup <- mapReads(
    mapFun = self$getSrMapFun(), maptag = maptag, reffile = reffile,
    allele = allele, readfile = self$getShortreads(), readtype = readtype,
    opts = opts, refname = "", includeDeletions = includeDeletions,
    includeInsertions = includeInsertions, callInsertions = TRUE,
    clip = FALSE, distributeGaps = FALSE, removeError = TRUE,
    topx = 0, outdir = outdir, force = force, clean = clean, ...) #minMapq = 50)

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
  #     readtype = self$getSrdType(),
  #     allele   = self$getReference(),
  #     refname  = "",
  #     force    = force,
  #     outdir   = outdir,
  #     opts     = opts,
  #   )
  # }

  ## Check if the coverage is somewhat equally distributed
  if (max(rowSums(consmat(pileup, freq = FALSE))) /
      quantile(rowSums(consmat(pileup, freq = FALSE)), 0.75) > 5) {
    plotFile <- self$absPath("plot.MapInit.SR.problem.pdf")
    .checkCoverage(pileup, forceMapping, plotFile, maptag)
  }

  ## calc initial consensus
  flog.info(" Construct initial consensus from shortreads", name = "info")
  ##
  baseLabel  <- sub(".bam", "", basename(path(pileup)))
  ## Get conseq
  conseqName <- "Init.consensus.1." %<<% baseLabel
  conseqPath <- file.path(outdir, conseqName %<<% ".fa")
  conseq <- .getWriteConseq(pileup = pileup, name = "mapInit1.0",
                            type = "prob",  threshold = threshold,
                            forceExcludeGaps = TRUE, conseqPath = conseqPath)

  if (microsat) {
    mapLabel <- "mapInit1.2"
    reffile  <- conseqPath
    allele   <- conseqName
    readtype <- self$getSrdType()
    maptag   <- paste(mapLabel, paste0(litArrows(c(conseqName, readtype,
                                                   self$getSrdMapper(),
                                                   optstring(opts))),
                                       collapse = " "))
    flog.info(" Refine microsatellites or repeats by extending the reference", name = "info")
    flog.info(" Remap shortreads to initial consensus from shortreads", name = "info")
    pileup <- mapReads(
      mapFun = self$getSrMapFun(), maptag = maptag, reffile = reffile,
      allele = allele, readfile = self$getShortreads(), readtype = readtype,
      opts = opts, refname = "", includeDeletions = includeDeletions,
      includeInsertions = includeInsertions, callInsertions = TRUE,
      clip = FALSE, distributeGaps = FALSE, removeError = TRUE, topx = 0,
      outdir = outdir, force = force, clean = clean, ...) #minMapq = 50)

    # Infer initial consensus
    flog.info(" Construct second consensus from shortreads " %<<%
                "with refined repeats", name = "info")
    conseqName <- "Init.consensus.2." %<<% baseLabel
    conseqPath <- file.path(outdir, conseqName %<<% ".fa")
    conseq <- .getWriteConseq(pileup, name = "mapInit1.2",
                              type = "prob", threshold = threshold,
                              forceExcludeGaps = TRUE, conseqPath = conseqPath)
  }

  mapInitSR1 = structure(
    list(
      reads   = self$relPath(self$getShortreads()),
      bamfile = self$relPath(path(pileup)),
      pileup  = pileup,
      tag     = maptag,
      conseq  = conseq,
      seqpath = self$relPath(conseqPath),
      ref     = conseqName,
      stats   = list(coverage = .coverage(pileup))
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
                                                 self$getSrdMapper(),
                                                 optstring(opts))),
                                     collapse = " "))

  flog.info(" Remap shortreads to consensus for SNP calling", name = "info")
  pileup <- mapReads(
    mapFun = self$getSrMapFun(), maptag = maptag, reffile = reffile,
    allele = allele, readfile = self$getShortreads(), readtype = readtype,
    opts = opts, refname = "", includeDeletions = TRUE,
    includeInsertions = FALSE, callInsertions = FALSE, clip = FALSE,
    distributeGaps = FALSE, removeError = TRUE, topx = 0,
    outdir = outdir, force = force, clean = clean, ...) #minMapq = 50)

  mapInitSR2 = structure(
    list(
      reads   = self$relPath(self$getShortreads()),
      bamfile = self$relPath(path(pileup)),
      pileup  = pileup,
      tag     = maptag,
      conseq  = conseq,
      seqpath = self$relPath(conseqPath),
      ref     = conseqName,
      stats   = list(coverage = .coverage(pileup))
    ),
    class = c("mapInit", "list")
  )

  list(
    mapInitSR1 = mapInitSR1,
    mapInitSR2 = mapInitSR2
  )
}

