# Method: mapInit ####
#' @export
mapInit.DR2S <- function(x,
                         opts = list(),
                         threshold = NULL,
                         includeDeletions = TRUE,
                         includeInsertions = TRUE,
                         microsatellite = FALSE,
                         filterScores = TRUE,
                         forceMapping = FALSE,
                         topx = 0,
                         createIgv = TRUE,
                         force = FALSE,
                         plot = TRUE,
                         ...) {
  x$runMapInit(opts = opts,
               threshold = threshold,
               includeDeletions = includeDeletions,
               includeInsertions = includeInsertions,
               microsatellite = microsatellite,
               filterScores = filterScores,
               forceMapping = forceMapping,
               topx = topx,
               createIgv = createIgv,
               force = force,
               plot = plot,
               ...)
  invisible(x)
}

DR2S_$set("public", "runMapInit", function(opts = list(),
                                           threshold = NULL,
                                           includeDeletions = TRUE,
                                           includeInsertions = TRUE,
                                           microsatellite = FALSE,
                                           filterScores = TRUE,
                                           forceMapping = FALSE,
                                           topx = 0,
                                           createIgv = TRUE,
                                           force = FALSE,
                                           plot = TRUE,
                                           ...) {
  # # debug
  # opts = list()
  # threshold = NULL
  # includeDeletions = TRUE
  # includeInsertions = TRUE
  # microsatellite = FALSE
  # filterScores = TRUE
  # forceMapping = FALSE
  # topx = 0
  # createIgv = TRUE
  # force = FALSE
  # plot = TRUE
  # library(assertthat)
  # library(ggplot2)
  # library(S4Vectors)
  # library(Rsamtools)
  # library(foreach)
  # library(futile.logger)
  # library(cowplot)
  # self <-dr2s
  # self <- mapper

  flog.info("# mapInit", name = "info")

  ## Collect starttime for mapInit runstats
  start.time <- Sys.time()
  ## Initiate indenter
  indent  <- indentation(1)
  indent2 <- incr(indent)
  ## Overide default arguments
  args <- self$getOpts("mapInit")
  if (!is.null(args)) {
    env  <- environment()
    list2env(args, envir = env)
  }
  if (is.null(threshold)) {
    threshold <- self$getThreshold()
  }
  if (!exists("callInsertionThreshold")) {
    ## an insertion needs to be at frequency <callInsertionThreshold> for it
    ## to be included in the pileup.
    callInsertionThreshold <- 0.2
  }
  if (!exists("minMapq")) {
    ## don't filter for mapping quality unless specified.
    ## for <shortreads> we hardcode <minMapq = 50>
    minMapq <- 0
  }

  ## Get options and prepare mapping
  outdir <- .dirCreateIfNotExists(self$absPath("mapInit"))
  .dirCreateIfNotExists(file.path(self$absPath(".plots")))
  clean <- TRUE
  igv <- list()
  SR <- list()

  if (self$hasShortreads()) {
    SR <- mapInitSR(
      self = self, threshold = threshold, opts = opts, includeDeletions = includeDeletions,
      includeInsertions = includeInsertions, callInsertions = TRUE,
      callInsertionThreshold = callInsertionThreshold, clip = FALSE,
      distributeGaps = FALSE, removeError = TRUE, topx = 0, outdir = outdir,
      force = force, clean = clean, minMapq = 50, indent = indent, ...)

    ### TODO wrap this command up
    if (createIgv)
      igv[["SR"]] <- createIgvJsFiles(
        reference = self$absPath(refpath(SR$SR2)),
        bamfile = self$absPath(bampath(SR$SR2)),
        outdir = outdir,
        sampleSize = 100)

    reffile  <- self$absPath(conspath(SR$SR1))
    refname  <- consname(SR$SR1)
  }
  else {
    maplabel <- "mapInit1"
    mapfun   <- self$getLrdMapFun()
    readtype <- self$getLrdType()
    readfile <- self$getLongreads()
    reffile  <- self$getRefPath()
    refname  <- self$getReference()
    pileup   <- mapReads(
      mapfun = mapfun, maplabel = maplabel, reffile = reffile, refname = refname,
      readfile = readfile, readtype = readtype, opts = opts, outdir = outdir,
      includeDeletions = TRUE, includeInsertions = TRUE, callInsertions = TRUE,
      callInsertionThreshold = callInsertionThreshold, clip = FALSE,
      distributeGaps = TRUE, removeError = TRUE, topx = topx, force = force,
      clean = clean, minMapq = minMapq, indent = indent, ...)

    if (!is.null(reads(pileup))) {
      fqfile <- paste(
        self$getSampleId(), self$getLrdType(), paste0("n", length(reads(pileup))),
        "fastq", "gz", sep = ".")
      topx_fqpath <- .fileDeleteIfExists(file.path(outdir, strip(fqfile, "_")))
      names(topx_fqpath) <- self$relPath(topx_fqpath)
      fq <- .extractFastq(bampath(pileup), reads(pileup))
      ShortRead::writeFastq(fq, topx_fqpath, compress = TRUE)
    }

    ## Construct initial longread consensus sequence
    refname <- refname %<<% ".consensus"
    self$setRefPath(file.path(outdir, strip(maplabel %<<% "." %<<% refname %<<% ".fa", "_")))
    reffile <- self$getRefPath()
    flog.info("%sConstruct consensus <%s>", indent2(), names(reffile), name = "info")
    conseq <- .writeConseq(x = pileup, name = refname, type = "prob",
                           threshold = threshold, suppressAllGaps = TRUE,
                           replaceIndel = "N", conspath = reffile)
  }

  maplabel <- "mapInit2"
  mapfun   <- self$getLrdMapFun()
  readtype <- self$getLrdType()
  readfile <- if (exists("topx_fqpath")) topx_fqpath else self$getLongreads()
  pileup <- mapReads(
    mapfun = mapfun, maplabel = maplabel, reffile = reffile, refname = refname,
    readfile = readfile, readtype = readtype, opts = opts, outdir = outdir,
    includeDeletions = TRUE, includeInsertions = FALSE, callInsertions = FALSE,
    clip = FALSE, distributeGaps = TRUE, removeError = TRUE, topx = topx,
    force = force, clean = clean, minMapq = minMapq, indent = indent, ...)

  if (createIgv)
    igv[["LR"]] <- createIgvJsFiles(
      reference = refpath(pileup),
      bamfile = bampath(pileup),
      outdir = outdir,
      sampleSize = 100,
      fragmentReads = TRUE)

  self$mapInit = MapList_(
    ## mapdata
    readpath  = self$relPath(readfile),
    refpath   = self$relPath(reffile),
    bampath   = self$relPath(bampath(pileup)),
    conspath  = NULL,
    pileup    = pileup,
    stats     = list(coverage = .coverage(pileup)),
    ## metadata
    maplabel  = maplabel,
    refname   = refname,
    mapper    = self$getLrdMapper(),
    opts      = opts,
    ## additional metadata
    SR1       = SR$SR1,
    SR2       = SR$SR2,
    igv       = igv)

  createIgvConfigs(x = self, map = "mapInit", open = "FALSE")

  if (plot) {
    flog.info("%sPlot MapInit summary", indent(), name = "info")
    ## Coverage and frequency of minor alleles
    p <- self$plotMapInitSummary(
      thin = 0.25,
      width = 2
    )
    plotRows <- ifelse(self$hasShortreads(), 2, 1)
    cowplot::save_plot(self$absPath("plot.MapInit.pdf"),
                       plot = p, ncol = 1, nrow = plotRows,
                       base_aspect_ratio = as.numeric(paste(5, plotRows, sep = ".")),
                       title = paste(self$getLocus(), self$getSampleId(), sep = "." ))
    cowplot::save_plot(self$absPath(".plots/plot.MapInit.svg"),
                       plot = p, ncol = 1, nrow = plotRows,
                       base_aspect_ratio = as.numeric(paste(5, plotRows, sep = ".")))
  }

  ## set mapInit runstats
  if (is(meta(self$mapInit, "SR2"), "MapList")) {
    .setRunstats(self, "mapInit",
                 list(Runtime = format(Sys.time() - start.time),
                      SRcoverage = stats(meta(self$mapInit, "SR2"), "coverage")[["50%"]],
                      LRcoverage = stats(self$mapInit, "coverage")[["50%"]]))
  } else {
    .setRunstats(self, "mapInit",
                 list(Runtime = format(Sys.time() - start.time),
                      LRcoverage = stats(self$mapInit, "coverage")[["50%"]]))
  }

  return(invisible(self))
})

# Method: partitionLongreads ####
#' @export
partitionLongreads.DR2S <- function(x,
                                    threshold         = NULL,
                                    skipGapFreq       = 2/3,
                                    distAlleles       = NULL,
                                    noGapPartitioning = FALSE,
                                    selectAllelesBy   = "count",
                                    minClusterSize    = 15,
                                    plot              = TRUE,
                                    ...) {
  ## Collect start time for partitionLongreads runstats
  start.time <- Sys.time()

  x$runPartitionLongreads(threshold         = threshold,
                          skipGapFreq       = skipGapFreq,
                          noGapPartitioning = noGapPartitioning,
                          selectAllelesBy   = selectAllelesBy,
                          minClusterSize    = minClusterSize,
                          distAlleles       = distAlleles,
                          plot              = plot,
                          ...)
  x$runSplitLongreadsByHaplotype(plot = plot)
  x$runExtractPartitionedLongreads()

  ## set partitionLongreads runstats
  .setRunstats(x, "partitionLongreads",
               list(Runtime = format(Sys.time() - start.time)))

  return(invisible(x))
}

DR2S_$set("public",
          "runPartitionLongreads",
          function(threshold = NULL,
                   skipGapFreq = 2/3,
                   distAlleles = NULL,
                   noGapPartitioning = FALSE,
                   selectAllelesBy = "count",
                   minClusterSize = 15,
                   plot = TRUE,
                   ...) {
    ## debug
    # threshold = NULL
    # skipGapFreq = 2/3
    # distAlleles = NULL
    # noGapPartitioning = TRUE
    # selectAllelesBy = "distance"
    # minClusterSize = 15
    # plot = TRUE
    # library(futile.logger)
    # library(assertthat)
    # self <- dr2s

    flog.info("# PartitionLongreads", name = "info")

    ## Initiate indenter
    indent <- indentation(1)
    ## Overide default arguments
    args <- self$getOpts("partitionLongreads")
    if (!is.null(args)) {
      env  <- environment()
      list2env(args, envir = env)
    }
    if (is.null(threshold)) {
      threshold <- self$getThreshold()
    }
    if (is.null(distAlleles)) {
      distAlleles <- self$getDistAlleles()
    }
    assert_that(
      self$hasMapInit(),
      is.double(skipGapFreq),
      is.double(threshold),
      is.count(distAlleles),
      is.logical(plot)
    )

    ## Get the reference sequence
    if (is(meta(self$mapInit, "SR1"), "MapList")) {
      useSR  <- TRUE
      flog.info("%sConstruct SNP matrix from shortreads", indent(), name = "info")
    } else {
      useSR  <- FALSE
      flog.info("%sConstruct SNP matrix from longreads", indent(), name = "info")
    }
    ppos <- self$polymorphicPositions(useSR = useSR)

    ## Spurious gaps, especially in longreads can hinder a correct clustering
    ## Remove gap positions for clustering
    if (noGapPartitioning) {
      allpp <- NROW(ppos)
      ppos <- dplyr::filter(ppos, a1 != "-" & a2 != "-")
      flog.info("%sUse %s non-indel polymorphisms out of %s for clustering", indent(),
                NROW(ppos), allpp, name = "info")
    }

    ## Check if already finished because it is a homozygous sample
    if (NROW(ppos) == 0) {
      flog.warn("%sNo polymorphic positions found for clustering", indent(), name = "info")
      flog.info("%sEntering single allele polish and report pipeline", indent(), name = "info")
      ## set runstats
      .setRunstats(self, "partitionLongreads",
                   list(foundPolymorphicPositions = 0L))
      return(invisible(.finishCn1(x = self, plot = plot)))
    }

    mat <- SNPmatrix(self$absPath(bampath(self$mapInit)), ppos)

    flog.info("%sPartition %s longreads over %s SNPs", indent(), NROW(mat), NCOL(mat), name = "info")
    prt <- partitionReads(x = mat,
                          skipGapFreq = skipGapFreq,
                          deepSplit = 1,
                          threshold = threshold,
                          distAlleles = distAlleles,
                          sortBy = selectAllelesBy,
                          minClusterSize = minClusterSize,
                          indent = incr(indent))

    ## set runstats
    .setRunstats(self, "partitionLongreads",
                 list(nLongreads = NROW(mat),
                      foundPolymorphicPositions = NCOL(mat),
                      usedPolymorphicPositions = K(prt),
                      foundClades = as.list(table(PRT(prt)))))

    ## Set sample haplotypes
    self$setHapTypes(levels(as.factor(PRT(prt))))

    # Check if we have only one cluster and finish the pipeline if so
    if (length(self$getHapTypes()) == 1) {
      flog.warn("%sOnly one allele left", indent(), name = "info")
      flog.info("%sEntering single allele polish and report pipeline", indent(), name = "info")
      return(invisible(.finishCn1(x = self, plot = plot)))
    }

    self$lrpartition = structure(list(
      mat = mat,
      prt = prt,
      hpl = NULL,
      lmt = NULL
    ),
    class = c("PartList", "list"))

    return(invisible(self))
  })

#' @export
print.PartList <- function(x, ...) {
  msg <- sprintf("An object of class '%s'.\n", class(x)[1])
  msg <- sprintf("%s [Matrix]    %s reads; %s polymorphic positions\n",
                 msg, NROW(x$mat), NCOL(x$mat))
  cat(msg)
  cat(" [Partition] ")
  print(x$prt, ncols = 4, nrows = 6)
  cat(" [Haplotypes] ")
  print(x$hpl)
}

DR2S_$set("public", "runSplitLongreadsByHaplotype", function(plot = TRUE) {

  ## Initiate indenter
  indent <- indentation(1)

  flog.info("%sSplit partitioned longreads by score", indent(), name = "info")
  ## Check if reporting is already finished and exit safely
  if (.checkReportStatus(self)) return(invisible(self))
  assert_that(self$hasLrdPartition())

  ## Overide default arguments
  args <- self$getOpts("partitionLongreads")
  if (!is.null(args)) {
    env  <- environment()
    list2env(args, envir = env)
  }

  prt <- partition(self$getPartition())
  haplotypes <- levels(prt$haplotype)
  tag <- self$getMapTag("init", "LR")

  # Set all limits to NULL
  self$setLimits(rep(NA, length(haplotypes)))
  prts <- lapply(haplotypes, function(x) prt[prt$haplotype == x,])
  names(prts) <- haplotypes
  scores <- lapply(prts, function(x) x$mcoef)
  lmts <- .optimalPartitionLimits(scores)
  self$setLimits(setNames(as.list(lmts$limits$c), haplotypes))
  self$lrpartition$lmt <- lmts$plt

  # Get only reads within the limit
  reads <- setNames(lapply(names(self$getLimits()), function(x) {
    dplyr::filter(prt, haplotype == x, mcoef >= self$getLimits()[x])
  }), names(self$getLimits()))

  ## Initiate indenter
  indent2 <- incr(indent)
  for (hp in haplotypes) {
    flog.info("%s<%s>: Using %s longreads with score > %.2f",
              indent2(), hp, nrow(reads[[hp]]), self$getLimits()[hp],
              name = "info")
  }

  ## set runstats
  .setRunstats(self, "partitionLongreads",
               list(reads = setNames(as.list(as.integer(lmts$limits$r)), haplotypes)))

  # results structure
  self$lrpartition$hpl <- structure(
    setNames(lapply(reads, function(x) {
      structure(
        x$read,
        q = x$mcoef,
        freq = NROW(x)/NROW(dplyr::bind_rows(reads)),
        limit = self$getLimits()[[x$haplotype[1]]]
      )}), haplotypes),
    class = c("HapList", "list"))

  if (plot) {
    ## Plot the consensus sequences from clustering
    .browseSeqs(SQS(self$getPartition()),
                file = self$absPath("partition.fa.html"),
                openURL = FALSE)

    p <- self$plotPartitionSummary(
      label = tag, limits = unlist(self$getLimits()))

    cowplot::save_plot(self$absPath("plot.Partition.pdf"), plot = p,
                       title = paste(self$getLocus(), self$getSampleId(), sep = "."),
                       base_height = 12, base_width = 10)
    cowplot::save_plot(self$absPath(".plots/plot.Partition.svg"), plot = p,
                       base_aspect_ratio = 1.2)

    outf  <- self$absPath("plot.Sequence")
    ppos <- SNP(self$getPartition())
    names(ppos) <- seq_along(ppos)
    pwm <- lapply(PWM(self$getPartition()), function(pwm) {
      pwm[pwm < 0.1] <- 0
      pwm
    })

    p <- suppressMessages(self$plotSeqLogo(ppos, pwm))
    cowplot::save_plot(filename = self$absPath("plot.Sequence.pdf"),
                       plot        = p,
                       base_width  = 0.4*length(ppos) + 1.4,
                       base_height = 2.5*length(pwm),
                       title       = paste(self$getLocus(), self$getSampleId(),
                                           sep = "." ),
                       units       = "cm",
                       limitsize   = FALSE)
    cowplot::save_plot(filename = self$absPath(".plots/plot.Sequence.svg"),
                       plot        = p,
                       base_width  = 0.4*length(ppos) + 1.4,
                       base_height = 2.5*length(pwm),
                       units       = "cm",
                       limitsize   = FALSE)

  }

  return(invisible(self))
})

#' @export
print.HapList <- function(x, ...) {
  msg <- sprintf("An object of class '%s'.\n", class(x)[1])
  for (haplotype in names(x)) {
    msg <- msg %<<%
      sprintf("%s: n %s; frequency %s; limit %s\n",
              haplotype,
              length(x[haplotype]),
              round(attr(x[[haplotype]], "freq"), 3),
              attr(x[[haplotype]], "limit"))
  }
  cat(msg)
}

DR2S_$set("public", "runExtractPartitionedLongreads", function() {
  ## Check if reporting has already been done and exit safely
  if (.checkReportStatus(self)) return(invisible(self))
  assert_that(self$hasHapList())

  ## Initiate indenter
  indent <- indentation(1)

  ## Overide default arguments
  args <- self$getOpts("partitionLongreads")
  if (!is.null(args)) {
    env  <- environment()
    list2env(args, envir = env)
  }

  ## extract msa from mapInit
  bamfile <- self$absPath(bampath(self$mapInit))
  mat <- .msaFromBam(Rsamtools::BamFile(bamfile), paddingLetter = ".")
  flog.info("%sExtract partitioned longreads and construct consensus" %<<%
            " sequences based on <%s>", indent(), names(bamfile), name = "info")

  ## For each <hptype>:
  ##    create the "mapIter/<hptype>" directory to hold the reads assigned to <hptype>
  ##    extract the reads assigned to <hptype> from the LR bam file from <mapInit>
  ##    construct a consensus sequence from the reads for each <hptype>
  ##    create a <mapIter0> object for each <hptype>
  ##    hptype <- "A"
  ##    hptype <- "B
  indent2 <- incr(indent)
  for (hptype in self$getHapTypes()) {
    flog.info("%sFor haplotype <%s>:", indent(), hptype, name = "info")
    ## create the "mapIter/<hptype>" directory
    fqdir <- .dirCreateIfNotExists(
      normalizePath(file.path(self$getOutdir(), "mapIter", hptype), mustWork = FALSE))
    ## get the ids of the reads assigned to <hptype>
    readIds  <- self$getHapList(hptype)
    ## extract these reads and write to file
    fq <- .extractFastq(bamfile, qnames = readIds)
    fqfile <- paste(
      "hap" %<<% hptype, self$getLrdType(), self$getLrdMapper(),
      "lim" %<<% as.character(floor(self$getLimits()[[hptype]])),
      "n" %<<% length(fq), "fastq", "gz", sep = ".")
    fqout <- .fileDeleteIfExists(file.path(fqdir, fqfile))
    nfq <- ShortRead::writeFastq(fq, file = fqout, full = FALSE, compress = TRUE)
    ## assert that all records got written to file
    assert_that(nfq == length(fq))
    ## Extract consensus matrix from mapInit for the clustered reads
    cmat <- .extractIdsFromMat(mat, readIds)
    ## Construct the initial longread consensus sequence for <hptype>
    maplabel <- "mapIter0"
    consname <- maplabel %<<% ".consensus." %<<% hptype
    conspath <- file.path(fqdir, consname %<<% ".fa")
    names(conspath) <- self$relPath(conspath)
    flog.info("%sConstruct consensus <%s>", indent2(), names(conspath), name = "info")
    conseq <- .writeConseq(x = cmat, name = consname, type = "prob",
                           threshold = threshold, suppressAllGaps = TRUE,
                           replaceIndel = "", conspath = conspath)
    ##
    self$mapIter$`0`[[hptype]] = MapList_(
      ## mapdata
      readpath  = self$relPath(fqout),
      refpath   = refpath(self$mapInit),
      bampath   = bampath(self$mapInit),
      conspath  = self$relPath(conspath),
      pileup    = self$mapInit$pileup,
      stats     = list(coverage = setNames(nfq, "50%")),
      ## required metadata
      maplabel  = maplabel,
      refname   = refname(self$mapInit),
      mapper    = meta(self$mapInit, "mapper"),
      opts      = meta(self$mapInit, "opts")
    )
  }

  return(invisible(self))
})

## Method: mapIter ####
#' @export
mapIter.DR2S <- function(x,
                         opts = list(),
                         iterations = 1,
                         columnOccupancy = 0.4,
                         force = FALSE,
                         plot = TRUE,
                         ...) {
  x$runMapIter(opts = opts,
               iterations = iterations,
               columnOccupancy = columnOccupancy,
               force = force,
               plot = plot,
               ...)
  return(invisible(x))
}

DR2S_$set(
  "public", "runMapIter",
  function(opts = list(),
           iterations = 1,
           columnOccupancy = 0.4,
           force = FALSE,
           plot = TRUE,
           ...) {

    # # debug
    # self <- dr2s
    # opts = list()
    # iterations = 2
    # columnOccupancy = 0.4
    # force = FALSE
    # plot = TRUE
    # ##

    flog.info("# MapIter", name = "info")

    ## Initiate indenter
    indent <- indentation(1)
    flog.info("%sIterative mapping of partitioned longreads", indent(), name = "info")

    ## Collect star t time for mapIter runstats
    start.time <- Sys.time()

    ## Check if reporting is already finished and exit safely
    if (.checkReportStatus(self)) return(invisible(self))

    ## Overide default arguments
    args <- self$getOpts("mapIter")
    if (!is.null(args)) {
      env  <- environment()
      list2env(args, envir = env)
    }
    if (!exists("callInsertionThreshold")) {
      callInsertionThreshold <- 1/5
    }

    threshold  <- self$getThreshold()
    iterations <- self$getIterations()
    baseoutdir <- self$absPath("mapIter")
    includeInsertions <- ifelse(self$hasShortreads(), FALSE, TRUE)
    callInsertions    <- ifelse(self$hasShortreads(), FALSE, TRUE)
    ## Mapper
    mapfun <- self$getLrdMapFun()

    # iterations <- 2
    # iteration <- 1
    # iteration <- 2
    for (iteration in seq_len(iterations)) {
      flog.info("%sIteration %s of %s", indent(), iteration, iterations, name = "info")
      iterationC <- toString(iteration)
      maplabel   <- "mapIter" %<<% iterationC
      prevIteration <- self$mapIter[[toString(iteration - 1)]]
      # hptype = "A"
      # hptype = "B"
      indent2 <- incr(indent)
      foreach(hptype = self$getHapTypes()) %do% {
        flog.info("%sFor haplotype <%s>:", indent2(), hptype, name = "info")
        readtype <- self$getLrdType()
        readfile <- self$absPath(readpath(prevIteration[[hptype]]))
        reffile  <- self$absPath(conspath(prevIteration[[hptype]]))
        refname  <- consname(prevIteration[[hptype]])
        outdir   <- file.path(baseoutdir, hptype)
        pileup <- mapReads(
          mapfun = mapfun, maplabel = maplabel, reffile = reffile, refname = refname,
          readfile = readfile, readtype = readtype, opts = opts, outdir = outdir,
          includeDeletions = TRUE, includeInsertions = includeInsertions,
          callInsertions = callInsertions, callInsertionThreshold = callInsertionThreshold,
          clip = FALSE, distributeGaps = TRUE, removeError = TRUE, topx = 0, force = force,
          clean = clean, indent = incr(indent2), ...)

        ## Construct consensus sequence
        consname <- maplabel %<<% ".consensus." %<<% hptype
        conspath <- file.path(outdir, consname %<<% ".fa")
        names(conspath) <- self$relPath(conspath)
        flog.info("%sConstruct consensus <%s>", indent2(), names(conspath), name = "info")
        conseq <- .writeConseq(x = pileup, name = consname, type = "prob",
                               threshold = threshold, suppressAllGaps = FALSE,
                               suppressInsGaps = TRUE, columnOccupancy = columnOccupancy,
                               replaceIndel = "", conspath = conspath)

        ## Initialize mapIter MapList
        self$mapIter[[iterationC]][[hptype]] = MapList_(
          ## mapdata
          readpath  = self$relPath(readfile),
          refpath   = self$relPath(reffile),
          bampath   = self$relPath(bampath(pileup)),
          conspath  = self$relPath(conspath),
          pileup    = pileup,
          stats     = list(coverage = .coverage(pileup)),
          ## required metadata
          maplabel  = maplabel,
          refname   = refname,
          mapper    = self$getLrdMapper(),
          opts      = opts
        )
      }
    }

    if (plot) {
      flog.info("%sPlot MapIter summary", indent(), name = "info")
      ## Coverage and base frequency
      plotlist <- foreach(iteration = seq_len(self$getIterations())) %do% {
        suppressWarnings(self$plotMapIterSummary(thin = 0.1, width = 4,
                         iteration = iteration, drop.indels = TRUE))
      }
      p <- cowplot::plot_grid(plotlist = plotlist, nrow = self$getIterations())
      cowplot::save_plot(p, filename = self$absPath("plot.MapIter.pdf"),
                         base_width = 12*length(self$getHapTypes()),
                         base_height = 3*self$getIterations(),
                         title = paste(self$getLocus(), self$getSampleId(), sep = "." ))
      cowplot::save_plot(p, filename = self$absPath(".plots/plot.MapIter.svg"),
                         base_width = 12*length(self$getHapTypes()),
                         base_height = 3*self$getIterations())

    }

    createIgvConfigs(x = self, map = "mapIter", open = "FALSE")

    ## set mapIter runstats
    .setRunstats(self, "mapIter",
                 list(Runtime = format(Sys.time() - start.time)))

    return(invisible(self))
  })

## Method: partitionShortreads ####
#' @export
partitionShortreads.DR2S <- function(x,
                                     opts = list(),
                                     force = FALSE,
                                     ...) {
  x$runPartitionShortreads(opts = opts,
                           force = force,
                           ...)
  invisible(x)
}
# TODO: look at arguments and make same
DR2S_$set("public", "runPartitionShortreads", function(opts = list(),
                                                       force = FALSE,
                                                       ...) {
  ## debug
  # opts = list()
  # force = FALSE

  flog.info("# PartitionShortreads ...", name = "info")

  ## Collect start time for partitionShortreads runstats
  start.time <- Sys.time()

  ## Initiate indenter
  indent <- indentation(1)

  ## exit savely if shortreads not provided
  if (!self$hasShortreads()) {
    flog.warn("%sCannot partition shortreads. No shortreads provided", indent(), name = "info")
    return(invisible(self))
  }

  ## exit savely if initial SR mapping not performed
  if (is.null(meta(self$mapInit, "SR2"))) {
    flog.warn("%sCannot partition shortreads. Run 'mapInit()' first", indent(), name = "info")
    return(invisible(self))
  }

  ## exit safely if reporting is already finished
  if (.checkReportStatus(self)) return(invisible(self))

  flog.info("%sPartition shortreads based on initial mapping and " %<<%
            "longread clustering", indent(), name = "info")

  ## Overide default arguments
  args <- self$getOpts("partitionShortreads")
  if (!is.null(args)) {
    env  <- environment()
    list2env(args, envir = env)
  }

  bamfile <- self$absPath(bampath(meta(self$mapInit, "SR2")))
  hptypes <- self$getHapTypes()
  prtMat <- self$lrpartition$mat
  seqs <- lapply(self$lrpartition$hpl, function(x) {
    .getSeqsFromMat(as.matrix(prtMat[x, ]))
  })
  mats <- lapply(seqs, function(x, names) {
    magrittr::set_colnames(createPWM(x), names)
  }, names = colnames(prtMat))

  # Run partitioning
  srpartition <- getSRPartitionScores(bamfile, mats)

  ## Assign read to haplotype with highest probability,
  ## i.e. product over probabilities of each haplotype and choose max
  flog.info("%sGet highest-scoring haplotype for each read", indent(), name = "info")
  srpartition$haplotypes <- scoreHighestSR(srpartition$srpartition, diffThreshold = 0.001)

  # Write fastqs
  # hp <- "A"
  indent2 <- incr(indent)
  foreach(hp = hptypes) %do% {
    srfilenames <- c()
    flog.info("%sWrite shortread fastq for haplotype <%s>", indent2(), hp, name = "info")
    fqs <- self$getShortreads()
    dontUse <- dplyr::filter(srpartition$haplotypes, haplotype != hp)$read
    # fq <- fqs[1]
    foreach(fq = fqs) %do% {
      fqPart <- self$absPath(
        file.path(
          dirname(readpath(self$mapIter$`0`[[hp]])),
          dot(c(strsplit1(basename(fq), ".", fixed = TRUE)[1], hp, "fastq.gz"))
        ))
      .writePartFq(fq, fqPart, dontUse = dontUse, indent = incr(indent2))
      srfilenames <- c(srfilenames, fqPart)
      self$srpartition[[hp]]$srpartition <- srpartition
      NULL
    }
    self$srpartition[[hp]]$SR <- self$relPath(srfilenames)
    NULL
  }

  ## set partitionShortreads runstats
  ## set mapIter runstats
  .setRunstats(self, "partitionShortreads",
               list(Runtime = format(Sys.time() - start.time)))

  return(invisible(self))
})

## Method: mapFinal ####
#' @export
mapFinal.DR2S <- function(x,
                          opts = list(),
                          includeDeletions = TRUE,
                          includeInsertions = TRUE,
                          force = FALSE,
                          plot = TRUE,
                          createIgv = TRUE,
                          clip = FALSE,
                          ...) {
  x$runMapFinal(opts = opts,
                includeDeletions = includeDeletions,
                includeInsertions = includeInsertions,
                force = force,
                plot = plot,
                createIgv = createIgv,
                clip = clip,
                ...)
  invisible(x)
}

DR2S_$set("public", "runMapFinal", function(opts = list(),
                                            includeDeletions = TRUE,
                                            includeInsertions = TRUE,
                                            force = FALSE,
                                            plot = TRUE,
                                            createIgv = TRUE,
                                            clip = FALSE,
                                            ...) {

  ## debug
  # opts = list()
  # includeDeletions = TRUE
  # includeInsertions = TRUE
  # force = FALSE
  # plot = TRUE
  # createIgv = TRUE
  # clip = FALSE
  # library(futile.logger)
  # library(foreach)
  # library(rlang)
  # self <- dr2s

  flog.info("# mapFinal", name = "info")

  ## Collect start time for mapFinal runstats
  start.time <- Sys.time()

  ## Initiate indenter
  indent <- indentation(1)
  flog.info("%sMap longreads and shortreads against mapIter consensus sequences", indent(), name = "info")

  ## Check if reporting is already finished and exit safely
  if (!force && .checkReportStatus(self)) return(invisible(self))

  ## Overide default arguments
  args <- self$getOpts("mapFinal")
  if (!is.null(args)) {
    env  <- environment()
    list2env(args, envir = env)
  }
  if (!exists("callInsertionThreshold")) {
    callInsertionThreshold <- 1/5
  }
  if (!exists("columnOccupancy")) {
    columnOccupancy <- 2/5
  }

  igv <- list()
  maplabel <- "mapFinal"
  outdir   <- .dirCreateIfNotExists(self$absPath(maplabel))
  lastIter <- self$mapIter[[max(names(self$mapIter))]]
  hptypes  <- self$getHapTypes()
  reffiles <- stats::setNames(lapply(hptypes, function(x)
    self$absPath(conspath(lastIter[[x]]))), hptypes)
  readfilesLR <- stats::setNames(lapply(hptypes, function(x)
    self$absPath(readpath(lastIter[[x]]))), hptypes)
  if (self$hasShortreads()) {
    readfilesSR <- stats::setNames(lapply(hptypes, function(x)
      self$absPath(self$srpartition[[x]]$SR)), hptypes)
  }
  ## hp = "A"
  ## hp = "B"
  for (hp in hptypes) {
    flog.info("%sFor haplotype <%s>", indent(), hp, name = "info" )
    reffile <- reffiles[[hp]]
    ##
    ## [1] Map longreads
    ##
    refname  <- "LR" %<<% hp
    mapfun   <- self$getLrdMapFun()
    readfile <- readfilesLR[[hp]]
    readtype <- self$getLrdType()
    pileup <- mapReads(
      mapfun = mapfun, maplabel = maplabel, reffile = reffile, refname = refname,
      readfile = readfile, readtype = readtype, opts = opts, outdir = outdir,
      includeDeletions = includeDeletions, includeInsertions = includeInsertions,
      callInsertions = FALSE, callInsertionThreshold = callInsertionThreshold,
      clip = FALSE, distributeGaps = TRUE, removeError = TRUE, topx = 0,
      force = force, clean = TRUE,  max_depth = 1e4, min_mapq = 0,
      indent = incr(indent), ...)
    ## Create igv
    if (createIgv) {
      igv <- createIgvJsFiles(
        refpath(pileup), bampath(pileup), self$getOutdir(), sampleSize = 100,
        fragmentReads = TRUE)
    }
    ## Construct consensus sequence
    consname <- maplabel %<<% ".consensus." %<<% refname
    conspath <- file.path(outdir, consname %<<% ".fa")
    names(conspath) <- self$relPath(conspath)
    flog.info("%sConstruct consensus <%s>", indent(), names(conspath), name = "info")
    conseq <- .writeConseq(x = pileup, name = consname, type = "ambig",
                           threshold = 1/4, suppressAllGaps = TRUE,
                           replaceIndel = "", conspath = conspath)
    ## Initialize mapFinal LR MapList
    self$mapFinal$LR[[hp]] = MapList_(
      ## mapdata
      readpath  = self$relPath(readfile),
      refpath   = self$relPath(reffile),
      bampath   = self$relPath(bampath(pileup)),
      conspath  = self$relPath(conspath),
      pileup    = pileup,
      stats     = list(coverage = .coverage(pileup)),
      ## required metadata
      maplabel  = maplabel,
      refname   = refname,
      mapper    = self$getLrdMapper(),
      opts      = opts,
      ## additional metadata
      igv       = igv
    )

    if (self$hasShortreads()) {
      ##
      ## [2] Map shortreads
      ##
      refname  <- "SR" %<<% hp
      mapfun   <- self$getSrdMapFun()
      readfile <- readfilesSR[[hp]]
      readtype <- self$getSrdType()
      pileup <- mapReads(
        mapfun = mapfun, maplabel = maplabel, reffile = reffile, refname = refname,
        readfile = readfile, readtype = readtype, opts = opts, outdir = outdir,
        includeDeletions = includeDeletions, includeInsertions = includeInsertions,
        callInsertions = TRUE, callInsertionThreshold = callInsertionThreshold,
        clip = FALSE, distributeGaps = TRUE, removeError = TRUE, topx = 0,
        force = force, clean = TRUE, max_depth = 1e5, min_mapq = 50,
        min_base_quality = 13, indent = incr(indent), ...)
      ## Create igv
      if (createIgv) {
        igv <- createIgvJsFiles(
          refpath(pileup), bampath(pileup), self$getOutdir(), sampleSize = 100)
      }
      ## Construct consensus sequence
      consname <- maplabel %<<% ".consensus." %<<% refname
      conspath <- file.path(outdir, consname %<<% ".fa")
      names(conspath) <- self$relPath(conspath)
      flog.info("%sConstruct consensus <%s>", indent(), names(conspath), name = "info")
      conseq <- .writeConseq(x = pileup, name = consname, type = "ambig",
                             threshold = 1/4, suppressAllGaps = TRUE,
                             replaceIndel = "", conspath = conspath)
      ## Initialize mapFinal SR MapList
      self$mapFinal$SR[[hp]] = MapList_(
        ## mapdata
        readpath  = self$relPath(readfile),
        refpath   = self$relPath(reffile),
        bampath   = self$relPath(bampath(pileup)),
        conspath  = self$relPath(conspath),
        pileup    = pileup,
        stats     = list(coverage = .coverage(pileup)),
        ## required metadata
        maplabel  = maplabel,
        refname   = refname,
        mapper    = self$getSrdMapper(),
        opts      = opts,
        ## additional metadata
        igv       = igv
      )
    }
  }

  if (plot) {
    flog.info("%sPlot MapFinal summary", indent(), name = "info")
    ## Coverage and base frequency
    readtypes <- if (self$hasShortreads()) c("LR", "SR") else "LR"
    plotRows  <- if (self$hasShortreads()) 2 else 1
    ## readtype = "LR"
    plotlist <- foreach(readtype = readtypes) %do% {
      suppressWarnings(self$plotMapFinalSummary(iteration = "final",
                       readtype = readtype, thin = 0.25, width = 20))
    }
    p <- cowplot::plot_grid(plotlist = plotlist, nrow = plotRows, labels = readtypes)
    cowplot::save_plot(p, filename = self$absPath("plot.MapFinal.pdf"),
                       base_width = 12*length(hptypes),
                       title = paste(self$getLocus(), self$getSampleId(), sep = "." ),
                       base_height = 3*length(readtypes))
    cowplot::save_plot(p, filename = self$absPath(".plots/plot.MapFinal.svg"),
                       base_width = 12*length(hptypes),
                       base_height = 3*length(readtypes))
  }

  ## set mapFinal runstats
  .setRunstats(self, "mapFinal",
               list(Runtime = format(Sys.time() - start.time)))

  return(invisible(self))
})

## Method: runPipeline ####
DR2S_$set("public", "runPipeline", function() {
  steps_ <- self$getPipeline()

  ## Collect start time for DR2Spipeline runstats
  start.time <- Sys.time()

  while (length(steps_) > 0) {
    step <- steps_[1]
    self$run_(step)
    steps_ <- steps_[-1]
  }

  ## set DR2Spipeline runstats
  .setRunstats(self, "DR2Spipeline",
               list(Runtime = format(Sys.time() - start.time)))

  return(invisible(self))
})

