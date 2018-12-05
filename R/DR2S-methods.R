# Method: mapInit ####
#' @export
mapInit.DR2S <- function(x, opts = list(), ...) {
  x$runMapInit(opts = opts, ...)
  invisible(x)
}

DR2S_$set("public", "runMapInit", function(opts = list(), ...) {
  # # debug
  # opts = list()
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

  ## Export mapInit config to function environment
  args <- self$getOpts("mapInit")
  list2env(args, envir = environment())
  assert_that(
    exists("includeDeletions") && is.logical(includeDeletions),
    exists("includeInsertions") && is.logical(includeInsertions),
    exists("callInsertionThreshold") && is.numeric(callInsertionThreshold),
    exists("minMapq") && is.numeric(minMapq),
    exists("topx"),
    exists("pickiness") && is.numeric(pickiness),
    exists("increasePickiness") && is.numeric(increasePickiness),
    exists("lowerLimit") && is.numeric(lowerLimit),
    exists("updateBackgroundModel") && is.logical(updateBackgroundModel),
    exists("createIgv") && is.logical(createIgv),
    exists("plot") && is.logical(plot)
  )

  ## Get options and prepare mapping
  outdir <- .dirCreateIfNotExists(self$absPath("mapInit"))
  .dirCreateIfNotExists(file.path(self$absPath(".plots")))
  clean <- TRUE
  pickedTopX <- FALSE
  igv <- list()
  SR <- list()

  if (self$hasShortreads()) {
    SR <- mapInitSR(
      self = self, opts = opts, includeDeletions = includeDeletions,
      includeInsertions = includeInsertions, callInsertions = TRUE,
      callInsertionThreshold = callInsertionThreshold, clip = FALSE,
      distributeGaps = FALSE, removeError = TRUE, topx = 0,
      outdir = outdir, clean = clean, minMapq = 50, indent = indent, ...)

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
      distributeGaps = TRUE, removeError = TRUE, topx = topx,
      updateBackgroundModel = updateBackgroundModel, clean = clean,
      minMapq = minMapq, pickiness = pickiness, lowerLimit = lowerLimit,
      indent = indent, ...)

    ## If topx = "auto" or set to non-zero generate a new subsampled
    ## read file
    if (!is.null(picked1 <- reads(pileup))) {
      fqFile <- dot(c(self$getSampleId(), readtype, "n" %<<% length(picked1), "fastq", "gz"))
      topxFqPath <- .fileDeleteIfExists(file.path(outdir, strip(fqFile, "_")))
      names(topxFqPath) <- self$relPath(topxFqPath)
      fqReads <- .extractFastq(bampath(pileup), picked1)
      ShortRead::writeFastq(fqReads, topxFqPath, compress = TRUE)
      ## picking plot 1
      plotpath <- file.path(outdir, "plot.readpicking1.png")
      cowplot::save_plot(plotpath, attr(picked1, "plot"), dpi = 150)
      ## set flag
      pickedTopX <- TRUE
    }

    ## Construct initial longread consensus sequence
    refname <- refname %<<% ".consensus"
    self$setRefPath(file.path(outdir, strip(maplabel %<<% "." %<<% refname %<<% ".fa", "_")))
    reffile <- self$getRefPath()
    flog.info("%sConstruct consensus <%s>", indent2(), names(reffile), name = "info")
    conseq <- .writeConseq(x = pileup, name = refname, type = "prob",
                           threshold = NULL, suppressAllGaps = TRUE,
                           replaceIndel = "N", conspath = reffile)
  }

  maplabel  <- "mapInit2"
  mapfun    <- self$getLrdMapFun()
  readtype  <- self$getLrdType()
  indelRate <- if (self$hasShortreads()) indelRate(SR$SR2$pileup) else indelRate(pileup)
  if (pickedTopX) {
    readfile  <- topxFqPath
    pickiness <- pickiness/increasePickiness
  } else {
    readfile <- self$getLongreads()
  }
  pileup <- mapReads(
    mapfun = mapfun, maplabel = maplabel, reffile = reffile, refname = refname,
    readfile = readfile, readtype = readtype, opts = opts, outdir = outdir,
    includeDeletions = TRUE, includeInsertions = FALSE, callInsertions = FALSE,
    clip = FALSE, distributeGaps = TRUE, removeError = TRUE, topx = topx,
    updateBackgroundModel = updateBackgroundModel, clean = clean,
    minMapq = minMapq, pickiness = pickiness, lowerLimit = lowerLimit,
    indent = indent, indelRate = indelRate, ...)

  if (!is.null(picked2 <- reads(pileup))) {
    ## picking plot 2
    plotpath <- file.path(outdir, "plot.readpicking2.png")
    cowplot::save_plot(plotpath, attr(picked2, "plot"), dpi = 150)
  }

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
    p <- self$plotMapInitSummary(thin = 0.25, width = 2)
    plotRows <- ifelse(self$hasShortreads(), 2, 1)
    cowplot::save_plot(self$absPath("plot.mapInit.png"),
                       plot = p, ncol = 1, nrow = plotRows, dpi = 150,
                       base_aspect_ratio = as.numeric(dot(c(5, plotRows))),
                       title = dot(c(self$getLocus(), self$getSampleId())))
    cowplot::save_plot(self$absPath(".plots/plot.mapInit.svg"),
                       plot = p, ncol = 1, nrow = plotRows,
                       base_aspect_ratio = as.numeric(dot(c(5, plotRows))))
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
partitionLongreads.DR2S <- function(x) {
  ## Collect start time for partitionLongreads runstats
  start.time <- Sys.time()
  x$runPartitionLongreads()
  x$runSplitLongreadsByHaplotype()
  x$runExtractPartitionedLongreads()
  ## set partitionLongreads runstats
  .setRunstats(x, "partitionLongreads",
               list(Runtime = format(Sys.time() - start.time)))
  return(invisible(x))
}

DR2S_$set("public", "runPartitionLongreads", function() {
    ## debug
    # library(futile.logger)
    # library(assertthat)
    # self <- dr2s

    flog.info("# PartitionLongreads", name = "info")

    ## Initiate indenter
    indent <- indentation(1)

    ## Export partitionLongreads config to function environment
    args <- self$getOpts("partitionLongreads")
    list2env(args, envir = environment())
    assert_that(
      self$hasMapInit(),
      exists("threshold") && is.double(threshold),
      exists("distAlleles") && is.count(distAlleles),
      exists("skipGapFreq") && is.numeric(skipGapFreq),
      exists("noGapPartitioning") && is.logical(noGapPartitioning),
      exists("restrictToCorrelatedPositions") && is.logical(restrictToCorrelatedPositions),
      exists("measureOfAssociation") && is.character(measureOfAssociation),
      exists("expectedAbsDeviation") && is.numeric(expectedAbsDeviation),
      exists("selectAllelesBy") && is.character(selectAllelesBy),
      exists("minClusterSize") && is.numeric(minClusterSize),
      exists("plot") && is.logical(plot)
    )

    ## Get the reference sequence
    if (self$hasShortreads()) {
      useSR <- TRUE
      flog.info("%sConstruct SNP matrix from shortreads", indent(), name = "info")
    } else {
      useSR <- FALSE
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
    base_height <- max(6, floor(sqrt(NCOL(mat))))
    if (restrictToCorrelatedPositions) {
      spos <- .selectAssociatedPolymorphicPositions(
        mat, method.assoc = measureOfAssociation, method.clust = "mclust",
        expectedAbsDeviation = expectedAbsDeviation, indent = indent)
      mat0 <- mat[, spos[order(as.numeric(spos))]]
      if (plot) {
        ## correlogram
        plotpath <- file.path(self$getOutdir(), "plot.correlogram.png")
        cowplot::save_plot(plotpath, attr(spos, "snp.correlogram"),
                           base_height = base_height, dpi = 150)
        ## associiation plot
        plotpath <- file.path(self$getOutdir(), "plot.association.png")
        cowplot::save_plot(plotpath, attr(spos, "snp.association"),
                           base_height = base_height/2.4, dpi = 150,
                           base_aspect_ratio = 3)
      }
    } else {
      mat0 <- mat
    }

    flog.info("%sPartition %s longreads over %s SNPs", indent(), NROW(mat0), NCOL(mat0), name = "info")
    prt <- partitionReads(x = mat0,
                          skipGapFreq = skipGapFreq,
                          deepSplit = 1,
                          threshold = threshold,
                          distAlleles = distAlleles,
                          sortBy = selectAllelesBy,
                          minClusterSize = minClusterSize,
                          indent = incr(indent))

    if (plot) {
      ppos <- SNP(prt)
      pwm <- lapply(PWM(prt), function(pwm) {
        pwm[pwm < 0.1] <- 0
        pwm
      })
      p <- suppressMessages(self$plotSeqLogo(ppos, pwm))
      cowplot::save_plot(filename = self$absPath("plot.sequence.png"),
                         plot = p, dpi = 150, units = "cm", limitsize = FALSE,
                         base_width = 0.4*length(ppos) + 1.4,
                         base_height = 2.5*length(pwm),
                         title = dot(c(self$getLocus(), self$getSampleId())))
      cowplot::save_plot(filename = self$absPath(".plots/plot.sequence.svg"),
                         plot = p, units = "cm", limitsize = FALSE,
                         base_width  = 0.4*length(ppos) + 1.4,
                         base_height = 2.5*length(pwm))
    }

    ## set runstats
    .setRunstats(self, "partitionLongreads",
                 list(nLongreads = NROW(mat),
                      foundPolymorphicPositions = NCOL(mat),
                      usedPolymorphicPositions = K(prt),
                      foundClades = OC(prt),
                      usedClades = as.list(table(PRT(prt)))))

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

DR2S_$set("public", "runSplitLongreadsByHaplotype", function() {
  ## If reporting is already done exit safely
  if (.checkReportStatus(self)) return(invisible(self))

  assert_that(self$hasLrdPartition())

  ## Initiate indenter
  indent <- indentation(1)
  flog.info("%sSplit partitioned longreads by score", indent(), name = "info")

  ## Export partitionLongreads config to function environment
  args <- self$getOpts("partitionLongreads")
  list2env(args, envir = environment())
  assert_that(
    exists("pickiness") && is.numeric(pickiness),
    exists("lowerLimit") && is.numeric(lowerLimit),
    exists("plot") && is.logical(plot)
  )

  prt <- partition(self$getPartition())
  haplotypes <- levels(prt$haplotype)
  tag <- self$getMapTag("init", "LR")

  # Set all limits to NULL
  self$setLimits(rep(NA, length(haplotypes)))
  prts <- lapply(haplotypes, function(hp) dplyr::filter(prt, .data$haplotype == hp))
  names(prts) <- haplotypes
  scores <- lapply(prts, function(x) x$mcoef)
  lmts <- .optimalPartitionLimits(scores, pickiness, lowerLimit)
  self$setLimits(setNames(as.list(lmts$limits$score), haplotypes))
  self$lrpartition$lmt <- lmts$plt

  # Get only reads within the limit
  reads <- stats::setNames(lapply(names(self$getLimits()), function(x) {
    dplyr::filter(prt, .data$haplotype == x, .data$mcoef >= self$getLimits()[x])
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
               list(reads = setNames(as.list(as.integer(lmts$limits$nreads)), haplotypes)))

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

    cowplot::save_plot(self$absPath("plot.partition.png"), plot = p, dpi = 150,
                       base_height = 12, base_width = 10,
                       title = dot(c(self$getLocus(), self$getSampleId())))
    cowplot::save_plot(self$absPath(".plots/plot.partition.svg"), plot = p,
                       base_aspect_ratio = 1.2)

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

  ## Export PartitionLongreadsconfig to function environment
  args <- self$getOpts("partitionLongreads")
  list2env(args, envir = environment())

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
    fqfile <- dot(c(
      "hap" %<<% hptype, self$getLrdType(), self$getLrdMapper(),
      "n" %<<% length(fq), "fastq", "gz"))
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
                           threshold = NULL, suppressAllGaps = TRUE,
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
mapIter.DR2S <- function(x, opts = list(), ...) {
  x$runMapIter(opts = opts, ...)
  return(invisible(x))
}

DR2S_$set("public", "runMapIter", function(opts = list(), ...) {
  # debug
  # self <- dr2s
  # opts = list()
  # plot = TRUE
  ## If reporting is already done exit safely
  if (.checkReportStatus(self)) return(invisible(self))

  ## Collect start time for mapIter runstats
  start.time <- Sys.time()

  flog.info("# MapIter", name = "info")

  ## Initiate indenter
  indent <- indentation(1)
  flog.info("%sIterative mapping of partitioned longreads", indent(), name = "info")

  ## Export mapIter config to function environment
  args <- self$getOpts("mapIter")
  list2env(args, envir = environment())
  assert_that(
    exists("iterations") && is.count(iterations),
    exists("columnOccupancy") && is.double(columnOccupancy),
    exists("callInsertionThreshold") && is.double(callInsertionThreshold),
    exists("plot") && is.logical(plot)
  )

  baseoutdir <- self$absPath("mapIter")
  clean <- TRUE
  includeInsertions <- ifelse(self$hasShortreads(), FALSE, TRUE)
  callInsertions <- ifelse(self$hasShortreads(), FALSE, TRUE)
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
        clip = FALSE, distributeGaps = TRUE, removeError = TRUE, topx = 0,
        clean = clean, indent = incr(indent2), ...)

      ## Construct consensus sequence
      consname <- maplabel %<<% ".consensus." %<<% hptype
      conspath <- file.path(outdir, consname %<<% ".fa")
      names(conspath) <- self$relPath(conspath)
      flog.info("%sConstruct consensus <%s>", indent2(), names(conspath), name = "info")
      conseq <- .writeConseq(x = pileup, name = consname, type = "prob",
                             threshold = NULL, suppressAllGaps = FALSE,
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
    cowplot::save_plot(p, dpi = 150, filename = self$absPath("plot.mapIter.png"),
                       base_width = 12*length(self$getHapTypes()),
                       base_height = 3*self$getIterations(),
                       title = dot(c(self$getLocus(), self$getSampleId())))
    cowplot::save_plot(p, filename = self$absPath(".plots/plot.mapIter.svg"),
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
partitionShortreads.DR2S <- function(x, opts = list(), ...) {
  x$runPartitionShortreads(opts = opts, ...)
  invisible(x)
}

# TODO: look at arguments and make same
DR2S_$set("public", "runPartitionShortreads", function(opts = list(), ...) {
  ## debug
  # opts = list()
  ## If reporting is already done exit safely
  if (.checkReportStatus(self)) return(invisible(self))

  ## Collect start time for partitionShortreads runstats
  start.time <- Sys.time()

  flog.info("# PartitionShortreads ...", name = "info")

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

  flog.info("%sPartition shortreads based on initial mapping and " %<<%
              "longread clustering", indent(), name = "info")

  ## Export partitionShortreads config to function environment
  args <- self$getOpts("partitionShortreads")
  list2env(args, envir = environment())

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
mapFinal.DR2S <- function(x, opts = list(), ...) {
  x$runMapFinal(opts = opts, ...)
  invisible(x)
}

DR2S_$set("public", "runMapFinal", function(opts = list(), ...) {
  ## debug
  # opts = list()
  # library(futile.logger)
  # library(foreach)
  # self <- dr2s
  ## If reporting is already done exit safely
  if (.checkReportStatus(self)) return(invisible(self))

  ## Collect start time for mapFinal runstats
  start.time <- Sys.time()

  flog.info("# mapFinal", name = "info")

  ## Initiate indenter
  indent <- indentation(1)
  flog.info("%sMap longreads and shortreads against mapIter consensus sequences", indent(), name = "info")

  ## Export mapFinal config to function environment
  args <- self$getOpts("mapFinal")
  list2env(args, envir = environment())
  assert_that(
    exists("includeDeletions") && is.logical(includeDeletions),
    exists("includeInsertions") && is.logical(includeInsertions),
    exists("callInsertionThreshold") && is.numeric(callInsertionThreshold),
    exists("trimPolymorphicEnds") && is.logical(trimPolymorphicEnds),
    exists("createIgv") && is.logical(createIgv),
    exists("plot") && is.logical(plot)
  )

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
      clean = TRUE, max_depth = 1e4, min_mapq = 0, indent = incr(indent), ...)
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
        clip = trimPolymorphicEnds, distributeGaps = TRUE, removeError = TRUE, topx = 0,
        clean = TRUE, max_depth = 1e5, min_mapq = 50, min_base_quality = 13,
        indent = incr(indent), ...)
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
    cowplot::save_plot(p, dpi = 150, filename = self$absPath("plot.mapFinal.png"),
                       base_width = 12*length(hptypes),
                       base_height = 3*length(readtypes),
                       title = dot(c(self$getLocus(), self$getSampleId())))
    cowplot::save_plot(p, filename = self$absPath(".plots/plot.mapFinal.svg"),
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

