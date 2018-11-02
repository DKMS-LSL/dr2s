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
  indent <- indentation(1)

  ## Overide default arguments
  args <- self$getOpts("mapInit")
  if (!is.null(args)) {
    env  <- environment()
    list2env(args, envir = env)
  }
  if (is.null(threshold)) {
    threshold <- self$getThreshold()
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
      includeInsertions = includeInsertions, callInsertions = TRUE, clip = FALSE,
      distributeGaps = FALSE, removeError = TRUE, topx = 0, outdir = outdir,
      force = force, clean = clean, minMapq = 50, indent = indent, ...)

    ### TODO wrap this command up
    if (createIgv)
      igv[["SR"]] <- createIgvJsFiles(
        reference = refpath(SR$mapInitSR2$pileup),
        bamfile = path(SR$mapInitSR2$pileup),
        outdir = outdir,
        sampleSize = 100)

    readfile <- self$getLongreads()
    reffile  <- self$absPath(SR$mapInitSR1$seqpath)
    allele   <- SR$mapInitSR1$ref
  }
  else {
    mapLabel <- "mapInit"
    reffile  <- self$getRefPath()
    allele   <- self$getReference()
    readfile <- self$getLongreads()
    readtype <- self$getLrdType()
    mapFun   <- self$getLrMapFun()
    maptag   <- paste(mapLabel, paste0(litArrows(c(allele, readtype,
                                                   self$getLrdMapper(),
                                                   optstring(opts))),
                                       collapse = " "))
    pileup <- mapReads(
      mapFun = mapFun, maptag = maptag, reffile = reffile,
      allele = allele, readfile = readfile, readtype = readtype, opts = opts,
      refname = "",  includeDeletions = TRUE, includeInsertions = TRUE,
      callInsertions = TRUE, clip = FALSE, distributeGaps = TRUE,
      removeError = TRUE, topx = topx, outdir = outdir, force = force,
      clean = clean, indent, ...)

    if (!is.null(reads(pileup))) {
      file <- paste(
        self$getSampleId(), self$getLrdType(), paste0("n", topx),
        "fastq", "gz", sep = ".")
      readfile <- .fileDeleteIfExists(file.path(outdir, file))
      names(readfile) <- self$relPath(readfile)
      fq <- .extractFastq(path(pileup), reads(pileup))
      ShortRead::writeFastq(fq, readfile, compress = TRUE)
    }

    ## Construct initial longread consensus sequence
    allele <- "Init.LRconsensus." %<<% sub(".bam", "", basename(path(pileup)))
    self$setRefPath(file.path(basename(outdir), allele %<<% ".fa"))
    reffile <- self$getRefPath()
    flog.info("%sConstruct consensus <%s>", indent(), names(reffile), name = "info")
    conseq <- .getWriteConseq(pileup = pileup, name = "mapInitLR", type = "prob",
                              threshold = threshold, suppressAllGaps = TRUE,
                              conseqPath = reffile)
  }

  mapLabel <- "mapInit"
  readtype <- self$getLrdType()
  mapFun   <- self$getLrMapFun()
  maptag   <- paste(mapLabel, paste0(litArrows(c(allele, readtype,
                                                 self$getLrdMapper(),
                                                 optstring(opts))),
                                     collapse = " "))
  pileup <- mapReads(
    mapFun = mapFun, maptag = maptag, reffile = reffile,
    allele = allele, readfile = readfile, readtype = readtype, opts = opts,
    refname = "", includeDeletions = TRUE, includeInsertions = FALSE,
    callInsertions = FALSE, clip = FALSE, distributeGaps = TRUE,
    removeError = TRUE, topx = topx, outdir = outdir, force = force,
    clean = clean, indent, ...)

  if (createIgv)
    igv[["LR"]] <- createIgvJsFiles(
      reference = refpath(pileup),
      bamfile = path(pileup),
      outdir = outdir,
      sampleSize = 100,
      fragmentReads = TRUE)

  self$mapInit = structure(
    list(
      ## mandatory fields
      reads   = self$relPath(readfile),
      bamfile = self$relPath(path(pileup)),
      pileup  = pileup,
      tag     = maptag,
      stats   = list(coverage = .coverage(pileup)),
      ## additional fields
      SR1     = SR$mapInitSR1,
      SR2     = SR$mapInitSR2,
      igv     = igv
    ),
    class  = c("mapInit", "list")
  )

  createIgvConfigs(x = self, map = "mapInit", open = "FALSE")

  if (plot) {
    flog.info("%sPlot MapInit summary", indent(), name = "info")
    ## Coverage and frequency of minor alleles
    p <- self$plotmapInitSummary(
      thin = 0.25,
      width = 2
    )
    plotRows <- ifelse(self$hasShortreads(), 2, 1)
    cowplot::save_plot(self$absPath("plot.MapInit.pdf"),
                       plot = p, ncol = 1, nrow = plotRows,
                       base_aspect_ratio = as.numeric(paste(5, plotRows,sep = ".")),
                       title = paste(self$getLocus(), self$getSampleId(), sep = "." ))
    cowplot::save_plot(self$absPath(".plots/plot.MapInit.svg"),
                       plot = p, ncol = 1, nrow = plotRows,
                       base_aspect_ratio = as.numeric(paste(5, plotRows, sep = ".")))
  }

  ## set mapInit runstats
  .setRunstats(self, "mapInit",
               list(Runtime = format(Sys.time() - start.time),
                    SRcoverage = self$mapInit$SR2$stats$coverage[["50%"]],
                    LRcoverage = self$mapInit$stats$coverage[["50%"]]))

  return(invisible(self))
})

#' @export
print.mapInit <- function(x, ...) {
  msg <- sprintf("An object of class '%s'\n", class(x)[1])
  msg <- sprintf(
    "%s [Tag]      %s\n [Reads]    %s\n [Bamfile]  %s\n [Coverage] %s\n",
    msg,
    x$tag,
    paste(basename(x$reads), collapse = "\n            "),
    basename(x$bamfile),
    x$stats$coverage[["50%"]]
  )
  cat(msg)
}

## Method: partitionLongReads ####
#' @export
partitionLongReads.DR2S <- function(x,
                                    threshold         = NULL,
                                    skipGapFreq       = 2/3,
                                    distAlleles       = NULL,
                                    noGapPartitioning = FALSE,
                                    selectAllelesBy   = "count",
                                    minClusterSize    = 15,
                                    plot              = TRUE,
                                    ...) {
  ## Collect start time for partitionLongReads runstats
  start.time <- Sys.time()

  x$runPartitionLongReads(threshold         = threshold,
                          skipGapFreq       = skipGapFreq,
                          noGapPartitioning = noGapPartitioning,
                          selectAllelesBy   = selectAllelesBy,
                          minClusterSize    = minClusterSize,
                          distAlleles       = distAlleles,
                          plot              = plot,
                          ...)
  x$runSplitLongReadsByHaplotype(plot = plot)
  x$runExtractLongReads()
  x$runGetPartitionedConsensus()

  ## set partitionLongReads runstats
  .setRunstats(x, "partitionLongReads",
               list(Runtime = format(Sys.time() - start.time)))

  return(invisible(x))
}

DR2S_$set("public",
          "runPartitionLongReads",
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

    flog.info("# PartitionLongReads", name = "info")

    ## Initiate indenter
    indent <- indentation(1)

    ## Overide default arguments
    args <- self$getOpts("partitionLongReads")
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
      self$hasPileup(),
      is.double(skipGapFreq),
      is.double(threshold),
      is.count(distAlleles),
      is.logical(plot)
    )

    ## Get the reference sequence
    if (!is.null(self$mapInit$SR1)) {
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
      flog.info("%sUse only non-gap positions for clustering", indent(), name = "info")
      ppos <- dplyr::filter(ppos, a1 != "-" & a2 != "-")
    }

    ## Check if already finished because it is a homozygous sample
    if (NROW(ppos) == 0) {
      flog.warn("%sNo polymorphic positions for clustering! Only single allele?", indent(), name = "info")
      flog.info("%sEntering polish and report pipeline", indent(), name = "info")
      ## set runstats
      .setRunstats(self, "partitionLongReads",
                   list(foundPolymorphicPositions = 0L))
      return(invisible(.finishCn1(self)))
    }

    mat <- if (tryCatch(
      !is(self$partition, "PartList"),
      error = function(e)
        TRUE
    ) ||
    !(all(ppos$position %in% colnames(self$partition$mat)) &&
      all(colnames(self$partition$mat) %in% ppos$position))) {
      SNPmatrix(bamfile = self$absPath(self$mapInit$bamfile), polymorphicPositions = ppos)
    } else {
      self$partition$mat
    }

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
    .setRunstats(self, "partitionLongReads",
                 list(nLongreads = NROW(mat),
                      foundPolymorphicPositions = NCOL(mat),
                      usedPolymorphicPositions = K(prt),
                      foundClades = as.list(table(PRT(prt)))))

    ## Set sample haplotypes
    self$setHapTypes(levels(as.factor(PRT(prt))))

    # Check if we have only one cluster and finish the pipeline if so
    if (length(self$getHapTypes()) == 1) {
      flog.warn("%sOnly one allele left!", indent(), name = "info")
      flog.info("%sEntering polish and report pipeline", indent(), name = "info")
      return(invisible(.finishCn1(self)))
    }

    self$partition = structure(list(
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

DR2S_$set("public", "runSplitLongReadsByHaplotype", function(plot = TRUE) {

  ## Initiate indenter
  indent <- indentation(1)

  flog.info("%sSplit partitioned longreads by score", indent(), name = "info")
  ## Check if reporting is already finished and exit safely
  if (.checkReportStatus(self)) return(invisible(self))
  assert_that(self$hasPartition())

  ## Overide default arguments
  args <- self$getOpts("partitionLongReads")
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
  self$partition$lmt <- lmts$plt

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
  .setRunstats(self, "partitionLongReads",
               list(reads = setNames(as.list(as.integer(lmts$limits$r)), haplotypes)))

  # results structure
  self$partition$hpl <- structure(
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

    p <- self$plotPartitionSummary(label = tag,
                                   limits = unlist(self$getLimits()))

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

DR2S_$set("public", "runExtractLongReads", function() {

  ## Initiate indenter
  indent <- indentation(1)

  flog.info("%sExtract partitioned longreads", indent(), name = "info")

  ## Check if reporting is already finished and exit safely
  if (.checkReportStatus(self)) return(invisible(self))
  assert_that(self$hasHapList())

  ## Overide default arguments
  args <- self$getOpts("partitionLongReads")
  if (!is.null(args)) {
    env  <- environment()
    list2env(args, envir = env)
  }

  ## do this for each haptype
  hptypes <- self$getHapTypes()
  for (hptype in hptypes) {
    dir <- .dirCreateIfNotExists(
      normalizePath(file.path(self$getOutdir(), "mapIter", (hptype)), mustWork = FALSE))
    qnames <- self$getHapList(hptype)
    fq  <- .extractFastq(self$absPath(self$mapInit$bamfile), qnames = qnames)
    file <- paste(
      "hap", hptype, self$getLrdType(), self$getLrdMapper(),
      "lim" %<<% as.character(floor(abs(attr(self$getHapList(hptype), "limit")))),
      "n" %<<% length(fq),
      "fastq", "gz", sep = ".")
    out <- .fileDeleteIfExists(file.path(dir, file))
    ShortRead::writeFastq(fq, out, compress = TRUE)
    self$mapIter$`0`[[hptype]] = structure(
      list(
        dir     = self$relPath(dir),
        reads   = self$relPath(out),
        bamfile = NULL,
        pileup  = NULL,
        tag     = NULL,
        conseq  = NULL,
        seqpath = NULL,
        ref     = NULL,
        params  = NULL,
        stats   = NULL
      ),
      class = c("mapIter", "list")
    )
  }

  return(invisible(self))
})


DR2S_$set(
  "public", "runGetPartitionedConsensus",
  function(opts = list()) {
    ## Check if reporting is already finished and exit safely
    if (.checkReportStatus(self)) return(invisible(self))
    ## Initiate indenter
    indent <- indentation(1)
    bamfile <- self$absPath(self$mapInit$bamfile)
    ## Construct consensus from initial mapping with the clustered reads
    flog.info("%sConstruct consensus sequences based on <%s>", indent(), names(bamfile), name = "info")
    mat <- .msaFromBam(Rsamtools::BamFile(bamfile), paddingLetter = ".")
    ## hptype = "A"
    indent2 <- incr(indent)
    foreach(hptype = self$getHapTypes()) %do% {
      flog.info("%sFor haplotype <%s>", indent(), hptype, name = "info")
      readIds <- self$getHapList(hptype)
      cmat <- .extractIdsFromMat(mat, readIds)

      ## Construct consensus sequence
      conseqName <- "consensus.mapIter.0." %<<% hptype
      conseqPath <- self$absPath(
        file.path(self$mapIter$`0`[[hptype]]$dir, conseqName %<<% ".fa"))
      flog.info("%sConstruct consensus <%s>", indent2(), names(conseqPath), name = "info")
      conseq <- conseq(cmat, name = conseqName, type = "prob", suppressInsGaps = FALSE)
      Biostrings::writeXStringSet(
        Biostrings::DNAStringSet(gsub("[-+]", "", conseq)),
        conseqPath)

      self$mapIter$`0`[[hptype]]$tag = "mapIter0"
      self$mapIter$`0`[[hptype]]$ref = "mapIter0"
      self$mapIter$`0`[[hptype]]$conseq  = conseq
      self$mapIter$`0`[[hptype]]$seqpath = self$relPath(conseqPath)
      NULL
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

    threshold <- self$getThreshold()
    hptypes <- self$getHapTypes()
    iterations <- self$getIterations()
    baseoutdir <- self$absPath("mapIter")
    includeInsertions <- ifelse(self$hasShortreads(), FALSE, TRUE)
    callInsertions <- ifelse(self$hasShortreads(), FALSE, TRUE)

    ## Mapper
    mapFun <- self$getLrMapFun()

    # iteration <- 1
    # iteration <- 2
    for (iteration in seq_len(iterations)) {
      flog.info("%sIteration %s of %s", indent(), iteration, iterations, name = "info")
      iterationC <- toString(iteration)
      prevIteration <- self$mapIter[[toString(iteration - 1)]]
      # hptype = "A"
      # hptype = "B"
      indent2 <- incr(indent)
      foreach(hptype = hptypes) %do% {
        flog.info("%sFor haplotype <%s>:", indent2(), hptype, name = "info")
        readtype <- self$getLrdType()
        reffile  <- self$absPath(prevIteration[[hptype]]$seqpath)
        readfile <- self$absPath(prevIteration[[hptype]]$reads)
        allele   <- "mapIter" %<<% iterationC
        refname  <- sprintf("%s.%s", prevIteration[[hptype]]$ref, hptype)
        outdir   <- file.path(baseoutdir, hptype)
        mapLabel <- "mapIter"
        maptag   <- paste(mapLabel,
                          paste0(litArrows(c(iteration, hptype,
                                             readtype, self$getLrdMapper())),
                                 collapse = " "))
        pileup <- mapReads(
          mapFun = mapFun, maptag = maptag, reffile = reffile, allele = allele,
          readfile = readfile, readtype = readtype, opts = opts, refname = refname,
          includeDeletions = TRUE, includeInsertions = includeInsertions,
          callInsertions = callInsertions, clip = FALSE, distributeGaps = TRUE,
          removeError = TRUE, topx = 0, outdir = outdir, force = force,
          clean = clean, indent = incr(indent2), ...)

        ## Construct consensus sequence
        conseqName <- "consensus.mapIter." %<<% iteration %<<% "." %<<% hptype
        conseqPath <- file.path(outdir, conseqName %<<% ".fa")
        flog.info("%sConstruct consensus <%s>", indent2(), self$relPath(conseqPath), name = "info")
        conseq <- conseq(pileup, name = conseqName, type = "prob", suppressInsGaps = TRUE, columnOccupancy = columnOccupancy)
        Biostrings::writeXStringSet(
          Biostrings::DNAStringSet(gsub("[-+]", "", conseq)),
          conseqPath)

        ## Initialize structure
        self$mapIter[[iterationC]][[hptype]] = structure(
          list(
            dir     = self$relPath(outdir),
            reads   = self$relPath(readfile),
            ref     = "mapIter" %<<% iteration,
            bamfile = self$relPath(path(pileup)),
            pileup  = pileup,
            conseq  = conseq,
            seqpath = self$relPath(conseqPath),
            params  = list(columnOccupancy = columnOccupancy),
            tag     = maptag
          ),
          class = c("mapIter", "list")
        )
      }
    }

    if (plot) {
      flog.info("%sPlot MapIter summary", indent(), name = "info")
      ## Coverage and base frequency
      plotlist <- foreach(iteration = seq_len(self$getIterations())) %do% {
        suppressWarnings(self$plotmapIterSummary(thin = 0.1, width = 4,
                         iteration = iteration, drop.indels = TRUE))
      }
      p <- cowplot::plot_grid(plotlist = plotlist, nrow = self$getIterations())
      cowplot::save_plot(p, filename = self$absPath("plot.MapIter.pdf"),
                         base_width = 12*length(hptypes),
                         title = paste(self$getLocus(),
                                       self$getSampleId(), sep = "." ),
                         base_height = 3*self$getIterations())
      cowplot::save_plot(p, filename = self$absPath(".plots/plot.MapIter.svg"),
                         base_width = 12*length(hptypes),
                         base_height = 3*self$getIterations())

    }

    createIgvConfigs(self, map = "mapIter", open = "FALSE")

    ## set mapIter runstats
    .setRunstats(self, "mapIter",
                 list(Runtime = format(Sys.time() - start.time)))

    return(invisible(self))
  })

#' @export
print.mapIter <- function(x, ...) {
  #x <- self$mapIter$`0`$A
  msg <- sprintf("An object of class '%s'.\n", class(x)[1])
  bamf <- ifelse(is.null(x$bamfile), " no bamfile", basename(x$bamfile %||% ""))
  msg <- sprintf("%s [Dir] %s\n", msg, x$dir)
  readmsg <- if (is.null(x$reads)) "" else sprintf(" [Reads] %s\n ",
                                                   basename(x$reads))
  refmsg <- if (is.null(x$ref)) "" else sprintf("[Reference] %s\n ",
                                                basename(x$ref))
  cat(msg%<<%readmsg%<<%refmsg)
}

## Method: partitionShortReads ####
#' @export
partitionShortReads.DR2S <- function(x,
                                     opts = list(),
                                     force = FALSE,
                                     ...) {
  x$runPartitionShortReads(opts = opts,
                           force = force,
                           ...)
  invisible(x)
}
# TODO: look at arguments and make same
DR2S_$set("public", "runPartitionShortReads", function(opts = list(),
                                                       force = FALSE,
                                                       ...) {
  ## debug
  # opts = list()
  # force = FALSE

  flog.info("# PartitionShortReads ...", name = "info")

  ## Collect start time for partitionShortReads runstats
  start.time <- Sys.time()

  ## Initiate indenter
  indent <- indentation(1)

  ## exit savely if shortreads not provided
  if (!self$hasShortreads()) {
    flog.warn("%sCannot partition shortreads. No shortreads provided", indent(), name = "info")
    return(invisible(self))
  }

  ## exit savely if initial SR mapping not performed
  if (is.null(self$mapInit$SR2)) {
    flog.warn("%sCannot partition shortreads. Run 'mapInit()' first", indent(), name = "info")
    return(invisible(self))
  }

  ## exit safely if reporting is already finished
  if (.checkReportStatus(self)) return(invisible(self))

  flog.info("%sPartition shortreads based on initial mapping and " %<<%
            "longread clustering", indent(), name = "info")

  ## Overide default arguments
  args <- self$getOpts("partitionSR")
  if (!is.null(args)) {
    env  <- environment()
    list2env(args, envir = env)
  }

  bamfile <- self$absPath(self$mapInit$SR2$bamfile)
  hptypes <- self$getHapTypes()
  prtMat  <- self$partition$mat
  seqs <- lapply(self$partition$hpl, function(x) .getSeqsFromMat(
    as.matrix(prtMat[x, ])))
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
  # hptype <- "A"
  indent2 <- incr(indent)
  foreach(hptype = hptypes) %do% {
    srfilenames <- c()
    flog.info("%sWrite shortread fastq for haplotype <%s>", indent2(), hptype, name = "info")
    fqs <- self$getShortreads()
    # dontUseReads <- srpartition$haplotypes$read[
    #   !srpartition$haplotypes$read %in% dplyr::filter(
    #     srpartition$haplotypes, haplotype == hptype)$read]
    dontUseReads <- dplyr::filter(srpartition$haplotypes, haplotype != hptype)$read
    # write fastq's
    # fq <- fqs[1]
    foreach(fq = fqs) %do% {
      srFastqHap <- self$absPath(
        file.path(
          self$mapIter$`0`[[hptype]]$dir,
          dot(c(strsplit1(basename(fq), "\\.")[1], hptype, "fastq.gz"))))
      .writePartFq(fq = fq, srFastqHap = srFastqHap, dontUseReads = dontUseReads, indent = incr(indent2))
      srfilenames <- c(srfilenames, srFastqHap)
      self$srpartition[[hptype]]$srpartition <- srpartition
      NULL
    }
    self$srpartition[[hptype]]$SR <- self$relPath(srfilenames)
    NULL
  }

  ## set partitionShortReads runstats
  ## set mapIter runstats
  .setRunstats(self, "partitionShortReads",
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
  flog.info("%sMap shortreads and longreads against mapIter consensus sequences", indent(), name = "info")

  ## Check if reporting is already finished and exit safely
  if (!force && .checkReportStatus(self)) return(invisible(self))

  ## Overide default arguments
  args <- self$getOpts("mapFinal")
  if (!is.null(args)) {
    env  <- environment()
    list2env(args, envir = env)
  }

  ## stop if no shortreads provided
  if (!self$hasShortreads()) {
    flog.warn("%sCannot run mapFinal. No shortreads provided", indent(), name = "info")
    return(invisible(self))
  }

  reftag <- "mapFinal"
  outdir <- .dirCreateIfNotExists(self$absPath(reftag))
  lastIter <- self$mapIter[[max(names(self$mapIter))]]
  hptypes  <- self$getHapTypes()
  reffiles <- set_names(lapply(hptypes, function(x)
    self$absPath(lastIter[[x]]$seqpath)), hptypes)
  readfilesLR <- set_names(lapply(hptypes, function(x)
    self$absPath(lastIter[[x]]$reads)), hptypes)
  readfilesSR <- set_names(lapply(hptypes, function(x)
    self$absPath(self$srpartition[[x]]$SR)), hptypes)

  self$mapFinal = structure(
    list(
      dir          = self$relPath(outdir),
      sreads       = lapply(readfilesSR, self$relPath),
      lreads       = lapply(readfilesLR, self$relPath),
      ref          = lapply(reffiles, self$relPath),
      bamfile      = list(),
      pileup       = list(),
      tag          = list(),
      seqpath      = list(),
      homopolymers = NULL
    ), class = c("mapFinal", "list")
  )

  ## hptype = "A"
  ## hptype = "B"
  indent2 <- incr(indent)
  for (hptype in hptypes) {
    flog.info("%sFor haplotype <%s>", indent(), hptype, name = "info" )
    reffile <- reffiles[[hptype]]
    ## Map longreads
    flog.info("%sMap longreads", indent2(), name = "info")
    readtype <- self$getLrdType()
    mapgroup <- "LR" %<<% hptype
    maptag   <- paste("mapFinal", mapgroup, readtype, self$getLrdMapper(), optstring(opts), sep = ".")
    readfile <- readfilesLR[[hptype]]
    pileup <- mapReads(
      mapFun = self$getLrMapFun(), maptag = maptag, reffile = reffile,
      allele = mapgroup, readfile = readfile, readtype = readtype,
      opts = opts, refname = hptype, includeDeletions = includeDeletions,
      includeInsertions = includeInsertions, callInsertions = FALSE,
      clip = FALSE, distributeGaps = TRUE, removeError = TRUE, topx = 0,
      outdir = outdir, force = force, clean = TRUE, max_depth = 1e4,
      min_mapq = 50, indent = incr(indent2), ...)
    self$mapFinal$bamfile[[mapgroup]] = self$relPath(path(pileup))
    self$mapFinal$pileup[[mapgroup]] = pileup
    self$mapFinal$tag[[mapgroup]] = maptag
    if (createIgv) {
      self$mapFinal$igv[[mapgroup]] = createIgvJsFiles(
        refpath(pileup), path(pileup), self$getOutdir(), sampleSize = 100,
        fragmentReads = TRUE)
    }

    ## Map shortreads
    flog.info("%sMap shortreads", indent2(), name = "info")
    readtype <- self$getSrdType()
    mapgroup <- "SR" %<<% hptype
    maptag   <- paste("mapFinal", mapgroup, readtype, self$getSrdMapper(), optstring(opts), sep = ".")
    readfile <- readfilesSR[[hptype]]
    pileup <- mapReads(
      mapFun = self$getSrMapFun(), maptag = maptag, reffile = reffile,
      allele = mapgroup, readfile = readfile, readtype = readtype,
      opts = opts, refname = hptype, includeDeletions = includeDeletions,
      includeInsertions = includeInsertions, callInsertions = TRUE,
      clip = FALSE, distributeGaps = TRUE, removeError = TRUE, topx = 0,
      outdir = outdir, force = force, clean = TRUE, max_depth = 1e5,
      min_mapq = 50, min_base_quality = 13, indent = incr(indent2), ...)
    self$mapFinal$bamfile[[mapgroup]] = self$relPath(path(pileup))
    self$mapFinal$pileup[[mapgroup]] = pileup
    self$mapFinal$tag[[mapgroup]] = maptag
    if (createIgv) {
      self$mapFinal$igv[[mapgroup]] = createIgvJsFiles(
        refpath(pileup), path(pileup), self$getOutdir(), sampleSize = 100)
    }

    ## Construct consensus sequence
    conseqName <- "consensus.mapFinal." %<<% hptype
    conseqPath <- file.path(outdir, conseqName %<<% ".fa")
    flog.info("%sConstruct consensus <%s>", indent2(), self$relPath(conseqPath), name = "info")
    conseq <- conseq(pileup, name = conseqName, type = "ambig", suppressAllGaps = TRUE, threshold = 1/4)
    Biostrings::writeXStringSet(
      Biostrings::DNAStringSet(gsub("[-+]", "", conseq)),
      conseqPath)

    self$mapFinal$seq[[hptype]] <- conseq
    self$mapFinal$seqpath[[hptype]] <- self$relPath(conseqPath)
  }

  if (plot) {
    flog.info("%sPlot MapFinal summary", indent(), name = "info")
    ## Coverage and base frequency
    if (!is.null(self$mapFinal$sreads$A)) {
      readtypes <- c("LR", "SR")
    } else {
      readtypes <- c("LR")
    }
    plotlist <- foreach(readtype = readtypes) %do% {
      suppressWarnings(self$plotmapFinalSummary(iteration = "final",
                       readtype = readtype, thin = 0.25, width = 20))
    }
    p <- cowplot::plot_grid(plotlist = plotlist, nrow = 2, labels = readtypes)
    cowplot::save_plot(p, filename = self$absPath("plot.MapFinal.pdf"),
                       base_width = 12*length(hptypes),
                       title = paste(self$getLocus(), self$getSampleId(),
                                         sep = "." ),
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

#' @export
print.mapFinal <- function(x, ...) {
  fmt <- "%s [Dir] %s\n [Longreads] %s\n [Shortreads] %s\n [References] %s\n [Bamfile] %s"
  msg  <- sprintf("An object of class '%s'.\n", class(x)[1])
  bamf <- comma(basename(unlist(x$bamfile) %||% ""))
  msg <- sprintf(fmt, msg, x$dir,
    comma(basename(unlist(x$lreads))),
    comma(basename(unlist(x$sreads))),
    comma(basename(unlist(x$ref))),
    bamf
  )
  cat(msg)
}

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

