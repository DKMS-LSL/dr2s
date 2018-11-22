#' @export
InitDR2S.DR2Sconf <- function(config, createOutdir = TRUE) {
  DR2S_$new(config, createOutdir)
}

#' @export
cache.DR2S <- function(x, outname, ...) {
  if (missing(outname)) {
    outname <- paste("DR2S", x$getLrdType(), x$getLrdMapper(), "rds", sep = ".")
  }
  x$cache(outname = outname)
  invisible(x)
}

#' @export
clear.DR2S <- function(x, ...) {
  x$clear()
  invisible(x)
}


# Class: DR2S -------------------------------------------------------------


#' Class \code{"DR2S"}
#'
#' @docType class
#' @usage InitDR2S(config, createOutdir = TRUE)
#' @field InitDR2S Initialize DR2S from a config.
#' @field mapInit \code{[MapList]}; the mapping of long reads to
#' \code{reference}.
#' @field partition \code{[PartList]}; the partitioning of full-length
#'   mapped long reads into different haplotypes.
#' @field mapIter \code{[MapList]}; the mapping of A and B reads to
#' consensus sequences produced in the previous step.
#' @field mapFinal \code{[MapList]}; the mapping of A, B, and short reads
#' to consensus sequences produced in the previous step.
#' @field consensus \code{[ConsList]}; final consensus sequences for A and
#' B.
#' @keywords data internal
#' @return Object of \code{\link[R6]{R6Class}} representing a DR2S analysis.
#' @section Public Methods:
#' \describe{
#' \item{\code{x$runMapInit( opts = list(), createIgv = TRUE, plot = TRUE, ...)}}
#' {Run the inital mapping step (long reads against the reference allele)}
#' \item{\code{x$runPartitionLongreadsplot = TRUE, ...)}}
#' {Partition mapped longreads into haplotypes}
#' \item{\code{x$runSplitLongreadsByHaplotype((plot = TRUE)}}{
#' Partition mapped longreads into haplotypes}
#' \item{\code{x$runExtractPartitionedLongreads()}}{
#' Extract FASTQs for partitioned reads.
#' }}
#' @section Internals:
#' \describe{
#' \item{\code{x$new(conf, createOutdir = TRUE)}}{Initialise object of class
#'  \code{DR2S}.}
#' \item{x$getConfig(name = NULL}{Extract configuration}
#' \item{x$setConfig(name, value}{Modify configuration}
#' }
DR2S_ <- R6::R6Class(
  classname = "DR2S",
  public = list(
    ## Public fields
    mapInit     = list(), ## <MapList>;
    lrpartition = list(), ## <PartList>;
    mapIter     = list(), ## <MapList>;
    srpartition = list(), ## <PartList>;
    mapFinal    = list(), ## <MapList>;
    consensus   = list(), ## <ConsList>;
    ## Public functions
    initialize  = function(conf, createOutdir = TRUE) {
      private$conf = initialiseDR2S(conf, createOutdir = createOutdir)
      private$runstats = NULL
      private$reportStatus = FALSE
      if (file.exists(refPath <- private$conf$reference)) {
        private$conf$extref = normalizePath(refPath, mustWork = TRUE)
        private$conf$reference = sub("\\.fa(s|sta)?$", "", basename(refPath))
        private$runstats$refPath = file.path("mapInit", basename(refPath))
        cPath <- .dirCreateIfNotExists(normalizePath(
          file.path(private$conf$outdir, "mapInit"), mustWork = FALSE))
        file.copy(refPath, file.path(cPath, basename(refPath)))
      } else {
        private$conf$reference = .expandAllele(conf$reference, conf$locus)
        private$runstats$refPath = .generateReferenceSequence(
          private$conf$reference,
          private$conf$locus,
          private$conf$outdir,
          "mapInit")
      }
      .confLog(outdir = private$conf$outdir, logName = "info")
      private$srStatus = self$checkGetShortreads()
      private$lrStatus = self$checkGetLongreads()
      writeDR2SConf(self)
      flog.info("Creating DR2S Object", name = "info")
    },
    cache = function(outname) {
      if (missing(outname)) {
        outname <- paste("DR2S", self$getLrdType(), self$getLrdMapper(), "rds",
                         sep = ".")
      }
      path <- file.path(self$getOutdir(), outname)
      flog.info("Caching <%s>", path, name = "info")
      saveRDS(self$clone(deep = TRUE), file = path, compress = TRUE)
      invisible(self)
    },
    clear = function() {
      unlink(self$getOutdir(), recursive = TRUE)
      .dirCreateIfNotExists(file.path(self$getOutdir()))
      if (!is.null(private$conf$extref) &&
          file.exists(refPath <- private$conf$extref)) {
        private$conf$reference = sub("\\.fa(s|sta)?$", "", basename(refPath))
        private$runstats$refPath = file.path("mapInit", basename(refPath))
        cPath <- .dirCreateIfNotExists(normalizePath(
          file.path(private$conf$outdir, "mapInit"), mustWork = FALSE))
        file.copy(refPath, file.path(cPath, basename(refPath)))
      } else {
        private$conf$reference = .expandAllele(private$conf$reference,
                                               private$conf$locus)
        outdir <- .dirCreateIfNotExists(normalizePath(
          file.path(private$conf$outdir, "mapInit"), mustWork = FALSE))
        private$runstats$refPath = .generateReferenceSequence(
          self$getReference(),
          self$getLocus(),
          private$conf$outdir,
          "mapInit")
      }
      .confLog(outdir = private$conf$outdir, logName = "info")
      private$srStatus = self$checkGetShortreads()
      private$lrStatus = self$checkGetLongreads()
      writeDR2SConf(self)
      invisible(self)
    },
    cleanup = function() {
      unlink(self$absPath(bampath(self$mapInit)))
      unlink(self$absPath(bampath(self$mapInit)) %<<% ".bai")
      foreach(hpt = self$getHapTypes()) %do% {
        unlink(self$absPath(hpt), recursive = TRUE)
      }
      unlink(self$absPath("final"), recursive = TRUE)
      invisible(self)
    },
    print = function() {
      fmt0 <- "DR2S driver for sample <%s> locus <%s>\n"
      cat(sprintf(fmt0, self$getSampleId(), self$getLocus()))
      fmt1 <- "Reference: <%s>\n" %<<%
              "Longreads: <%s> Shortreads: <%s>\n" %<<%
              "Mapper: <%s>\nDatadir: <%s>\nOutdir: <%s>\n"
      cat(sprintf(fmt1,
                  self$getReference(), self$getLrdType(),
                  self$getSrdType() %||% "", self$getLrdMapper(),
                  self$getDatadir(), self$getOutdir()
      ))
      invisible(self)
    },
    run_ = function(step) {
      switch(
        step,
        clear    = self$clear(),
        cache    = self$cache(),
        partitionLongreads = {
          self$runPartitionLongreads()
          self$runSplitLongreadsByHaplotype()
          self$runExtractPartitionedLongreads()
        },
        partitionShortreads = {
          self$runPartitionShortreads()
        },
        mapInit  = self$runMapInit(),
        mapIter  = self$runMapIter(),
        mapFinal = self$runMapFinal(),
        polish   = DR2S::polish(self),
        report   = DR2S::report(self),
        stop("<", step, "> is not a valid step in the mapping pipeline")
      )
    },
    ##
    ## Config getters and setters ####
    ##
    getConfig = function(name = NULL) {
      if (is.null(name))
        private$conf
      else
        private$conf[[name]]
    },
    ##
    setConfig = function(name, value) {
      private$conf[name] = value
      invisible(self)
    },
    ##
    getDatadir = function() {
      self$getConfig("datadir")
    },
    ##
    getOutdir = function() {
      self$getConfig("outdir")
    },
    ##
    getThreshold = function() {
      self$getConfig("threshold")
    },
    ##
    setThreshold = function(threshold) {
      assert_that(threshold > 0 && threshold <= 1)
      self$setConfig("threshold", threshold)
      invisible(self)
    },
    ##
    getIterations = function() {
      self$getConfig("iterations")
    },
    ##
    setIterations = function(iterations) {
      assert_that(iterations < 10 && iterations > 0 && iterations %% 1 == 0)
      self$setConfig("iterations", iterations)
      invisible(self)
    },
    ##
    getMicrosatellite = function() {
      self$getConfig("microsatellite")
    },
    ##
    setMicrosatellite = function(microsatellite) {
      assert_that(is.logical(microsatellite))
      self$setConfig("microsatellite", microsatellite)
      invisible(self)
    },
    ##
    getFilterScores = function() {
      self$getConfig("filterScores")
    },
    ##
    setFilterScores = function(filterScores) {
      assert_that(is.logical(filterScores))
      self$setConfig("filterScores", filterScores)
      invisible(self)
    },
    ##
    getForceMapping = function() {
      self$getConfig("forceMapping")
    },
    ##
    setForceMapping = function(forceMapping) {
      assert_that(is.logical(forceMapping))
      self$setConfig("forceMapping", forceMapping)
      invisible(self)
    },
    ##
    getLrdMapper = function() {
      self$getConfig("longreads")$mapper
    },
    ##
    getLrdType = function() {
      self$getConfig("longreads")$type
    },
    ##
    getLrdDir = function() {
      lrddir <- self$getConfig("longreads")$dir
      if (is.null(lrddir)) {
        self$getDatadir()
      } else {
        file.path(self$getDatadir(), lrddir)
      }
    },
    ##
    getLongreads = function() {
      lrdfile <- self$getConfig("longreads")$file
      if (!is.null(lrdfile)) {
        readpath <- file.path(self$getDatadir(), lrdfile)
        names(readpath) <- lrdfile
      } else {
        lrddir <- self$getLrdDir()
        readpath <- findReads(datadir = lrddir,
                              sampleId = self$getSampleId(),
                              locus = self$getLocus())
        names(readpath) <- .cropPath(self$getDatadir(), readpath)
      }
      if (is.null(readpath) || length(readpath) == 0 || !file.exists(readpath)) {
        flog.error("No reads available for readtype <%s>", self$getLrdType(), name = "info")
        stop("No reads available for readtype <", self$getLrdType(), ">")
      }
      readpath
    },
    ##
    getSrdMapper = function() {
      self$getConfig("shortreads")$mapper
    },
    ##
    getSrdType = function() {
      self$getConfig("shortreads")$type
    },
    ##
    getSrdDir = function() {
      filtered <- self$getConfig("filteredShortreads")
      if (!is.null(filtered)) {
        return(self$absPath(filtered))
      }
      srddir <- self$getConfig("shortreads")$dir
      if (!is.null(srddir)) {
        file.path(self$getDatadir(), srddir)
      } else {
        NULL
      }
    },
    ##
    getShortreads = function() {
      srddir <- self$getSrdDir()
      if (is.null(srddir)) {
        return(NULL)
      }
      readpath <- findReads(srddir, self$getSampleId(), self$getLocus())
      names(readpath) <- .cropPath(self$getDatadir(), readpath)
      if (is.null(readpath) || length(readpath) == 0 || !file.exists(readpath)) {
        flog.error("No reads available for readtype <%s>", self$getSrdType(), name = "info")
        stop("No reads available for readtype <", self$getSrdType(), ">")
      }
      readpath
    },
    ##
    getPipeline = function() {
      self$getConfig("pipeline")
    },
    ##
    getOpts = function(name = NULL) {
      if (is.null(name))
        mergeList(self$getConfig("opts"), self$getConfig("longreads")$opts, update = TRUE)
      else {
        opts <- mergeList(self$getConfig("opts"),
                          self$getConfig("longreads")$opts, update = TRUE)
        if (is.list(opts))
          opts[[name]]
        else
          opts
      }
    },
    ##
    getSampleId = function() {
      self$getConfig("sampleId")
    },
    ##
    getLocus = function() {
      self$getConfig("locus")
    },
    ##
    getReference = function() {
      self$getConfig("reference")
    },
    ##
    getPlatform = function() {
      self$getDetails("platform") %||%
        self$getConfig("longreads")$platform  %||%
        self$getConfig("longreads")$type
    },
    ##
    getDistAlleles = function() {
      self$getConfig("distAlleles")
    },
    ##
    setDistAlleles = function(distAlleles) {
      assert_that(is.numeric(distAlleles))
      self$setConfig("distAlleles", distAlleles)
      invisible(self)
    },
    ##
    getFormat = function() {
      self$getConfig("format")
    },
    ##
    setFormat = function(format) {
      assert_that(is.character(format) ,format %in% c("yaml", "json"))
      self$setConfig("format", format)
      invisible(self)
    },
    ##
    getDetails = function(name = NULL) {
      if (is.null(name)) {
        self$getConfig("details")
      } else {
        assert_that(is.character(name))
        self$getConfig("details")[[name]]
      }
    },
    ##
    setDetails = function(details) {
      assert_that(is.character(details))
      self$setConfig("details", details)
      invisible(self)
    },
    ##
    getSampleDetails = function() {
      details <- semicolon(
        vapply(seq_along(self$getDetails()),
               function(item) paste(names(
                 self$getDetails()[item]),
                 underscore(litQuote(self$getDetails()[item])),
                 sep = "="), FUN.VALUE = character(1)))
      sr <- self$hasShortreads()
      lr <- self$hasLongreads()
      paste("locus=" %<<% litQuote(self$getLocus()),
            "ref=" %<<% litQuote(self$getReference()),
            details,
            "short_read_data=" %<<% litQuote(ifelse(sr, "yes", "no")),
            "short_read_type=" %<<% litQuote(ifelse(sr, self$getSrdType(), "")),
            "long_read_data=" %<<% litQuote(ifelse(lr, "yes", "no")),
            "long_read_type=" %<<% litQuote(ifelse(lr, self$getLrdType(), "")),
            "software=\"DR2S\"",
            "version=" %<<% litQuote(packageVersion("DR2S")),
            sep = ";")
    },
    ##
    ## Runstats getters and setters ####
    ##
    getStats = function(name = NULL) {
      if (is.null(name))
        private$runstats
      else
        private$runstats[[name]]
    },
    ##
    setStats = function(name, value) {
      private$runstats[name] = value
      invisible(self)
    },
    ##
    getLimits = function() {
      self$getStats("partitionLongreads")$limits
    },
    ##
    setLimits = function(lmts) {
      # stopifnot(is.null(lmts) || all(lmts >=0 && lmts <=1))
      x <- list()
      private$runstats$partitionLongreads$limits = lmts
      invisible(self)
    },
    ##
    getHapTypes = function() {
      self$getStats("partitionLongreads")$haplotypes
    },
    ##
    setHapTypes = function(x) {
      ## Todo: add test for data consistency!
      private$runstats$partitionLongreads$haplotypes = x
      invisible(self)
    },
    ##
    getRefPath = function() {
      self$absPath(self$getStats("refPath"))
    },
    ##
    setRefPath = function(refPath) {
      self$setStats("refPath", self$relPath(refPath))
    },
    ##
    ## Report Status ####
    ##
    getReportStatus = function() {
      private$reportStatus
    },
    ##
    setReportStatus = function(report) {
      assert_that(is.logical(report))
      private$reportStatus <- report
      invisible(self)
    },
    ##
    ## Internal getters/setters ####
    ##
    getRefSeq = function() {
      if (!is.null(self$getRefPath())) {
        Biostrings::readDNAStringSet(self$getRefPath())
      } else  {
        if (startsWith(self$getLocus, "KIR")) {
          .ipdKir()$getClosestComplete(self$getReference(),
                                       self$getLocus())
        } else {
          .ipdHla()$getClosestComplete(self$getReference(),
                                       self$getLocus())
        }
      }
    },
    ##
    getPileup = function(useSR = FALSE) {
      if (useSR)
        meta(self$mapInit, "SR2")$pileup
      else
        self$mapInit$pileup
    },
    ##
    getSnpMatrix = function() {
      self$lrpartition$mat
    },
    ##
    getPartition = function() {
      self$lrpartition$prt
    },
    ##
    getHapList = function(group) {
      if (missing(group)) {
        self$lrpartition$hpl
      } else {
        self$lrpartition$hpl[[group]]
      }
    },
    ##
    getLatestRefPath = function() {
      ##
      if (self$hasMapFinal() && self$hasShortreads()) {
        return(vapply(self$mapFinal$SR, function(x) self$absPath(conspath(x)),
                      FUN.VALUE = character(1)))
      }
      ##
      if (self$hasMapFinal() && self$hasLongreads()) {
        return(vapply(self$mapFinal$LR, function(x) self$absPath(conspath(x)),
                      FUN.VALUE = character(1)))
      }
      ##
      if (self$hasMapIter()) {
        latest <- self$mapIter[[
          max(names(self$mapIter))]][self$getHapTypes()]
        return(vapply(latest, function(x) self$absPath(conspath(x)),
                      FUN.VALUE = character(1)))
      }
      ##
      return(self$getRefPath())
    },
    ##
    getLatestRef = function() {
      filepath <- self$getLatestRefPath()
      Biostrings::readDNAStringSet(filepath, use.names = TRUE)
    },
    ##
    getMapTag = function(iter = "init", group = NULL, ref = "LR") {
      if (is.numeric(iter)) {
        group <- match.arg(group, self$getHapTypes())
        return(tag(self$mapIter[[as.character(iter)]][[group]]))
      }
      if (iter == "init") {
        ref <- match.arg(ref, c("LR", "SR"))
        if (ref == "SR") {
          return(tag(meta(self$mapInit, "SR1")))
        }
        return(tag(self$mapInit))
      }
      if (iter == "final") {
        group <- match.arg(group, self$getHapTypes())
        ref <- match.arg(ref, c("LR", "SR"))
        return(tag(self$mapFinal[[ref]][[group]]))
      }
    },
    ##
    getGroupPileup = function(group = NULL, ref = NULL, iteration = 0) {
      group  <- match.arg(group, self$getHapTypes())
      ref    <- match.arg(ref, c("LR", "SR"))
      if (is.numeric(iteration)) {
        self$mapIter[[as.character(iteration)]][[group]]$pileup
      } else {
        self$mapFinal[[ref]][[group]]$pileup
      }
    },
    ##
    setConsensus = function(x) {
      stopifnot(is(x, "ConsList"))
      self$consensus = x
      invisible(self)
    },
    ##
    getLrdMapFun = function() {
      match.fun("run" %<<% self$getLrdMapper())
    },
    ##
    getSrdMapFun = function() {
      match.fun("run" %<<% self$getSrdMapper())
    },
    ## Get the absolut path
    absPath = function(filename) {
      assert_that(is.character(filename))
      ## set relpath as name attribute
      names(filename) <- filename
      vapply(filename, function(x) {
        normalizePath(
          file.path(self$getOutdir(), x),
          mustWork = FALSE)
      }, FUN.VALUE = character(1), USE.NAMES = TRUE)
    },
    ## Get the relative path
    relPath = function(filepath) {
      ## TODO add assertives
      if (all(startsWith(filepath, self$getOutdir()))) {
        .cropPath(base = self$getOutdir(), path = filepath)
      } else if (all(startsWith(filepath, self$getDatadir()))) {
        .cropPath(base = self$getDatadir(), path = filepath)
      } else {
        filepath
      }
    },
    ##
    ## Predicate methods ####
    ##
    hasMapInit = function() {
      tryCatch(
        is(self$mapInit, "MapList"),
        error = function(e)
          FALSE
      )
    },
    ##
    hasLrdPartition = function() {
      tryCatch(
        is(self$lrpartition$prt, "HapPart"),
        error = function(e)
          FALSE
      )
    },
    ##
    hasHapList = function() {
      tryCatch(
        is(self$lrpartition$hpl, "HapList"),
        error = function(e)
          FALSE
      )
    },
    ##
    hasMapIter = function() {
      tryCatch(
        all(vapply(self$mapIter[[self$getIterations() + 1]],
                   is, "MapList", FUN.VALUE = logical(1))),
        error = function(e)
          FALSE
      )
    },
    ##
    hasMapFinal = function() {
      tryCatch(
        all(vapply(self$mapFinal$LR, is, "MapList", FUN.VALUE = logical(1))),
        error = function(e)
          FALSE
      )
    },
    ##
    checkGetShortreads = function() {
      sr <- try(self$getShortreads(), silent = TRUE)
      !(is(sr, "try-error") || is.null(sr))
    },
    ##
    hasShortreads = function() {
      private$srStatus
    },
    ##
    checkGetLongreads = function() {
      lr <- try(self$getLongreads(), silent = TRUE)
      !(is(lr, "try-error") || is.null(lr))
    },
    ##
    hasLongreads = function() {
      private$lrStatus
    },
    ##
    ## Plotting methods ####
    ##
    ##
    plotCoverage = function(threshold,
                            range = NULL,
                            thin = 0.1,
                            width = 1,
                            label = "",
                            drop.indels = FALSE,
                            readtype = "LR") {
      tag <- self$getMapTag(ref = readtype)
      if (missing(threshold)) {
        threshold <- self$getThreshold()
      }
      pileup <- switch(
        readtype,
        LR = self$getPileup(useSR = FALSE),
        SR = self$getPileup(useSR = TRUE)
      )
      plotPileupCoverage(
        pileup,
        threshold,
        range,
        thin,
        width,
        label %|ch|% tag,
        drop.indels = drop.indels
      )
    },
    ##
    plotGroupCoverage = function(group = NULL,
                                 ref = NULL,
                                 iteration = "init",
                                 threshold = NULL,
                                 range = NULL,
                                 thin = 0.1,
                                 width = 1,
                                 label = "",
                                 drop.indels = FALSE) {
      pileup <- self$getGroupPileup(group, ref, iteration)
      lbl <- self$getMapTag(iteration, group)
      if (is.null(threshold)) {
        threshold <- self$getThreshold()
      }
      plotPileupCoverage(
        x = pileup,
        threshold = threshold,
        range = range,
        thin = thin,
        width = width,
        label = label %|ch|% lbl,
        drop.indels = drop.indels
      )
    },
    ##
    plotBasecallFrequency = function(threshold, label = "",
                                     drop.indels = FALSE) {
      tag <- self$getMapTag()
      if (missing(threshold)) {
        threshold <- self$getThreshold()
      }
      plotPileupBasecallFrequency(
        self$getPileup(),
        threshold,
        label %|ch|% tag,
        drop.indels = drop.indels
      )
    },
    ##
    plotGroupBasecallFrequency = function(group = NULL,
                                          ref = NULL,
                                          iteration = "init",
                                          threshold = NULL,
                                          label = "",
                                          drop.indels = FALSE) {
      pileup <- self$getGroupPileup(group, ref, iteration)
      lbl <- self$getMapTag(iteration, group)
      if (is.null(threshold)) {
        threshold <- self$getThreshold()
      }
      plotPileupBasecallFrequency(
        pileup,
        threshold,
        label %|ch|% lbl,
        drop.indels = drop.indels
      )
    },
    ##
    plotPartitionHistogram = function(label = "", limits = NULL) {
      tag <- self$getMapTag()
      if (length(self$getHapTypes()) == 2) {
        plotPartitionHistogram(x = self$getPartition(), label %|ch|% tag,
                                 limits = limits)
      } else {
        plotPartitionHistogramMulti(x = self$getPartition(),
                                       limits = limits,
                                       label %|ch|% tag)
      }
    },
    ##
    plotPartitionTree = function() {
      plotPartitionTree(x = self$getPartition())
    },
    ##
    plotPartitionRadar = function() {
      plotRadarPartition(x = self$getPartition())
    },
    ##
    ## Summary methods ####
    ##
    polymorphicPositions = function(threshold, useSR = FALSE) {
      if (missing(threshold)) {
        threshold <- self$getThreshold()
      }
      pileup <- self$getPileup(useSR = useSR)
      polymorphicPositions(pileup, threshold)
    },
    ##
    plotMapInitSummary = function(thin = 0.2, width = 4, label = "") {
      readtypes <- if (self$hasShortreads()) {
        c("SR", "LR")
      } else {
        c("LR")
      }
      plotlist <- foreach(readtype = readtypes) %do% {
        tag <- self$getMapTag(ref = readtype)
        self$plotCoverage(
          thin  = thin,
          width = width,
          label = label %|ch|% tag,
          readtype = readtype
        )
      }
      cowplot::plot_grid(plotlist = plotlist, labels = readtypes,
                         nrow = length(plotlist))
    },
    ##
    plotPartitionSummary = function(label = "", limits = NULL) {
      tag <- self$getMapTag()
      p1 <-  self$plotPartitionHistogram(label = label %||% tag,
                                         limits = limits) +
        ggplot2::theme(legend.position = "none")
      p2 <- self$plotPartitionTree()
      p3 <- self$plotPartitionRadar()
      if (is.null(p2)) {
        cowplot::plot_grid(p1, p3, ncol = 1)
      } else {
        cowplot::plot_grid(p1, p2, p3, ncol = 1, rel_heights = c(1,2,1))
      }
    },
    ##
    plotMapIterSummary = function(thin = 0.2, width = 10, iteration = 0,
                                  drop.indels = TRUE) {
      hptypes <- self$getHapTypes()
      if (iteration == self$getIterations()) {
        drop.indels <- FALSE
        width <- 10
      }
      plotlist <- foreach(hptype = hptypes) %do% {
        tag <- self$getMapTag(iteration, hptype)
        self$plotGroupCoverage(hptype, ref = NULL, iteration = iteration,
                               threshold = NULL, range = NULL, thin, width,
                               label = tag, drop.indels = drop.indels)
      }
      cowplot::plot_grid(plotlist = plotlist, labels = hptypes)
    },
    ##
    plotMapFinalSummary = function(readtype, thin = 0.2, width = 10,
                                   iteration = "final") {
      hptypes <- self$getHapTypes()
      ## hp = "A"
      plotlist <- foreach(hp = hptypes) %do% {
        tag <- self$getMapTag(iter = iteration, group = hp, ref = readtype)
        ref <- readtype %<<% hp
        self$plotGroupCoverage(group = hp, ref = readtype,
                               iteration = iteration, threshold = NULL,
                               range = NULL, thin = thin,
                               width = width, label = tag)
      }

      cowplot::plot_grid(plotlist = plotlist, labels = hptypes)
    },
    ##
    plotSeqLogo = function(ppos = NULL, pwm = NULL) {
      if (is.null(ppos)) {
        ppos <- SNP(self$getPartition())
        names(ppos) <- seq_along(ppos)
      }
      if (is.null(pwm)) {
        pwm <- lapply(PWM(self$getPartition()), function(pwm) {
          pwm[pwm < 0.1] <- 0
          pwm
        })
      }
      ggplot() +
        scale_x_continuous(labels = ppos, breaks = seq_along(ppos)) +
        geom_logo(pwm, method = "bits", seq_type = "dna", stack_width = 0.9) +
        facet_wrap(~seq_group, ncol = 1, strip.position = "left") +
        theme_logo() +
        theme(axis.text.x  = ggplot2::element_text(size = 10, angle = 60),
              axis.title.y = ggplot2::element_blank(),
              axis.text.y  = ggplot2::element_blank(),
              axis.ticks.y = ggplot2::element_blank(),
              strip.text.y = ggplot2::element_text(face = "bold",
                                                   size = 42, angle = 180))
    }
  ),
  ##
  ## Private fields ####
  ##
  private = list(
    conf         = NULL,
    runstats     = NULL,
    srStatus     = FALSE,
    lrStatus     = FALSE,
    reportStatus = FALSE
  )
)


# Helpers -----------------------------------------------------------------

findReads <- function(datadir, sampleId, locus) {
  locus <- sub("^HLA-", "", toupper(locus))
  locus <- sub("^KIR-", "", toupper(locus))
  filePattern <- sampleId %<<% "_" %<<% locus %<<% ".+" %<<% "fast(q|a)(\\.gz)?$"
  readPath <- dir(datadir, pattern = filePattern, full.names = TRUE)
  readPath <- readPath[grep(pattern = "^((?!_trimmed.fastq).)*$", readPath, perl = TRUE)]
  if (length(readPath) > 0) {
    normalizePath(readPath, mustWork = TRUE)
  } else readPath
}
