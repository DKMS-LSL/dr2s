#' @export
InitDR2S.DR2Sconf <- function(config, createOutdir = TRUE) {
  DR2S_$new(config, createOutdir)
}

#' @export
cache.DR2S <- function(x, outname, ...) {
  if (missing(outname)) {
    outname <- paste("DR2S", x$getLrdType(), x$getLrMapper(), "rds", sep = ".")
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
#' @field mapInit \code{[mapInit]}; the mapping of long reads to 
#' \code{reference}.
#' @field partition \code{[PartList]}; the partitioning of full-length 
#'   mapped long reads into different haplotypes.
#' @field mapIter \code{[mapIter]}; the mapping of A and B reads to 
#' consensus sequences produced in the previous step.
#' @field mapFinal \code{[mapFinal]}; the mapping of A, B, and short reads 
#' to consensus sequences produced in the previous step.
#' @field consensus \code{[ConsList]}; final consensus sequences for A and 
#' B.
#' @keywords data internal
#' @return Object of \code{\link[R6]{R6Class}} representing a DR2S analysis.
#' @section Public Methods:
#' \describe{
#' \item{\code{x$runMapInit(opts = list(), optsname = "",  threshold = 0.20, 
#' iterations = 1, microsatellite = FALSE, distAlleles = 2, filterScores = TRUE,
#' partSR = TRUE, minBaseQuality = 3, minMapq = 0, maxDepth = 1e4,
#' minNucleotideDepth = 3, includeDeletions = FALSE, includeInsertions = FALSE,
#' force = FALSE, fullname = TRUE, plot = TRUE)}}{Run the inital mapping step 
#' (long reads against the reference allele)}
#' \item{\code{x$runHaplotypePartitioning(maxDepth = 1e4, shuffle = TRUE,
#' skipGapFreq = 2/3, plot = TRUE)}}{
#' Partition mapped long reads into haplotypes}
#' \item{\code{x$splitReadsByHaplotype(limit = list(), #' plot = TRUE)}}{
#' Partition mapped long reads into haplotypes}
#' \item{\code{x$extractFastq()}}{
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
    mapInit     = list(),
    partition   = list(),
    mapIter     = list(),
    srpartition = list(),
    mapFinal    = list(),
    ## Final consensus sequences for A and B
    ## class: ConsList
    consensus   = list(),
    initialize  = function(conf, createOutdir = TRUE) {
      ## Public fields
      self$mapInit   = list()
      self$partition = list()
      self$mapIter   = list()#A = list(), B = list())
      self$mapFinal  = list()
      self$consensus = list()
      ## Private fields
      private$conf   = initialiseDR2S(conf, createOutdir = createOutdir)
      private$reportStatus = FALSE

      if (file.exists(refPath <- private$conf$reference)) {
        private$conf$reference = basename(refPath)
        private$conf$refPath  = basename(refPath)
      } else {
        private$conf$reference = .expandAllele(conf$reference, conf$locus)
        private$conf$refPath  = generateReferenceSequence(
          private$conf$reference,
          private$conf$locus,
          private$conf$outdir,
          "mapInit",
          fullname = FALSE)
      }
      .confLog(outdir = private$conf$outdir, logName = "info")
      writeDR2SConf(self)
      flog.info("Creating DR2S Object", name = "info")
    },
    cache = function(outname) {
      if (missing(outname)) {
        outname <- paste("DR2S", self$getLrdType(), self$getLrMapper(), "rds", 
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

      if (file.exists(refPath <- private$conf$reference)) {
        private$conf$reference = basename(refPath)
        private$conf$refPath  = basename(refPath)
      } else {
        private$conf$reference = .expandAllele(private$conf$reference, 
                                               private$conf$locus)
        outdir <- .dirCreateIfNotExists(normalizePath(
          file.path(private$conf$outdir, "mapInit"), mustWork = FALSE))
        private$conf$refPath  = generateReferenceSequence(
          self$getReference(),
          self$getLocus(),
          private$conf$outdir,
          "mapInit",
          fullname = FALSE)
      }
      invisible(self)
    },
    cleanup = function() {
      unlink(self$absPath(self$mapInit$bamfile))
      unlink(paste0(self$absPath(self$mapInit$bamfile), ".bai"))
      foreach(hpt = self$getHapTypes()) %do% {
        unlink(self$absPath(hpt), recursive = TRUE)
      }
      unlink(self$absPath("final"), recursive = TRUE)
      invisible(self)
    },
    print = function() {
      fmt0 <- "DR2S mapper for sample <%s> locus <%s>\n"
      cat(sprintf(fmt0, self$getSampleId(), self$getLocus()))
      fmt1 <- paste0("Reference alleles: <%s>\n", 
                     "Longreads: <%s> Shortreads: <%s>\n", 
                     "Mapper: <%s>\nDatadir: <%s>\nOutdir: <%s>\n")
      cat(sprintf(fmt1,
                  self$getReference(),
                  self$getLrdType(), self$getSrdType(), self$getLrMapper(), 
                  self$getDatadir(), self$getOutdir()
      ))
      invisible(self)
    },
    run_         = function(step) {
      switch(
        step,
        clear          = self$clear(),
        cache          = self$cache(),
        partitionLongReads         = {
          self$runPartitionLongReads()
          self$runSplitLongReadsByHaplotype()
          self$runExtractLongReads()
        },
        partitionShortReads         = self$runPartitionShortReads(),
        mapInit        = self$runMapInit(),
        mapIter        = self$runMapIter(),
        mapFinal       = self$runMapFinal(),
        polish         = DR2S::polish(self),
        report         = DR2S::report(self),
        stop("<", step, "> is not a valid step in the mapping pipeline")
      )
    },
    ##
    ## Getters and Setters ####
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
      stopifnot(threshold > 0 && threshold <= 1)
      self$setConfig("threshold", threshold)
      invisible(self)
    },
    ##
    getIterations = function() {
      self$getConfig("iterations")
    },
    ##
    setIterations = function(iterations) {
      stopifnot(iterations < 10 && iterations > 0 && iterations %% 1 == 0)
      self$setConfig("iterations", iterations)
      invisible(self)
    },
    ##
    getPartSR = function() {
      self$getConfig("partSR")
    },
    ##
    setPartSR = function(partSR) {
      stopifnot(is.logical(partSR))
      self$setConfig("partSR", partSR)
      invisible(self)
    },
    ##
    getReportStatus = function(){
      private$reportStatus
    },
    ##
    setReportStatus = function(report){
      stopifnot(is.logical(report))
      private$reportStatus <- report
      invisible(self)
    },
    ##
    getDetails = function(){
      self$getConfig("details")
    },
    ##
    setDetails = function(details){
      stopifnot(is.character(details))
      self$setConfig("details", details)
      invisible(self)
    },
    ##
    getForceMapping = function(forceMapping){
      self$getConfig("forceMapping")
    },
    ##
    setForceMapping = function(forceMapping){
      stopifnot(is.logical(forceMapping))
      self$setConfig("forceMapping", forceMapping)
      invisible(self)
    },
    ##
    getFilterScores = function() {
      self$getConfig("filterScores")
    },
    ##
    setFilterScores = function(filterScores) {
      stopifnot(is.logical(filterScores))
      self$setConfig("filterScores", filterScores)
      invisible(self)
    },
    ##
    getMicrosatellite = function() {
      self$getConfig("microsatellite")
    },
    ##
    setMicrosatellite = function(microsatellite) {
      stopifnot(is.logical(microsatellite))
      self$setConfig("microsatellite", microsatellite)
      invisible(self)
    },
    ##
    getSampleDetails = function() {
      details <- semicolon(
        sapply(seq_along(self$getDetails()), 
               function(item) paste(names(
                 self$getDetails()[item]), 
                 underscore(litQuote(self$getDetails()[item])), 
                 sep = "=")))
      sr <- !is(try(self$getShortreads(), silent = TRUE), "try-error")
      lr <- !is(try(self$getShortreads(), silent = TRUE), "try-error")
      paste(paste0("locus=", litQuote(self$getLocus())),
            paste0("ref=", litQuote(self$getReference())),
            details,
            paste0("short_read_data=", litQuote(ifelse(sr, "yes", "no"))),
            paste0("short_read_type=", litQuote(ifelse(sr, 
                                                       self$getSrdType(), 
                                                       ""))),
            paste0("long_read_data=", litQuote(ifelse(sr, "yes", "no"))),
            paste0("long_read_type=", litQuote(ifelse(lr, 
                                                       self$getLrdType(), 
                                                       ""))),
            "software=\"DR2S\"", 
            paste0("version=", 
                   litQuote(packageVersion("DR2S"))),
            sep = ";")
    },
    ##
    getDistAlleles = function() {
      self$getConfig("distAlleles")
    },
    ##
    setDistAlleles = function(distAlleles) {
      stopifnot(is.numeric(distAlleles))
      self$setConfig("distAlleles", distAlleles)
      invisible(self)
    },
    ##
    getLongreads = function() {
      dir <- self$getLrdDir()
      readpath <- findReads(dir, self$getSampleId(), self$getLocus())
      if (is.null(readpath) || length(readpath) == 0) {
        flog.error("No reads available for readtype <%s>", 
                   self$getLrdType(), name = "info")
        stop("No reads available for readtype <", self$getLrdType(), ">")
      }
      readpath
    },
    ##
    getShortreads = function() {
      dir <- self$getSrdDir()
      if (is.null(dir)) {
        return(NULL)
      }
      readpath <- findReads(dir, self$getSampleId(), self$getLocus())
      if (is.null(readpath) || length(readpath) == 0) {
        flog.error("No reads available for readtype <%s>", 
                   self$getSrdType(), name = "info")
        stop("No reads available for readtype <", self$getSrdType(), ">")
      }
      readpath
    },
    ##
    getLrMapper = function() {
      self$getConfig("lrmapper")
    },
    ##
    getSrMapper = function() {
      self$getConfig("srmapper")
    },
    ##
    getLimits = function() {
      self$getConfig("limits")
    },
    ##
    setLimits = function(lmts) {
      # stopifnot(is.null(lmts) || all(lmts >=0 && lmts <=1))
      private$conf$limits = lmts
      invisible(self)
    },
    ##
    getHapTypes = function() {
      self$getConfig("haptypes")
    },
    ##
    setHapTypes = function(x) {
      ## Todo: add test for data consistency!
      private$conf$haptypes = x
      invisible(self)
    },
    ##
    getPipeline = function() {
      self$getConfig("pipeline")
    },
    ##
    getLrdType = function() {
      self$getConfig("longreads")$type
    },
    ##
    getLrdDir = function() {
      file.path(self$getDatadir(), self$getConfig("longreads")$dir)
    },
    ##
    getSrdType = function() {
      self$getConfig("shortreads")$type
    },
    ##
    getSrdDir = function() {
      if (!is.null(filtered <- self$getConfig("filteredShortreads"))) {
        return(self$absPath(filtered))
      }
      if (!is.null(srddir <- self$getConfig("shortreads")$dir)) {
        file.path(self$getDatadir(), srddir)
      } else {
        NULL
      }
    },
    ##
    getReftype = function() {
      self$getConfig("reftype")
    },
    ##
    getReference = function() {
      self$getConfig("reference")
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
    getRefPath = function() {
      self$absPath(self$getConfig("refPath"))
    },
    ##
    getRefSeq = function() {
      if (!is.null(self$getRefPath())) {
        Biostrings::readDNAStringSet(self$getRefPath())
      } else  {
        ipd.Hsapiens.db::getClosestComplete(self$getReference(),
                                            self$getLocus())
      }
    },
    ##
    getOpts = function(name = NULL) {
      if (is.null(name))
        .mergeList(self$getConfig("opts"), self$getConfig("longreads")$opts, 
                   update = TRUE)
      else {
        opts <- .mergeList(self$getConfig("opts"), 
                           self$getConfig("longreads")$opts, update = TRUE)
        if (is.list(opts))
          opts[[name]]
        else
          opts
      }
    },
    ##
    getPileup = function() {
      self$mapInit$pileup
    },
    ##
    getSnpMatrix = function() {
      self$partition$mat
    },
    ##
    getPartition = function() {
      self$partition$prt
    },
    ##
    getHapList = function(group) {
      if (missing(group)) {
        self$partition$hpl
      } else {
        self$partition$hpl[[group]]
      }
    },
    ##
    getLatestRefPath = function() {
      if (length(self$mapFinal) > 0) {
        return(self$absPath(self$mapFinal$seq))
      } else if (length(self$mapIter) > 0) {
        latest <-  self$mapIter[max(names(self$mapIter))][self$getHapTypes()]
        return(sapply(latest, function(x) self$absPath(x$seqpath)))
      } else {
        return(self$getRefPath())
      }
    },
    ##
    getLatestRef = function() {
      if (!is.null(self$consensus$seq)) {
        return(self$consensus$seq)
      }else if (length(self$mapFinal) > 0) {
        return(self$mapFinal$seq)
      } else if (length(self$mapIter) > 0) {
        latest <- self$mapIter[max(names(self$mapIter))][self$getHapTypes()]
        return(sapply(latest, function(x) x$conseq))
      }  else {
        return(self$getRefSeq())
      }
    },
    ##
    getConseqs = function(group = NULL, mapn = 1) {
      group <- match.arg(group, self$getHapTypes)
      if (mapn == 1) {
        ref <- self$mapIter[["0"]][[group]]$conseq$reference %||% 
          Biostrings::BStringSet()
          Biostrings::BStringSet()
      } else if (mapn == 2) {
        ref <- Biostrings::BStringSet()
      }
      seqs <- tryCatch({
        Biostrings::DNAStringSet(ref)
      }, error = function(e) {
        warning("Possibly indels in consensus")
        Biostrings::DNAStringSet(gsub("[-+]", "", ref ))
      })
      metadata(seqs) <- compact(list(
        ref = tryCatch(
          metadata(ref)$zscore,
          error = function(e)
            NULL
        )
      ))
      seqs
    },
    ##
    getMapTag = function(iter = "init", group = NULL, ref = "LR") {

      if (is.numeric(iter)) {
        group <- match.arg(group, self$getHapTypes())
        return(self$mapIter[[as.character(iter)]][[group]]$tag)
      }
      if (iter == "init") {
        ref <- match.arg(ref, c("LR", "SR"))
        if (ref == "SR") {
          return(self$mapInit$SR1$tag)
        }
        return(self$mapInit$tag)
      }
      if (iter == "final") {
        group = match.arg(group, self$getHapTypes())
        ref = match.arg(ref, c("LR", "SR"))
        return(self$mapFinal$tag[[paste0(ref, group)]])
      }
    },
    ##
    getGroupPileup = function(group = NULL, ref = NULL, iteration = 0) {
      group  <- match.arg(group, self$getHapTypes())
      ref    <- match.arg(ref, c("LR", "SR"))
      if (is.numeric(iteration)) {
        self$mapIter[[as.character(iteration)]][[group]]$pileup
      } else {
        self$mapFinal$pileup[[paste0(ref, group)]]
      }
    },
    ##
    setConsensus = function(x) {
      stopifnot(is(x, "ConsList"))
      self$consensus = x
      invisible(self)
    },
    ##
    getLrMapFun = function() {
      match.fun(paste0("run", self$getLrMapper()))
    },
    ##
    getSrMapFun = function() {
      match.fun(paste0("run", self$getSrMapper()))
    },
    ##
    ## Predicate methods ####
    ##
    hasPileup = function() {
      tryCatch(
        is(self$mapInit$pileup, "pileup"),
        error = function(e)
          FALSE
      )
    },
    ##
    hasPartition = function() {
      tryCatch(
        is(self$partition$prt, "HapPart"),
        error = function(e)
          FALSE
      )
    },
    ##
    hasHapList = function() {
      tryCatch(
        is(self$partition$hpl, "HapList"),
        error = function(e)
          FALSE
      )
    },
    ## Get the absolut path
    absPath = function(filename) {
      sapply(filename, function(x) {
        normalizePath(
          file.path(self$getOutdir(), x),
          mustWork = FALSE)
      })
    },
    ## Get the relative path
    relPath = function(filepath) {
      gsub("^/", "", gsub(self$getOutdir(),"",filepath))
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
        LR = self$getPileup(),
        SR = self$mapInit$SR2$pileup
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
        pileup,
        threshold,
        range,
        thin,
        width,
        label %|ch|% lbl,
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
      }else {
        plotPartitionHistogramMulti(x = self$getPartition(), 
                                       limits = self$getLimits(), 
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
    plotPartitionHaplotypes = function(thin = 1, label = "") {
      tag <- self$getMapTag()
      plotPartitionHaplotypes(x = self$getPartition(), thin, label %|ch|% tag)
    },
    ##
    ## Summary methods ####
    ##
    polymorphicPositions = function(threshold, useSR = FALSE) {
      if (missing(threshold)) {
        threshold <- self$getThreshold()
      }
      if (useSR) {
        pileup <- self$mapInit$SR2$pileup
      } else {
        pileup <- self$getPileup()
      }
      polymorphicPositions(pileup, threshold)
    },
    ##
    plotmapInitSummary = function(thin = 0.2, width = 4, label = "") {
      readtypes <- if (self$getPartSR()) {
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
    plotmapIterSummary = function(thin = 0.2, width = 10, iteration = 0, 
                                  drop.indels = TRUE) {
      hptypes <- self$getHapTypes()
      if (iteration == self$getIterations()) {
        drop.indels <- FALSE
        width <- 10
      }
      plotlist <- foreach(hp = hptypes) %do% {
        tag <- self$getMapTag(iteration, hp)
        self$plotGroupCoverage(hp, ref = NULL, iteration = iteration,
                                    threshold = NULL, range = NULL, thin, width,
                                    label = tag, drop.indels = drop.indels)
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
        geom_logo(pwm, method = "bits", seq_type = "dna", 
                             stack_width = 0.9) +
        facet_wrap(~seq_group, ncol = 1, strip.position = "left") +
        scale_x_continuous(labels = ppos, breaks = seq_along(ppos)) +
        theme_logo() +
        theme(axis.text.x  = ggplot2::element_text(size = 10, angle = 60),
              axis.title.y = ggplot2::element_blank(),
              axis.text.y  = ggplot2::element_blank(),
              axis.ticks.y = ggplot2::element_blank(),
              strip.text.y = ggplot2::element_text(face = "bold", 
                                                   size = 42, angle = 180))
    },
    ##
    plotmapFinalSummary = function(readtype, thin = 0.2, width = 10, 
                                   iteration = "final") {
      hptypes <- self$getHapTypes()
      plotlist <- foreach(hp = hptypes) %do% {
        tag <- self$getMapTag(iteration, hp, readtype)
        ref <- paste0(readtype, hp)
        # tag <- self$mapFinal$tag[[ref]]
        self$plotGroupCoverage(group = hp, ref = readtype, 
                               iteration = iteration, threshold = NULL, 
                               range = NULL, thin = thin,
                               width = width, label = tag)
      }
      
      cowplot::plot_grid(plotlist = plotlist, labels = hptypes)
    },
    ##
    getmapIterConseqPileup = function(group) {
      conseq <- self$getConseqs(group, 2)
      pileup <- self$getGroupPileup(group, reads = "pacbio2")
      cutoff <- self$mapIter[[group]]$params$pruningCutoff
      cm <- consmat.pileup(pileup, freq = FALSE)
      cm <- cm[, which(.colSums(cm, NROW(cm), NCOL(cm)) > 0)]
      cm <- .pruneConsensusMatrix(cm, cutoff = cutoff, verbose = FALSE)
      if (!is.null(attr(cm, "pruningPosition"))) {
        cm <- cm[-attr(cm, "pruningPosition"), ]
      }
      dplyr::bind_cols(dplyr::data_frame(
        seq = strsplit1(toString(conseq), "")
      ),
      as.data.frame(cm))
    }
  ),
  ##
  ## Private fields ####
  ##
  private = list(
    conf         = NULL,
    reportStatus = NULL
  )
)


# Helpers -----------------------------------------------------------------

findReads <- function(datadir, sampleId, locus) {
  locus <- sub("^HLA-", "", toupper(locus))
  locus <- sub("^KIR-", "", toupper(locus))
  filePattern <- paste0(sampleId, ".+", "fast(q|a)(\\.gz)?$")
  readPath <- dir(datadir, pattern = filePattern, full.names = TRUE)
  readPath <- readPath[grep(pattern = "^((?!_trimmed.fastq).)*$", 
                              readPath, perl = TRUE)]
  if (length(readPath) > 0) {
    normalizePath(readPath, mustWork = TRUE)
  } else readPath
}
