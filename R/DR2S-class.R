#' @export
DR2Smap.default <- function(sample,
                            locus,
                            longreads = list(type = "pacbio", dir = "pacbio"),
                            shortreads = list(type = "illumina", dir = "illumina"),
                            datadir = ".",
                            outdir = ".",
                            reference = NULL,
                            consensus = "multialign",
                            threshold = 0.20,
                            iterations = 1,
                            microsatellite = FALSE,
                            partSR = TRUE,
                            fullname = TRUE,
                            forceBadMapping = FALSE,
                            filterScores = TRUE,
                            dist_alleles = 2,
                            create_outdir = TRUE,
                            ...) {
  conf <- create_dr2s_conf(
    sample = sample,
    locus = locus,
    longreads = longreads,
    shortreads = shortreads,
    datadir = datadir,
    outdir = outdir,
    reference = reference,
    consensus = consensus,
    threshold = threshold,
    iterations = iterations,
    microsatellite = microsatellite,
    filterScores = filterScores,
    partSR = partSR,
    fullname = fullname,
    dist_alleles = dist_alleles,
    forceBadMapping = forceBadMapping,
    ...
  )
  DR2S_$new(conf, create_outdir = create_outdir)
}

#' @export
DR2Smap.DR2Sconf <- function(sample) DR2S_$new(conf)

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
#' @usage DR2Smap(sample, locus, longreads = list(type = "pacbio", dir = "pacbio"),
#' shortreads = list(type = "illumina", dir = "illumina"), datadir = ".",
#' outdir = "./output", reference = NULL, consensus = "multialign",
#' threshold = 0.20, iterations = 1, microsatellite = FALSE, dist_alleles = 2, filterScores = TRUE, partSR = TRUE, fullname = TRUE, create_outdir = TRUE, ...)
#' @field mapInit \code{[mapInit]}; the mapping of long reads to \code{reference}.
#' @field partition \code{[PartList]}; the partitioning of full-length mapped
#'   long reads into different haplotypes.
#' @field mapIter \code{[mapIter]}; the mapping of A and B reads to consensus
#'   sequences produced in the previous step.
#' @field mapFinal \code{[mapFinal]}; the mapping of A, B, and short reads to consensus
#'   sequences produced in the previous step.
#' @field consensus \code{[ConsList]}; final consensus sequences for A and B.
#'
#' @keywords data internal
#' @return Object of \code{\link{R6Class}} representing a DR2S analysis.
#' @section Public Methods:
#' \describe{
#' \item{\code{x$runMapInit(opts = list(), optsname = "", pct = 100, threshold = 0.20, iterations = 1,
#' microsatellite = FALSE, dist_alleles = 2, filterScores = TRUE, partSR = TRUE, min_base_quality = 3, min_mapq = 0, max_depth = 1e4, min_nucleotide_depth = 3,
#' include_deletions = FALSE, include_insertions = FALSE, force = FALSE,
#' fullname = TRUE, plot = TRUE)}}{
#' Run the inital mapping step (long reads against the reference allele)}
#' \item{\code{x$runHaplotypePartitioning(max_depth = 1e4, shuffle = TRUE,
#' skip_gap_freq = 2/3, plot = TRUE)}}{
#' Partition mapped long reads into haplotypes}
#' \item{\code{x$splitReadsByHaplotype(limit = list(), #' plot = TRUE)}}{
#' Partition mapped long reads into haplotypes}
#' \item{\code{x$extractFastq(nreads = NULL, replace = FALSE, nalign = 30)}}{
#' Extract FASTQs for partitioned reads and if the consensus method is
#' \code{multialign} contruct a multiple sequence alignment from the best
#' \code{nalign} reads.}
#' }
#' @section Internals:
#' \describe{
#' \item{\code{x$new(conf, create_outdir = TRUE)}}{Initialise object of class
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
    initialize  = function(conf, create_outdir = TRUE) {
      ## Public fields
      self$mapInit   = list()
      self$partition = list()
      self$mapIter   = list()#A = list(), B = list())
      self$mapFinal  = list()
      self$consensus = list()
      ## Private fields
      private$conf   = initialise_dr2s(conf, create_outdir = create_outdir)
      # private$ipd    = findLocus(private$conf$locus)
      private$reportStatus = FALSE

      if (file.exists(ref_path <- private$conf$reference)) {
        private$conf$reference = basename(ref_path)
        private$conf$ref_path  = normalizePath(ref_path, mustWork = TRUE)
      } else {
        private$conf$reference = expand_allele(conf$reference, conf$locus)
        # private$conf$ref_path  = generate_reference_sequence(private$ipd, private$conf$reference, private$conf$outdir, fullname = FALSE)
        private$conf$ref_path  = generate_reference_sequence(private$conf$reference,
                                                             private$conf$locus,
                                                             private$conf$outdir,
                                                             fullname = FALSE)
      }
      if (is.null(private$conf$alternate)) {
        private$conf["alternate"] = list(NULL)
      }
      else {
        if (file.exists(alt_path <- private$conf$alternate)) {
          private$conf$alternate = basename(alt_path)
          private$conf$alt_path  = normalizePath(alt_path, mustWork = TRUE)
        } else {
          private$conf$alternate = expand_allele(private$conf$alternate, private$conf$locus)
          # private$conf$alt_path  = generate_reference_sequence(private$ipd, private$conf$alternate, private$conf$outdir, fullname = FALSE)
          private$conf$ref_path  = generate_reference_sequence(private$conf$alternate,
                                                             private$conf$locus,
                                                             private$conf$outdir,
                                                             fullname = FALSE)
        }
      }
      conf_log(outdir = private$conf$outdir, logName = "info")
      flog.info("Creating DR2S Object", name = "info")
    },
    cache = function(outname) {
      if (missing(outname)) {
        outname <- paste("DR2S", self$getLrdType(), self$getLrMapper(), "rds", sep = ".")
      }
      path <- file.path(self$getOutdir(), outname)
      message("\nCaching ", dQuote(path), "\n")
      saveRDS(self$clone(deep = TRUE),
              file = path,
              compress = TRUE)
      invisible(self)
    },
    clear = function() {
      unlink(self$getOutdir(), recursive = TRUE)
      dir_create_if_not_exists(file.path(self$getOutdir()))

      if (file.exists(ref_path <- private$conf$reference)) {
        private$conf$reference = basename(ref_path)
        private$conf$ref_path  = normalizePath(ref_path, mustWork = TRUE)
      } else {
        private$conf$reference = expand_allele(private$conf$reference, private$conf$locus)
        private$conf$ref_path  = generate_reference_sequence(self$getReference(),
                                                             self$getLocus(),
                                                             self$getOutdir(),
                                                             fullname = FALSE)
      }
      if (is.null(private$conf$alternate)) {
        private$conf["alternate"] = list(NULL)
      }
      else {
        if (file.exists(alt_path <- private$conf$alternate)) {
          private$conf$alternate = basename(alt_path)
          private$conf$alt_path  = normalizePath(alt_path, mustWork = TRUE)
        } else {
          private$conf$alternate = expand_allele(private$conf$alternate, private$conf$locus)
          private$conf$ref_path  = generate_reference_sequence(self$getAlternate(),
                                                               self$getLocus(),
                                                               self$getOutdir(),
                                                               fullname = FALSE)
        }
      }


        ## Todo: rm if works
        # private$conf$ref_path = generate_reference_sequence(self$getIPD(), self$getReference(), self$getOutdir(), fullname = FALSE)
        # private$conf$alt_path = generate_reference_sequence(self$getIPD(), self$getAlternate(), self$getOutdir(), fullname = FALSE)
      invisible(self)
    },
    cleanup = function() {
      unlink(self$mapInit$bamfile)
      unlink(paste0(self$mapInit$bamfile, ".bai"))
      foreach(hpt = self$getHapTypes()) %do% {
        unlink(file.path(self$getOutdir(), hpt), recursive = TRUE)
      }
      # unlink(file.path(self$getOutdir(), "A"), recursive = TRUE)
      # unlink(file.path(self$getOutdir(), "B"), recursive = TRUE)
      unlink(file.path(self$getOutdir(), "final"), recursive = TRUE)
      invisible(self)
    },
    print = function() {
      fmt0 <- "DR2S mapper for sample <%s> locus <%s>\n"
      cat(sprintf(fmt0, self$getSampleId(), self$getLocus()))
      fmt1 <- "Reference alleles: <%s%s>\nLongreads: <%s> Shortreads: <%s>\nMapper: <%s>\nDatadir: <%s>\nOutdir: <%s>\n"
      cat(sprintf(fmt1,
        self$getReference(),
        if (is.null(self$getAlternate())) "" else paste0(", ", self$getAlternate()),
        self$getLrdType(), self$getSrdType(), self$getLrMapper(), self$getDatadir(), self$getOutdir()
      ))
      invisible(self)
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
    getForceBadMapping = function(forceBadMapping){
      self$getConfig("forceBadMapping")
    },
    ##
    setForceBadMapping = function(forceBadMapping){
      stopifnot(is.logical(forceBadMapping))
      self$setConfig("forceBadMapping", forceBadMapping)
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
    getDistAlleles = function() {
      self$getConfig("dist_alleles")
    },
    ##
    setDistAlleles = function(dist_alleles) {
      stopifnot(is.numeric(dist_alleles))
      self$setConfig("dist_alleles", dist_alleles)
      invisible(self)
    },
    ##
    getLongreads = function() {
      dir <- self$getLrdDir()
      readpath <- findReads(dir, self$getSampleId(), self$getLocus())
      if (is.null(readpath) || length(readpath) == 0) {
        flog.error("No reads available for readtype <%s>", self$getLrdType(), name = "info")
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
        flog.error("No reads available for readtype <%s>", self$getSrdType(), name = "info")
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
        return(filtered)
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
    getAlternate = function() {
      self$getConfig("alternate")
    },
    ##
    getSampleId = function() {
      self$getConfig("sample_id")
    },
    ##
    getLocus = function() {
      self$getConfig("locus")
    },
    ##
    getRefPath = function() {
      self$getConfig("ref_path")
    },
    ##
    getRefSeq = function() {
      if (!is.null(self$getRefPath())) {
        Biostrings::readDNAStringSet(self$getRefPath())
      } else  {
        ipd.Hsapiens.db::getClosestComplete(ipd.Hsapiens.db::ipd.Hsapiens.db,
                                            self$getReference())
      }
      #   if (self$getReference() == "consensus") {
      #     self$getIPD()$cons
      #   } else {
      #
      #     self$getIPD()$get_reference_sequence(self$getReference(), "a")
      #   }
      # }
    },
    ##
    getAltPath = function() {
      self$getConfig("alt_path")
    },
    ##
    getAltSeq = function() {
      if (!is.null(self$getAltPath())) {
        Biostrings::readDNAStringSet(self$getRefPath())
      } else if (self$getAlternate() %in%
                 ipd.Hsapiens.db::getAlleles(ipd.Hsapiens.db::ipd.Hsapiens.db,
                                             self$getLocus()) ){
        ipd.Hsapiens.db::getClosestComplete(ipd.Hsapiens.db::ipd.Hsapiens.db,
                                            self$getAlternate())
      } else NULL
      # if (!is.null(self$getAlternate())) {
      #   if (!is.null(self$getIPD())) {
      #     self$getIPD()$get_reference_sequence(self$getAlternate())
      #   } else {
      #     Biostrings::readDNAStringSet(self$getAltPath())
      #   }
      # } else NULL
    },
    ##
    getNreads = function() {
      self$getConfig("nreads")
    },
    ##
    setNreads = function(nreads) {
      stopifnot(is.null(nreads) || nreads >= 1)
      self$setConfig("nreads", nreads)
      invisible(self)
    },
    ##
    getOpts = function(name = "mapInit") {
      if (is.null(name))
        merge_list(self$getConfig("opts"), self$getConfig("longreads")$opts, update = TRUE)
      else {
        opts <- merge_list(self$getConfig("opts"), self$getConfig("longreads")$opts, update = TRUE)
        if (is.list(opts))
          opts[[name]]
        else
          opts
      }
    },
    ##
    # getIPD = function() {
    #   private$ipd
    # },
    # ##
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
        return(self$mapFinal$seq)
      } else if (length(self$mapIter) > 0) {
        return(sapply(self$mapIter[max(names(self$mapIter))][self$getHapTypes()], function(x) x$seqpath))
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
        return(sapply(self$mapIter[max(names(self$mapIter))][self$getHapTypes()], function(x) x$conseq))
      }  else {
        return(self$getRefSeq())
      }
    },
    ##
    getConseqs = function(group = NULL, mapn = 1) {
      group <- match.arg(group, self$getHapTypes)
      if (mapn == 1) {
        ref <- self$mapIter[["0"]][[group]]$conseq$reference %||% Biostrings::BStringSet()
        alt <- self$mapIter[["0"]][[group]]$conseq$alternate %||% Biostrings::BStringSet()
        merged <- self$mapIter[["0"]][[group]]$conseq$merged %||% Biostrings::BStringSet()
      } else if (mapn == 2) {
        ref <- alt <- Biostrings::BStringSet()
        merged <- self$mapIter[[group]]$conseq %||% Biostrings::BStringSet()
      }
      seqs <- tryCatch({
        Biostrings::DNAStringSet(c(ref, alt, merged))
      }, error = function(e) {
        warning("Possibly indels in consensus")
        Biostrings::DNAStringSet(gsub("[-+]", "", c(ref, alt, merged)))
      })
      metadata(seqs) <- compact(list(
        ref = tryCatch(
          metadata(ref)$zscore,
          error = function(e)
            NULL
        ),
        alt =  tryCatch(
          metadata(alt)$zscore,
          error = function(e)
            NULL
        ),
        merged = tryCatch(
          metadata(merged)$zscore,
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
      match.fun(paste0("run_", self$getLrMapper()))
    },
    ##
    getSrMapFun = function() {
      match.fun(paste0("run_", self$getSrMapper()))
    },
    ##
    ## Predicate methods ####
    ##
    hasMultialign = function() {
      self$getConfig("consensus") == "multialign"
    },
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
      plot_pileup_coverage(
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
      plot_pileup_coverage(
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
    plotBasecallFrequency = function(threshold, label = "", drop.indels = FALSE) {
      tag <- self$getMapTag()
      if (missing(threshold)) {
        threshold <- self$getThreshold()
      }
      plot_pileup_basecall_frequency(
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
      plot_pileup_basecall_frequency(
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
         plot_partition_histogram(x = self$getPartition(), label %|ch|% tag, limits = limits)
      }else {
         plot_partition_histogram_multi(x = self$getPartition(), limits = self$getLimits(), label %|ch|% tag)
      }
    },
    ##
    plotPartitionTree = function() {
      plot_partition_tree(x = self$getPartition())
    },
    ##
    plotPartitionRadar = function() {
      plot_radar_partition(x = self$getPartition())
    },
    ##
    plotPartitionHaplotypes = function(thin = 1, label = "") {
      tag <- self$getMapTag()
      plot_partition_haplotypes(x = self$getPartition(), thin, label %|ch|% tag)
    },
    ##
    plotConseqProbability = function(ref = NULL,
                                     iteration = "init",
                                     threshold = "auto",
                                     text_size = 3,
                                     point_size = 1) {
      hptypes <- self$getHapTypes()
      cseqs <- foreach(hp = hptypes) %do% {
        seqs <- list(
          label = self$getMapTag(iteration, hp),
          if(is.numeric(iteration)) {
            cseq = self$mapIter[[as.character(iteration)]][[hp]]$conseq
          } else {
            cseq = self$mapFinal$conseq[[hp]]
          } )
        names(seqs) <- c("label", "cseq")
        seqs
      }

      # label = self$getMapTag(iteration, hp)

      names(cseqs) <- hptypes
      plot_conseq_probability(cseqs, threshold, text_size, point_size)

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
        polymorphic_positions(pileup, threshold)
    },
    ##
    summarisePartitions = function() {
      stopifnot(self$hasHapList())
      summary(self$getHapList())
    },
    ##
    plotmapInitSummary = function(thin = 0.2, width = 4, label = "", readtype = "LR") {
      # debug
      # readtype = "SR"
      tag <- self$getMapTag(ref = readtype)
      multiplot(
        self$plotCoverage(
          thin = thin,
          width = width,
          label = label %|ch|% tag,
          readtype = readtype
        ),
        self$plotBasecallFrequency(threshold = self$getThreshold(), label = " "),
        layout = matrix(c(1, 2, 2, 2), ncol = 1)
      )
    },
    ##
    plotPartitionSummary = function(label = "", limits = NULL) {
      tag <- self$getMapTag()
      p1 <-  self$plotPartitionHistogram(label = label %||% tag, limits = limits) +
                    ggplot2::theme(legend.position = "none")
      p2  <- self$plotPartitionTree()
      p3 <- self$plotPartitionRadar()
      if (!is.null(p2)) {
        multiplot(p1, p2, p3, layout = matrix(c(1, 2, 2, 2, 3), ncol = 1))
      } else {
        multiplot(p1,p3, layout = matrix(c(1,2)))
      }
    },
    ##
    plotmapIterSummary = function(thin = 0.2, width = 10, iteration = 0, drop.indels = TRUE) {
      hptypes <- self$getHapTypes()
      if (iteration > 1) {
        drop.indels <- FALSE
        width <- 10
      }
      plotlist <- foreach(hp = hptypes) %do% {
        tag <- self$getMapTag(iteration, hp)
        list(self$plotGroupCoverage(hp, ref = NULL, iteration = iteration,
                                    threshold = NULL, range = NULL, thin, width,
                                    label = tag, drop.indels = drop.indels),
        self$plotGroupBasecallFrequency(hp, NULL, iteration, NULL, " ", drop.indels = drop.indels))
      }
      plotlist <- unlist(plotlist, recursive = FALSE)
      layout_vector <- rep(1:length(plotlist), rep(c(1,3), length(hptypes)))
      layout = matrix(layout_vector, ncol = length(hptypes))

      multiplot(
        plotlist = plotlist,
        layout = layout)
    },
    ##
    plotmapFinalSummary = function(readtype, thin = 0.2, width = 10, iteration = "final") {
      hptypes <- self$getHapTypes()
      plotlist <- foreach(hp = hptypes) %do% {
        tag <- self$getMapTag(iteration, hp, readtype)
        ref <- paste0(readtype, hp)
        # tag <- self$mapFinal$tag[[ref]]
          self$plotGroupCoverage(group = hp, ref = readtype, iteration = iteration,
                                 threshold = NULL, range = NULL, thin = thin,
                                 width = width, label = tag)
        list(
          self$plotGroupCoverage(group = hp, ref = readtype, iteration = iteration,
                                 threshold = NULL, range = NULL, thin = thin,
                                 width = width, label = tag),
          self$plotGroupBasecallFrequency(group = hp, ref = readtype, iteration = iteration,
                                          threshold = NULL, label = " "))
      }
          self$plotGroupBasecallFrequency(group = hp, ref = readtype, iteration = iteration,
                                          threshold = NULL, label = " ")
      plotlist <- unlist(plotlist, recursive = FALSE)
      layout_vector <- rep(1:length(plotlist), rep(c(1,3), length(hptypes)))
      layout = matrix(layout_vector, ncol = length(hptypes))
      multiplot(
        plotlist = plotlist,
        layout = layout)
    },
    ##
    plotmapIterSummaryConseqProb = function(iteration = 0,
                                         text_size = 1.75,
                                         point_size = 0.75,
                                         threshold = "auto") {
      # debug
      # reads = "pacbio2"
      # ref = NULL
      print(
        self$plotConseqProbability(
          iteration = iteration,
          text_size = text_size,
          point_size = point_size,
          threshold = threshold
        )
      )
    },
    ##
    plotmapFinalSummaryConseqProb = function(iteration = "final",
                                         text_size = 1.75,
                                         point_size = 0.75,
                                         threshold = "auto") {
      print(
        self$plotConseqProbability(
          iteration = iteration,
          text_size = text_size,
          point_size = point_size,
          threshold = threshold
        )
      )
    },
    ##
    getmapIterConseqPileup = function(group) {
      conseq <- self$getConseqs(group, 2)
      pileup <- self$getGroupPileup(group, reads = "pacbio2")
      cutoff <- self$mapIter[[group]]$params$pruning_cutoff
      cm <- consmat.pileup(pileup, freq = FALSE)
      cm <- cm[, which(.colSums(cm, NROW(cm), NCOL(cm)) > 0)]
      cm <- prune_consensus_matrix(cm, cutoff = cutoff, verbose = FALSE)
      if (!is.null(attr(cm, "pruning_position"))) {
        cm <- cm[-attr(cm, "pruning_position"), ]
      }
      dplyr::bind_cols(dplyr::data_frame(
        seq = strsplit(toString(conseq), "")[[1]],
        z   = S4Vectors::metadata(conseq)$merged,
        p   = pnorm(S4Vectors::metadata(conseq)$merged)
      ),
      as.data.frame(cm))
    }
  ),
  ##
  ## Private fields ####
  ##
  private = list(
    conf      = NULL,
    # ipd       = NULL,
    reportStatus = NULL,
    run_      = function(step) {
      switch(
        step,
        clear          = self$clear(),
        mapInit        = self$runMapInit(),
        partition      = self$runHaplotypePartitioning(),
        partShortReads = self$runPartitionShortReads(),
        split          = self$splitReadsByHaplotype(),
        extract        = self$extractFastq(),
        mapIter        = self$runMapIter(),
        mapFinal       = self$runMapFinal(),
        cache          = self$cache(),
        polish         = DR2S::polish(self),
        report         = DR2S::report(self),
        stop("<", step, "> is not a valid step in the mapping pipeline")
      )
    }
  )
)


# Helpers -----------------------------------------------------------------

findReads <- function(datadir, sample_id, locus) {
  locus <- sub("^HLA-", "", toupper(locus))
  locus <- sub("^KIR-", "", toupper(locus))
  file_pattern <- paste0(sample_id, ".+", "fast(q|a)(\\.gz)?$")
  read_path <- dir(datadir, pattern = file_pattern, full.names = TRUE)
  read_path <- read_path[grep(pattern = "^((?!_trimmed.fastq).)*$", read_path, perl = TRUE)]
  if (length(read_path) > 0) {
    normalizePath(read_path, mustWork = TRUE)
  } else read_path
}


# findLocus <- function(locus) {
#   stopifnot(requireNamespace("ipd.Hsapiens.db", quietly = TRUE))
#   IPDdata::loadIPDdata(locus)
# }
