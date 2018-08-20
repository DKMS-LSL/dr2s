#' @export
createDR2SConf <- function(sample,
                           locus,
                           longreads       = list(type = "pacbio",
                                                  dir = "pacbio"),
                           shortreads      = list(type = "illumina",
                                                  dir = "illumina"),
                           datadir         = ".",
                           outdir          = "./output",
                           reference       = NULL,
                           threshold       = 0.20,
                           iterations      = 1,
                           microsatellite  = FALSE,
                           distAlleles     = 2,
                           filterScores    = TRUE,
                           partSR          = TRUE,
                           forceMapping    = FALSE,
                           fullname        = TRUE,
                           details         = NULL,
                           ...) {
  conf0 <- list(...)
  conf1 <- list(
    datadir        = normalizePath(datadir, mustWork = TRUE),
    outdir         = outdir,
    threshold      = threshold,
    iterations     = iterations,
    microsatellite = microsatellite,
    distAlleles    = distAlleles,
    filterScores   = filterScores,
    partSR         = partSR,
    forceMapping   = forceMapping,
    lrmapper       = conf0$lrmapper %||% "minimap",
    srmapper       = conf0$srmapper %||% "bwamem",
    limits         = conf0$limits   %||% NULL,
    haptypes       = conf0$haptypes %||% NULL,
    pipeline       = conf0$pipeline %||% c("clear", "mapInit",
                                           "partitionLongReads", "mapIter",
                                           "partitionShortReads", "mapFinal",
                                           "polish", "report"),
    longreads      = longreads,
    shortreads     = shortreads,
    opts           = conf0$opts     %||% NULL,
    sampleId       = sample,
    locus          = locus,
    reference      = reference,
    details        = gsub(";", ",", details)        %||% NULL
  )
  structure(conf1, class = c("DR2Sconf", "list"))
}

#' Read a DR2S config file in yaml format
#' @param configFile The path to the config file.
#' @details DR2S config files can be created manually or by the
#' \code{\link{writeDR2SConf}} function.
#' @export
readDR2SConf <- function(configFile) {
  conf <- yaml::yaml.load_file(configFile)
  ## set defaults if necessary
  conf$datadir        <- conf$datadir        %||% normalizePath(".",
                                                                mustWork = TRUE)
  conf$outdir         <- conf$outdir         %||% file.path(conf$datadir,
                                                            "output")
  conf$threshold      <- conf$threshold      %||% 0.2
  conf$iterations     <- conf$iterations     %||% 2
  conf$microsatellite <- conf$microsatellite %||% FALSE
  conf$filterScores   <- conf$filterScores   %||% TRUE
  conf$partSR         <- conf$partSR         %||% TRUE
  conf$forceMapping   <- conf$forceMapping   %||% FALSE
  conf$lrmapper       <- conf$lrmapper       %||% "minimap"
  conf$srmapper       <- conf$srmapper       %||% "bwamem"
  conf$limits         <- conf$limits         %||% list(NULL)
  conf$haptypes       <- conf$haptypes       %||% list(NULL)
  conf$distAlleles    <- conf$distAlleles    %||% 2
  conf$pipeline       <- conf$pipeline       %||% c("clear", "mapInit",
                                                    "partitionLongReads",
                                                    "mapIter",
                                                    "partitionShortReads",
                                                    "mapFinal",
                                                    "polish", "report")
  conf$longreads       <- conf$longreads     %||% list(type = "pacbio",
                                                       dir = "pacbio")
  conf$shortreads      <- conf$shortreads    %||% list(type = "illumina",
                                                       dir = "illumina")
  conf$details         <- gsub(";", ",", conf$details)       %||% list(NULL)

  if (length(conf$shortreads) == 1 && is.list(conf$shortreads[[1]]))
    conf$shortreads <- conf$shortreads[[1]]

  if (is.null(conf$opts))
    conf["opts"] <-  list(NULL)
  expandDR2SConf(structure(conf, class = c("DR2Sconf", "list")))
}

expandDR2SConf <- function(conf) {
  ## we can have more than one sample
  samples <- conf$samples
  conf$samples <- NULL
  sampleIds <- names(samples)
  ## we can have more than one longread type
  lrds <- conf$longreads
  if (!is.null(names(lrds))) {
    lrds <- list(lrds)
  }
  conf$longreads <- NULL

  foreach(sample = samples, sampleId = sampleIds, .combine = "c") %:%
    foreach(lrd = lrds, .combine = "c") %:%
    foreach(dst = sample$distAlleles, .combine = "c") %:%
    foreach(ref = sample$reference, .combine = "c") %do% {
      updateDR2SConf(conf, lrd, sampleId, sample, ref,
                     dst)
    }
}

updateDR2SConf <- function(conf0, lrd, sampleId, locus, reference,
                          dst) {
  conf0$datadir   <- normalizePath(conf0$datadir, mustWork = TRUE)
  conf0$outdir    <- normalizePath(conf0$outdir, mustWork = FALSE)
  conf0$longreads <- lrd
  conf0$distAlleles  <- dst
  conf0$sampleId <- sampleId
  conf0["reference"] <- reference %|ch|% list(NULL)
  locus$reference    <- NULL
  ## add overides if they exist
  list(.mergeList(conf0, locus, update = TRUE))
}

initialiseDR2S <- function(conf, createOutdir = TRUE) {
  conf <- validateDR2SConf(conf)
  if (createOutdir) {
    conf$outdir <- .dirCreateIfNotExists(path = gsub("//+", "/", file.path(
      conf$outdir,
      conf$sampleId#,
      ## Use only sample id
      # conf$longreads$name %||% conf$longreads$dir %||% "",
      # paste0(.normaliseLocus(conf$locus), ".", conf$longreads$type, ".",
    )))
  }
  conf
}

validateDR2SConf <- function(conf) {
  fields <- c("datadir", "outdir", "threshold", "iterations", "microsatellite",
              "distAlleles", "filterScores", "partSR", "forceMapping",
              "lrmapper", "srmapper", "limits", "haptypes", "pipeline",
              "longreads", "shortreads", "opts", "sampleId", "locus",
              "reference", "details")
  if (!all(fields %in% names(conf))) {
    stop("Missing fields <", comma(fields[!fields %in% names(conf)]),
         "> in config", call. = FALSE)
  }
  conf <- structure(conf[fields], class = c("DR2Sconf", "list"))
  conf$datadir <- normalizePath(conf$datadir, mustWork = TRUE)
  conf$outdir  <- normalizePath(conf$outdir, mustWork = FALSE)
  if (!is.numeric(conf$threshold) ||
      (conf$threshold <= 0 ||
       conf$threshold >= 1)) {
    stop("The polymorphism threshold must be a number between 0 ans 1",
         call. = FALSE)
  }
  # check iterations
  if (!is.numeric(conf$iterations) ||
      !conf$iterations&&1 ==  0 ||
      conf$iterations > 10) {
    stop("The number of iterations must be integers between 1 and 10",
         call. = FALSE)
  }
  # check partSR
  if (!is.logical(conf$partSR)){
    stop("partSR must be logical", call. = FALSE)
  }
  # check filterScores
  if (!is.logical(conf$filterScores)){
    stop("filterScores must be logical", call. = FALSE)
  }
  # check microsatellite
  if (!is.logical(conf$microsatellite)){
    stop("microsatellite must be logical", call. = FALSE)
  }
  # check number of distinct alleles
  if (!is.numeric(conf$distAlleles)){
    stop("Number of distinct alleles (distAlleles) must be numeric",
         call. = FALSE)
  }
  # check forceMapping
  if (!is.logical(conf$forceMapping)){
    stop("forceMapping must be logical", call. = FALSE)
  }
  ## Check mapper
  conf$lrmapper <- match.arg(conf$lrmapper, c("bwamem", "graphmap", "minimap"))
  conf$srmapper <- match.arg(conf$srmapper, c("bwamem", "graphmap", "minimap"))
  ## Check pipeline
  pipesteps <- c("clear", "cache", "mapInit", "mapIter", "mapFinal",
                 "partitionLongReads", "partitionShortReads", "polish",
                 "report")
  if (!all(conf$pipeline %in% pipesteps)) {
    stop("Invalid pipeline step <",
         comma(conf$pipeline[!conf$pipeline %in% pipesteps]), ">",
         call. = FALSE)
  }
  ## Long reads must have a "type" and a "dir" field.
  lrdfields <- c("type", "dir")
  if (!all(lrdfields %in% names(conf$longreads))) {
    stop("Missing fields <",
         comma(lrdfields[!lrdfields %in% names(conf$longreads)]),
         "> in longreads config", call. = FALSE)
  }
  ## Type must be one of "pacbio" or "nanopore"
  if (!conf$longreads$type %in% c("pacbio", "nanopore")) {
    stop("Unsupported longread type <", conf$longreads$type, ">", call. = FALSE)
  }
  longreaddir <- file.path(conf$datadir, conf$longreads$dir)
  if (!file.exists(longreaddir)) {
    stop("No long read directory <", longreaddir, ">", call. = FALSE)
  }
  if (length(dir(longreaddir, pattern = ".+fast(q|a)(.gz)?")) == 0) {
    stop("No FASTQ files in long read directory <", longreaddir, ">",
         call. = FALSE)
  }
  ## Check short reads (short reads are optional i.e. the field can be NULL)
  if (!is.null(conf$shortreads)) {
    srdfields <- c("type", "dir")
    if (!all(srdfields %in% names(conf$shortreads))) {
      stop("Missing fields <",
           comma(srdfields[!srdfields %in% names(conf$shortreads)]),
           "> in shortreads config", call. = FALSE)
    }
    if (!conf$shortreads$type %in% c("illumina")) {
      stop("Unsupported shortread type <", conf$shortreads$type, ">",
           call. = FALSE)
    }
    shortreaddir <- file.path(conf$datadir, conf$shortreads$dir)
    if (!file.exists(longreaddir)) {
      warning("No short read directory <", shortreaddir, ">", call. = FALSE,
              immediate. = TRUE)
    }
    if (length(dir(shortreaddir, pattern = ".+fastq(.gz)?")) == 0) {
      warning("No FASTQ files in short read directory <",
              shortreaddir, ">", call. = FALSE, immediate. = TRUE)
    }
  }
  ## Normalise locus
  ##conf$locus   <- sub("^(HLA-|KIR)", "", .normaliseLocus(conf$locus))
  conf$locus   <- sub("^HLA-", "", .normaliseLocus(conf$locus))
  ## Reference type
  conf$reftype <- if (is.null(conf$reference)) {
    stop("Need to provide a reference, either as fasta or as ipd identifier")
  } else {
    "ref"
  }
  conf
}
#' @export
print.DR2Sconf <- function(x, ...) {
  cat(yaml::as.yaml(x, indent = 4))
  invisible(x)
}

#' Write the DR2S config as yaml.
#' @param x A \code{\link{DR2S}} object.
#' @param outFile The destination file.
#' @export
writeDR2SConf <- function(x, outFile = NULL) {
  if (is.null(outFile)) {
    outFile <- file.path(x$getOutdir(), "config.yaml")
  }
  sample <- list(list(reference = x$getReference(),
                      locus = x$getLocus(),
                      distAlleles = x$getDistAlleles(),
                      details = x$getDetails()))
  names(sample) <- x$getSampleId()
  conf <- list(
    datadir         = x$getDatadir(),
    outdir          = dirname(x$getOutdir()),
    threshold       = x$getThreshold(),
    forceMapping    = x$getForceMapping(),
    iterations      = x$getIterations(),
    microsatellites = x$getMicrosatellite(),
    srmapper        = x$getSrMapper(),
    lrmapper        = x$getLrMapper(),
    partSR          = x$getPartSR(),
    filterScores    = x$getFilterScores(),
    pipeline        = x$getPipeline(),
    opts            = x$getConfig("opts"),
    longreads       = x$getConfig("longreads"),
    shortreads      = x$getConfig("shortreads"),
    samples         = sample
  )
  yaml::write_yaml(conf, outFile)
}
