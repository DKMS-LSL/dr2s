createDR2SConf <- function(sample,
                             locus,
                             longreads       = list(type = "pacbio", 
                                                    dir = "pacbio"),
                             shortreads      = list(type = "illumina", 
                                                    dir = "illumina"),
                             datadir         = ".",
                             outdir          = "./output",
                             reference       = NULL,
                             consensus       = NULL,
                             threshold       = 0.20,
                             iterations      = 1,
                             microsatellite  = FALSE,
                             distAlleles    = 2,
                             filterScores    = TRUE,
                             partSR          = TRUE,
                             forceMapping = FALSE,
                             fullname        = TRUE,
                             note            = NULL,
                             ...) {
  conf0 <- list(...)
  if (length(reference) == 0) {
    ref <- NULL
    alt <- NULL
  } else   if (length(reference) == 1) {
    ref <- reference
    alt <- NULL
  } else if (length(reference) == 2) {
    ref <- reference[1]
    alt <- reference[2]
  }
  conf1 <- list(
    datadir    = normalizePath(datadir, mustWork = TRUE),
    outdir     = outdir,
    threshold  = threshold,
    iterations = iterations,
    microsatellite  = microsatellite,
    distAlleles    = distAlleles,
    filterScores    = filterScores,
    partSR          = partSR,
    forceMapping = forceMapping,
    lrmapper   = conf0$lrmapper %||% "minimap",
    srmapper   = conf0$srmapper %||% "bwamem",
    limits     = conf0$limits   %||% NULL,
    haptypes   = conf0$haptypes %||% NULL,
    pipeline   = conf0$pipeline %||% c("clear", "mapInit", "partitionLongReads",
                                       "mapIter", "partitionShortReads",
                                       "mapFinal", "polish", "report"),
    longreads  = longreads,
    shortreads = shortreads,
    opts       = conf0$opts     %||% NULL,
    sampleId  = sample,
    locus      = locus,
    reference  = ref,
    alternate  = alt,
    consensus  = consensus      %||% "mapping",
    note       = note           %||% NULL
  )
  structure(conf1, class = c("DR2Sconf", "list"))
}

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
  conf$distAlleles   <- conf$distAlleles   %||% 2
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
  conf$note            <- conf$note          %||% NULL
  
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
  ## we can have different consensus calling methods
  cnss <- conf$consensus %||% list("mapping")
  conf$consensus <- NULL

  foreach(sample = samples, sampleId = sampleIds, .combine = "c") %:%
    foreach(nrd = nrds, .combine = "c") %:%
    foreach(lrd = lrds, .combine = "c") %:%
    foreach(cns = cnss, .combine = "c") %:%
    foreach(dst = sample$distDlleles, .combine = "c") %:%
    foreach(ref = sample$reference, .combine = "c") %do% {
      updateDR2SConf(conf, lrd, nrd, sampleId, sample, ref, alternate = "", 
                     cns, dst)
    }
}

updateDR2SConf <- function(conf0, lrd, nrd, sampleId, locus, reference, 
                           alternate, consensus, dst) {
  conf0$datadir   <- normalizePath(conf0$datadir, mustWork = TRUE)
  conf0$outdir    <- normalizePath(conf0$outdir, mustWork = FALSE)
  conf0$longreads <- lrd
  conf0$distAlleles  <- dst
  conf0$sampleId <- sampleId
  conf0["reference"] <- reference %|ch|% list(NULL)
  locus$reference    <- NULL
  conf0["alternate"] <- alternate %|ch|% list(NULL)
  locus$alternate    <- NULL
  conf0["consensus"] <- consensus %|ch|% list(NULL)
  locus$consensus    <- NULL
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
      # conf$reftype, ".", conf$consensus)
    )))
  }
  conf
}

validateDR2SConf <- function(conf) {
  fields <- c("datadir", "outdir", "threshold", "iterations", "microsatellite", 
              "distAlleles", "filterScores", "partSR", "forceMapping", 
              "lrmapper", "srmapper", "limits", "haptypes", "pipeline", 
              "longreads", "shortreads", "opts", "sampleId", "locus", 
              "reference", "alternate", "consensus", "note")
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
  ## Check consensus method
  conf$consensus <- match.arg(conf$consensus, c("mapping"))
  ## Reference type
  conf$reftype <- if (is.null(conf$reference) && is.null(conf$alternate)) {
    stop("Need to provide a reference, either as fasta or as ipd identifier")
  } else if (!is.null(conf$reference) && is.null(conf$alternate)) {
    "ref"
  } else if (!is.null(conf$reference) && !is.null(conf$alternate)) {
    "ref.alt"
  }
  conf
}
#' @export
print.DR2Sconf <- function(x, ...) {
  cat(yaml::as.yaml(x$getConfig, indent = 4))
  invisible(x)
}

writeDR2SConf <- function(x, outFile = NULL) {
  if (is.null(outFile)) {
    outFile <- file.path(x$getOutdir(), "config.yaml")
  }
  yaml::write_yaml(yaml::as.yaml(x$getConfig(), indent = 4), outFile)
}
