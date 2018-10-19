#' @export
createDR2SConf <- function(sample,
                           locus,
                           longreads = list(type = "pacbio", dir = "pacbio"),
                           shortreads = NULL,
                           datadir = ".",
                           outdir = "./output",
                           reference = NULL,
                           threshold = 0.20,
                           iterations = 1,
                           microsatellite = FALSE,
                           distAlleles = 2,
                           filterScores = FALSE,
                           forceMapping = FALSE,
                           details = NULL,
                           ...) {
  conf0 <- list(...)
  conf1 <- list(
    datadir = normalizePath(datadir, mustWork = TRUE),
    outdir = outdir,
    threshold = threshold,
    iterations = iterations,
    microsatellite = microsatellite,
    distAlleles = distAlleles,
    filterScores = filterScores,
    forceMapping = forceMapping,
    lrmapper = conf0$lrmapper %||% "minimap",
    srmapper = conf0$srmapper %||% "bwamem",
    limits = conf0$limits   %||% NULL,
    haptypes = conf0$haptypes %||% NULL,
    pipeline = conf0$pipeline %||% if (is.null(shortreads)) {
      LR_PIPELINE()
    } else {
      SR_PIPELINE()
    },
    longreads = longreads,
    shortreads = shortreads,
    opts = conf0$opts %||% NULL,
    sampleId = sample,
    locus = locus,
    reference = reference,
    details = lapply(details, gsub, pattern = ";", replacement = ",") %||% NULL
  )
  structure(conf1, class = c("DR2Sconf", "list"))
}

#' Read a DR2S config file in yaml format
#' @param configFile The path to the valid DR2S config file.
#' @details DR2S config files can be created manually or by the
#' \code{\link{writeDR2SConf}} function.
#' @return A \code{DR2Sconf} object or a list of \code{DR2Sconf} objects.
#' @export
readDR2SConf <- function(configFile) {
  conf <- yaml::yaml.load_file(configFile)
  ## set defaults if necessary
  conf["datadir"] <- conf$datadir %||% normalizePath(".", mustWork = TRUE)
  conf["outdir"] <- conf$outdir %||% file.path(conf$datadir, "output")
  conf["threshold"] <- conf$threshold %||% 0.2
  conf["iterations"] <- conf$iterations %||% 2
  conf["microsatellite"] <- conf$microsatellite %||% FALSE
  conf["filterScores"] <- conf$filterScores %||% TRUE
  conf["forceMapping"] <- conf$forceMapping %||% FALSE
  conf["lrmapper"] <- conf$lrmapper %||% "minimap"
  conf["srmapper"] <- conf$srmapper %||% "bwamem"
  conf$pipeline <- conf$pipeline %||% if (is.null(conf$shortreads)) {
    LR_PIPELINE()
  } else {
    SR_PIPELINE()
  }
  conf$longreads <- conf$longreads %||% list(type = "pacbio", dir = "pacbio")
  if (is.null(conf$shortreads))
    conf["shortreads"] <-  list(NULL)
  if (is.null(conf$opts))
    conf["opts"] <- list(NULL)
  if (is.null(conf$limits))
    conf["limits"] <- list(NULL)
  if (is.null(conf$haptypes))
    conf["haptypes"] <- list(NULL)

  conf <- expandDR2SConf(conf)
  if (length(conf) == 1) {
    conf[[1]]
  } else {
    conf
  }
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
  # lrd <- lrds[[1]]
  # sample <- samples[[1]]
  # sampleId <- sampleIds[[1]]
  foreach(sample = samples, sampleId = sampleIds, .combine = "c") %:%
    foreach(lrd = lrds, .combine = "c") %do% {
      updateDR2SConf(conf, lrd, sampleId, sample)
    }
}

updateDR2SConf <- function(conf0, lrd, sampleId, sample) {
  conf0$datadir <- normalizePath(conf0$datadir, mustWork = TRUE)
  conf0$outdir <- normalizePath(conf0$outdir, mustWork = FALSE)
  conf0$longreads <- lrd
  conf0$sampleId <- sampleId
  conf0$distAlleles <- sample$distAlleles
  sample$distAlleles <- NULL
  conf0["reference"] <- sample$reference %|ch|% list(NULL)
  sample$reference <- NULL
  ## add overides if they exist
  conf1 <- mergeList(conf0, sample, update = TRUE)
  list(validateDR2SConf(conf1))
}

initialiseDR2S <- function(conf, createOutdir = TRUE) {
  conf <- validateDR2SConf(conf)
  if (createOutdir) {
    conf$outdir <- .dirCreateIfNotExists(path = gsub("//+", "/", file.path(
      conf$outdir,
      conf$sampleId
    )))
  }
  conf
}

validateDR2SConf <- function(conf) {
  fields <- c("datadir", "outdir", "threshold", "iterations", "microsatellite",
              "filterScores", "forceMapping", "lrmapper", "srmapper",
              "limits", "haptypes", "pipeline", "longreads", "shortreads",
              "opts", "sampleId", "locus", "reference", "distAlleles",
              "details")
  assert_that(all(fields %in% names(conf)),
              msg = paste("Missing fields <",
                          comma(fields[!fields %in% names(conf)]),
                          "> in config"))
  conf <- structure(conf[fields], class = c("DR2Sconf", "list"))
  conf$datadir <- normalizePath(conf$datadir, mustWork = TRUE)
  conf$outdir  <- normalizePath(conf$outdir, mustWork = FALSE)
  assert_that(
    is.numeric(conf$threshold),
    conf$threshold >= 0,
    conf$threshold <= 1,
    msg = "The polymorphism threshold must be a number between 0 ans 1")

  # check iterations
  assert_that(
    is.count(conf$iterations),
    conf$iterations > 0,
    conf$iterations < 10,
    msg = "The number of iterations must be integers between 1 and 10"
  )
  # check all logicals
  assert_that(
    is.logical(conf$filterScores),
    is.logical(conf$microsatellite),
    is.logical(conf$forceMapping)
  )
  # check number of distinct alleles
  assert_that(is.count(conf$distAlleles),
              msg = "Number of distinct alleles (distAlleles) must be numeric")

  ## Check mapper
  conf$lrmapper <- match.arg(conf$lrmapper, c("bwamem", "minimap"))
  conf$srmapper <- match.arg(conf$srmapper, c("bwamem", "minimap"))

  ## Check pipeline
  pipesteps <- c("clear", "cache", "mapInit", "mapIter", "mapFinal",
                 "partitionLongReads", "partitionShortReads", "polish",
                 "report")
  assert_that(all(conf$pipeline %in% pipesteps),
              msg = paste(
                "Invalid pipeline step <",
                comma(conf$pipeline[!conf$pipeline %in% pipesteps]), ">") )

  ## Check long reads
  ## Long reads must have a "type" and a "dir" field.
  readfields <- c("type", "dir")
  assert_that(
    all(readfields %in% names(conf$longreads)),
    msg = paste("Missing fields <",
                comma(readfields[!readfields %in% names(conf$longreads)]),
                "> in longreads config"))
  ## Type must be one of "pacbio" or "nanopore"
  assert_that(
    conf$longreads$type %in% c("pacbio", "nanopore"),
    msg = paste("Unsupported longread type", conf$longreads$type)
  )
  longreaddir <- file.path(conf$datadir, conf$longreads$dir)
  assert_that(is.dir(longreaddir))
  assert_that(length(dir(longreaddir, pattern = ".+fast(q|a)(.gz)?")) > 0,
              msg = "No fasta or fastq file in the longread directory"
  )

  ## Check short reads (short reads are optional i.e. the field can be NULL)
  if (!is.null(conf$shortreads)) {
    assert_that(
      all(readfields %in% names(conf$shortreads)),
      msg = paste("Missing fields <",
                  comma(readfields[!readfields %in% names(conf$shortreads)]),
                  "> in shortreads config"))
    ## Type must be "illumina"
    assert_that(
      conf$shortreads$type == "illumina",
      msg = paste("Unsupported shortread type", conf$shortreads$type)
    )
    shortreaddir <- file.path(conf$datadir, conf$shortreads$dir)
    assert_that(is.dir(shortreaddir))
    assert_that(length(dir(shortreaddir, pattern = ".+fast(q|a)(.gz)?")) > 0,
                msg = "No fasta or fastq file in the shortread directory"
    )
  }
  ## Normalise locus
  conf$locus   <- sub("^HLA-", "", .normaliseLocus(conf$locus))
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
  sample <- list(
    list(
      locus = x$getLocus(),
      reference = x$getReference(),
      distAlleles = x$getDistAlleles(),
      details = x$getDetails()))
  names(sample) <- x$getSampleId()
  conf <- list(
    datadir        = x$getDatadir(),
    outdir         = dirname(x$getOutdir()),
    threshold      = x$getThreshold(),
    iterations     = x$getIterations(),
    microsatellite = x$getMicrosatellite(),
    filterScores   = x$getFilterScores(),
    forceMapping   = x$getForceMapping(),
    lrmapper       = x$getLrMapper(),
    srmapper       = x$getSrMapper(),
    limits         = x$getLimits(),
    haptypes       = x$getHapTypes(),
    pipeline       = x$getPipeline(),
    longreads      = x$getConfig("longreads"),
    shortreads     = x$getConfig("shortreads"),
    opts           = x$getConfig("opts"),
    samples        = sample
  )
  yaml::write_yaml(conf, outFile)
}

