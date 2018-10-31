#' @export
createDR2SConf <- function(sample,
                           locus,
                           longreads = list(dir = "pacbio", type = "pacbio", mapper = "minimap"),
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
                           details = list(),
                           ...) {
  conf0 <- list(...)
  datadir0 <- normalizePath(datadir, mustWork = TRUE)
  conf1 <- list(
    datadir = datadir0,
    outdir = outdir,
    threshold = threshold,
    iterations = iterations,
    microsatellite = microsatellite,
    distAlleles = distAlleles,
    filterScores = filterScores,
    forceMapping = forceMapping,
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
  conf["threshold"] <- conf$threshold %||% 0.2
  conf["iterations"] <- conf$iterations %||% 2
  conf["microsatellite"] <- conf$microsatellite %||% FALSE
  conf["filterScores"] <- conf$filterScores %||% TRUE
  conf["forceMapping"] <- conf$forceMapping %||% FALSE
  conf["distAlleles"] <- conf$distAlleles %||% 2
  conf["datadir"] <- conf$datadir %||% normalizePath(".", mustWork = TRUE)
  conf["outdir"] <- conf$outdir %||% file.path(conf$datadir, "output")
  conf$longreads <- conf$longreads %||% list(dir = "pacbio", type = "pacbio", mapper = "minimap")
  conf$pipeline  <- conf$pipeline  %||% if (is.null(conf$shortreads)) {
    LR_PIPELINE()
  } else {
    SR_PIPELINE()
  }
  if (is.null(conf$shortreads))
    conf["shortreads"] <- list(NULL)
  if (is.null(conf$opts))
    conf["opts"] <- list(NULL)
  if (is.null(conf$details))
    conf["details"] <- list(NULL)
  conf$runstats <- NULL
  conf <- expandDR2SConf(conf)
  if (length(conf) == 1) {
    conf[[1]]
  } else {
    conf
  }
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
  assert_that(is(x, "DR2S"))
  if (is.null(outFile)) {
    outFile <- file.path(x$getOutdir(), "config.yaml")
  }
  conf <- x$getConfig()
  conf$runstats <- x$getStats()
  conf$runstats$reported <- x$getReportStatus()
  yaml::write_yaml(conf, outFile)
}

expandDR2SConf <- function(conf) {
  ## we can have more than one sample
  samples <- conf$samples
  if (is.null(samples)) {
    return(list(validateDR2SConf(conf)))
  }
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


# Validator ---------------------------------------------------------------


MANDATORY_CONF_FIELDS <- function() {
  c(
    # Sample
    "sampleId", "locus",
    # Run parameters
    "threshold", "iterations", "microsatellite", "filterScores", "forceMapping", "distAlleles",
    # Run data
    "datadir", "outdir", "reference", "longreads",
    # Options
    "pipeline"
  )
}

ORDERED_CONF_FIELDS <- function() {
  c(
    # Sample
    "sampleId", "locus", "details",
    # Run parameters
    "threshold", "iterations", "microsatellite", "filterScores", "forceMapping", "distAlleles",
    # Run data
    "datadir", "outdir", "reference", "longreads", "shortreads",
    # Options
    "pipeline", "opts"
  )
}

normaliseLongreads <- function(lrd) {
  ## Set defaults
  lrd$platform <- tolower(lrd$platform) %||% tolower(lrd$type)
  lrd$mapper   <- tolower(lrd$mapper) %||% "minimap"

  ## Long reads must have a "dir" or "file", and a "platform" and "mapper" field.
  fields0 <- c("dir", "file", "platform", "mapper")
  lrd <- compact(lrd[fields0])
  fields1 <- c("file", "platform", "mapper")
  fields2 <- c("dir", "platform", "mapper")
  assert_that(
    all(fields1 %in% names(lrd)) || all(fields2 %in% names(lrd)),
    msg = paste0("Missing fields <",
                 comma(fields0[!fields0 %in% names(lrd)]),
                 "> in longreads config"))

  ## Platform must be one of "pacbio" or "nanopore"
  assert_that(
    lrd$platform %in% c("pacbio", "nanopore"),
    msg = paste0("Unsupported longread type <", lrd$platform, ">")
  )

  ## Mapper must be one of "bwamem" or "minimap"
  assert_that(
    lrd$mapper %in% c("bwamem", "minimap"),
    msg = paste0("Unsupported longread mapper <", lrd$mapper, ">")
  )

  lrd
}

normaliseShortreads <- function(srd) {
  ## Set defaults
  srd$platform <- tolower(srd$platform) %||% tolower(srd$type)
  srd$mapper   <- tolower(srd$mapper) %||% "bwamem"

  ## short reads must have a "dir" a "platform" and and "mapper" field.
  fields0 <- c("dir", "platform", "mapper")
  srd <- compact(srd[fields0])
  assert_that(
    all(fields0 %in% names(srd)),
    msg = paste0("Missing fields <",
                 comma(fields0[!fields0 %in% names(srd)]),
                 "> in shortreads config"))

  ## Platform must be "illumina"
  assert_that(
    srd$platform %in% "illumina",
    msg = paste0("Unsupported longread type <", srd$platform, ">")
  )

  ## Mapper must be one of "bwamem" or "minimap"
  assert_that(
    srd$mapper %in% c("bwamem", "minimap"),
    msg = paste("Unsupported longread mapper <", srd$mapper, ">")
  )

  srd
}

validateDR2SConf <- function(conf) {
  fields <- MANDATORY_CONF_FIELDS()
  assert_that(all(fields %in% names(conf)),
              msg = paste("Missing fields <",
                          comma(fields[!fields %in% names(conf)]),
                          "> in config"))
  conf <- structure(compact(conf[ORDERED_CONF_FIELDS()]), class = c("DR2Sconf", "list"))
  conf$datadir <- normalizePath(conf$datadir, mustWork = TRUE)
  conf$outdir  <- normalizePath(conf$outdir, mustWork = FALSE)

  # Assert polymorphism threshold
  assert_that(
    is.numeric(conf$threshold),
    conf$threshold >= 0,
    conf$threshold <= 1,
    msg = "The polymorphism threshold must be a number between 0 ans 1")

  # Assert number of mapIter iterations
  assert_that(
    is.count(conf$iterations),
    conf$iterations > 0,
    conf$iterations < 10,
    msg = "The number of mapIter() iterations must fall between 1 and 10"
  )

  # Assert logicals
  assert_that(
    is.logical(conf$filterScores),
    is.logical(conf$microsatellite),
    is.logical(conf$forceMapping)
  )

  # Assert number of distinct alleles
  assert_that(is.count(conf$distAlleles),
              msg = "Number of distinct alleles (distAlleles) must be numeric")

  ## Assert pipeline
  pipesteps <- c("clear", "cache", "mapInit", "mapIter", "mapFinal",
                 "partitionLongReads", "partitionShortReads", "polish",
                 "report")
  assert_that(all(conf$pipeline %in% pipesteps),
              msg = paste(
                "Invalid pipeline step <",
                comma(conf$pipeline[!conf$pipeline %in% pipesteps]), ">") )

  ## Assert longreads
  conf$longreads <- normaliseLongreads(conf$longreads)

  ## Assert that longreads are available
  if (!is.null(conf$longreads$dir)) {
    longreaddir <- file.path(conf$datadir, conf$longreads$dir)
    assert_that(is.dir(longreaddir))
    assert_that(length(dir(longreaddir, pattern = ".+fast(q|a)(.gz)?")) > 0,
                msg = "No fasta or fastq file in the longread directory")
  } else {
    longreadfile <- file.path(conf$datadir, conf$longreads$file)
    assert_that(file.exists(longreadfile), msg = "No longreads available")
  }

  ## Check short reads (short reads are optional i.e. the field can be NULL)
  if (!is.null(conf$shortreads)) {
    ## Assert shortreads
    conf$shortreads <- normaliseShortreads(conf$shortreads)
    assert_that(!is.null(conf$shortreads$dir))
    shortreaddir <- file.path(conf$datadir, conf$shortreads$dir)
    assert_that(is.dir(shortreaddir))
    assert_that(length(dir(shortreaddir, pattern = ".+fast(q|a)(.gz)?")) > 0,
                msg = "No fasta or fastq file in the shortread directory"
    )
  } else {
    ## set shortreads field in config to NULL
    conf$shortreads = NULL
  }

  ## Normalise locus
  conf$locus <- sub("^HLA-", "", .normaliseLocus(conf$locus))

  conf
}


# Initialiser -------------------------------------------------------------


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

