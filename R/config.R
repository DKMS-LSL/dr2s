#' @include configOpts.R
NULL

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
    format = conf0$format %||% "yaml",
    sampleId = sample,
    locus = locus,
    reference = reference,
    details = lapply(details, gsub, pattern = ";", replacement = ",") %||% NULL
  )
  structure(conf1, class = c("DR2Sconf", "list"))
}

#' Read a DR2S config file in yaml or json format
#' @param configFile The path to the valid DR2S config file or to a DR2S project
#' directrory containing a \file{config.json} or \file{config.yaml} file.
#' @param format Input format ("auto", "yaml" or "json")
#' @details DR2S config files can be created manually or by the
#' \code{\link{writeDR2SConf}} function.
#' @return A \code{DR2Sconf} object or a list of \code{DR2Sconf} objects.
#' @export
readDR2SConf <- function(configFile, format = "auto") {
  format0 <- match.arg(format, c("auto", "yaml", "json"))
  conf0 <- configFile
  if (format0 == "auto") {
    if (endsWith(tolower(configFile), ".yaml")) {
      format0 <- "yaml"
    } else if (endsWith(tolower(configFile), ".json")) {
      format0 <- "json"
    } else if (length(dir(configFile, pattern = "config\\.yaml")) == 1) {
      conf0 <- file.path(configFile, "config.yaml")
      format0 <- "yaml"
    } else if (length(conf0 <- dir(configFile, pattern = "config\\.json")) == 1) {
      conf0 <- file.path(configFile, "config.json")
      format0 <- "json"
    } else {
      stop("Config file format not recognised")
    }
  }
  conf <- switch(
    format0,
    yaml = yaml::yaml.load_file(con0),
    json = jsonlite::fromJSON(conf0)
  )
  ## set eferenc if extref exists
  if (!is.null(conf$extref)) {
    tryCatch(
      conf$reference <- normalizePath(conf$extref, mustWork = TRUE),
      error = function(e) {
        flog.error("External reference file <%s> not found", conf$extref, name = "info")
        stop("External reference file <", conf$extref, "> not found")
      })
  }
  ## set defaults if necessary
  conf["threshold"] <- conf$threshold %||% 0.2
  conf["iterations"] <- conf$iterations %||% 2
  conf["microsatellite"] <- conf$microsatellite %||% FALSE
  conf["filterScores"] <- conf$filterScores %||% TRUE
  conf["forceMapping"] <- conf$forceMapping %||% FALSE
  conf["distAlleles"] <- conf$distAlleles %||% 2
  conf["datadir"] <- conf$datadir %||% normalizePath(".", mustWork = TRUE)
  conf["outdir"] <- .cropOutdir(conf)
  conf["format"] <- if (format == "auto") format0 else format
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

.cropOutdir <- function(conf) {
  outdir <- conf$outdir %||% file.path(conf$datadir, "output")
  pathext <- gsub("//", "/", file.path("",
                                       conf$details$platform %||% "",
                                       conf$locus,
                                       conf$sampleId))
  sub(pathext, "", outdir, fixed = TRUE)
}

#' @export
print.DR2Sconf <- function(x, ...) {
  if (x$format == "yaml") {
    cat(yaml::as.yaml(x, indent = 4))
  } else if (x$format == "json") {
    print(jsonlite::toJSON(x, auto_unbox = TRUE, digits = 4, pretty = TRUE))
  }
  invisible(x)
}

#' Write the DR2S config as yaml or jsoon.
#' @param x A \code{\link{DR2S}} object.
#' @param outFile The destination file.
#' @param format Output format ("auto", "yaml" or "json")
#' @export
writeDR2SConf <- function(x, outFile = NULL, format = "auto") {
  assert_that(is(x, "DR2S"))
  format <- match.arg(format, c("auto", "yaml", "json"))
  if (format == "auto") {
    format <- match.arg(x$getConfig("format"), c("yaml", "json"))
  }
  if (is.null(outFile)) {
    outFile <- switch(
      format,
      yaml = file.path(x$getOutdir(), "config.yaml"),
      json = file.path(x$getOutdir(), "config.json")
    )
  }
  assert_that({
    if (format == "yaml")
      endsWith(outFile, ".yaml")
    else if (format == "json")
      endsWith(outFile, ".json")
  }, msg = paste0("Config file format <", format, "> does not match file extension"))
  conf <- x$getConfig()
  conf$runstats <- x$getStats()
  conf$runstats$reported <- x$getReportStatus()
  switch(
    format,
    yaml = yaml::write_yaml(conf, outFile),
    json = jsonlite::write_json(conf, path = outFile, auto_unbox = TRUE, pretty = TRUE)
  )
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
  conf0$outdir  <- normalizePath(conf0$outdir, mustWork = FALSE)
  conf0$longreads <- lrd
  conf0$sampleId  <- sampleId
  conf0$distAlleles <- conf0$distAlleles %||% sample$distAlleles
  sample$distAlleles <- NULL
  conf0["reference"] <- sample$reference %||% conf0$reference %||% list(NULL)
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
    "pipeline", "format"
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
    "pipeline", "opts", "format"
  )
}

normaliseLongreads <- function(lrd) {
  ## Set defaults
  lrd$type <- tolower(lrd$type) %||% "pacbio"
  lrd$mapper <- tolower(lrd$mapper) %||% "minimap"
  ## Long reads optional fields
  fields0 <- c("dir", "file", "type", "platform", "mapper", "hpc")
  lrd <- compact(lrd[fields0])
  ## Long reads mandatory fields
  fields1 <- c("file", "type", "mapper")
  fields2 <- c("dir", "type", "mapper")
  assert_that(
    all(fields1 %in% names(lrd)) || all(fields2 %in% names(lrd)),
    msg = paste0("Missing fields <",
                 comma(fields0[!fields0 %in% names(lrd)]),
                 "> in longreads config"))

  ## Type must be one of "pacbio" or "nanopore"
  assert_that(
    lrd$type %in% c("pacbio", "nanopore"),
    msg = paste0("Unsupported longread type <", lrd$type, ">")
  )

  ## Mapper must be one of "bwamem" or "minimap"
  assert_that(
    lrd$mapper %in% c("bwamem", "minimap"),
    msg = paste0("Unsupported longread mapper <", lrd$mapper, ">")
  )

  ## Homopolymer compression level 4, 3, 2, 1, 0
  if (!is.null(lrd$hpc)) {
    lrd$hpc <- as.integer(lrd$hpc)
    assert_that(
      lrd$hpc %in% 0:4,
      msg = paste0("Homopolymer compression level <", lrd$hpc, "> not supported")
    )
  }

  lrd
}

normaliseShortreads <- function(srd) {
  ## Set defaults
  srd$type   <- tolower(srd$type) %||% "illumina"
  srd$mapper <- tolower(srd$mapper) %||% "bwamem"

  ## short reads optional fields.
  fields0 <- c("dir", "type", "platform", "mapper")
  srd <- compact(srd[fields0])
  ## short reads mandatory fields
  fields1 <- c("dir", "type", "mapper")
  assert_that(
    all(fields1 %in% names(srd)),
    msg = paste0("Missing fields <",
                 comma(fields1[!fields1 %in% names(srd)]),
                 "> in shortreads config"))

  ## Platform must be "illumina"
  assert_that(
    srd$type %in% "illumina",
    msg = paste0("Unsupported longread type <", srd$type, ">")
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
                 "partitionLongreads", "partitionShortreads", "polish",
                 "report")
  assert_that(all(conf$pipeline %in% pipesteps),
              msg = paste(
                "Invalid pipeline step <",
                comma(conf$pipeline[!conf$pipeline %in% pipesteps]), ">") )

  ## Assert longreads
  conf$longreads <- normaliseLongreads(lrd = conf$longreads)

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

  ## Normalise options
  pipeline <- if (!is.null(conf$shortreads)) "SR" else "LR"
  conf$opts <- normaliseOpts(conf$opts, pipeline)

  conf
}


# Initialiser -------------------------------------------------------------


initialiseDR2S <- function(conf, createOutdir = TRUE) {
  conf <- validateDR2SConf(conf)
  if (createOutdir) {
    conf$outdir <- .dirCreateIfNotExists(path = gsub("//+", "/", file.path(
      conf$outdir,
      conf$details$platform %||% "",
      conf$locus,
      conf$sampleId
    )))
  }
  conf
}

