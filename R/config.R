create_dr2s_conf <- function(sample,
                             locus,
                             longreads = list(type = "pacbio", dir = "pacbio"),
                             shortreads = list(type = "illumina", dir = "illumina"),
                             datadir = ".",
                             outdir = "./output",
                             reference = NULL,
                             consensus = NULL,
                             threshold = 0.20,
                             iterations = 1,
                             fullname = TRUE,
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
    mapper     = conf0$mapper   %||% "bwamem",
    # limitA     = conf0$limitA   %||% NULL,
    # limitB     = conf0$limitB   %||% NULL,
    limits     = conf0$limits   %||% NULL,
    haptypes   = conf0$haptypes   %||% NULL,
    pipeline   = conf0$pipeline %||% c("clear", "mapInit", "partition", "split", "extract", "mapIter", "mapFinal", "polish", "report"),
    longreads  = longreads,
    shortreads = shortreads,
    nreads     = conf0$nreads   %||% NULL,
    opts       = conf0$opts     %||% NULL,
    sample_id  = sample,
    locus      = locus,
    reference  = ref,
    alternate  = alt,
    consensus  = consensus      %||% "multialign"
  )
  structure(conf1, class = c("DR2Sconf", "list"))
}

read_dr2s_conf <- function(config_file) {
  conf <- yaml::yaml.load_file(config_file)
  ## set defaults if necessary
  conf$datadir    <- conf$datadir    %||% normalizePath(".", mustWork = TRUE)
  conf$outdir     <- conf$outdir     %||% file.path(conf$datadir, "output")
  conf$threshold  <- conf$threshold  %||% 0.2
  conf$iterations <- conf$iterations %||% 1
  conf$mapper     <- conf$mapper     %||% "bwamem"
  # conf["limitA"]  <- conf$limitA     %||% list(NULL)
  # conf["limitB"]  <- conf$limitN     %||% list(NULL)
  conf$limits     <- conf$limits     %||% list(NULL)
  conf$haptypes   <- conf$haptypes   %||% list(NULL)
  conf$pipeline   <- conf$pipeline   %||% c("clear", "mapInit", "partition", "split", "extract", "mapIter", "mapFinal", "polish", "report")
  conf$longreads  <- conf$longreads  %||% list(type = "pacbio", dir = "pacbio")
  conf$shortreads <- conf$shortreads %||% list(type = "illumina", dir = "illumina")

  if (length(conf$shortreads) == 1 && is.list(conf$shortreads[[1]]))
    conf$shortreads <- conf$shortreads[[1]]

  if (is.null(conf$nreads))
    conf["nreads"] <-  list(NULL)

  if (is.null(conf$opts))
    conf["opts"] <-  list(NULL)

  structure(conf, class = c("DR2Sconf", "list"))
}

expand_dr2s_conf <- function(conf) {
  ## we can have more than one sample
  samples <- conf$samples
  conf$samples <- NULL
  sample_ids <- names(samples)
  ## we can have more than one longread type
  lrds <- conf$longreads
  if (!is.null(names(lrds))) {
    lrds <- list(lrds)
  }
  conf$longreads <- NULL
  ## we can have subsampling of longreads
  nrds <- conf$nreads %||% list(NULL)
  conf$nreads <- NULL
  ## we can have different consensus calling methods
  cnss <- conf$consensus %||% list("multialign")
  conf$consensus <- NULL
  foreach(sample = samples, sample_id = sample_ids, .combine = "c") %:%
    foreach(locus = sample, .combine = "c") %:%
    foreach(nrd = nrds, .combine = "c") %:%
    foreach(lrd = lrds, .combine = "c") %:%
    foreach(cns = cnss, .combine = "c") %:%
    foreach(ref = iter(locus$reference %||% ""), alt = iter(locus$alternate %||% ""), .combine = "c") %do% {
      update_conf(conf, lrd, nrd, sample_id, locus, ref, alt, cns)
    }
}

update_conf <- function(conf0, lrd, nrd, sample_id, locus, reference, alternate, consensus) {
  conf0$datadir   <- normalizePath(conf0$datadir, mustWork = TRUE)
  conf0$outdir    <- normalizePath(conf0$outdir, mustWork = FALSE)
  conf0$longreads <- lrd
  conf0['nreads'] <- nrd %||% list(NULL)
  conf0$sample_id <- sample_id
  conf0["reference"] <- reference %|ch|% list(NULL)
  locus$reference    <- NULL
  conf0["alternate"] <- alternate %|ch|% list(NULL)
  locus$alternate    <- NULL
  conf0["consensus"] <- consensus %|ch|% list(NULL)
  locus$consensus    <- NULL
  ## add overides if they exist
  list(merge_list(conf0, locus, update = TRUE))
}

initialise_dr2s <- function(conf, create_outdir = TRUE) {
  conf <- validate_dr2s_conf(conf)
  if (create_outdir) {
    conf$outdir <- dir_create_if_not_exists(path = gsub("//+", "/", file.path(
      conf$outdir,
      conf$sample_id,
      sprintf("N%03i", conf$nreads) %||% "",
      conf$longreads$name %||% conf$longreads$dir %||% "",
      paste0(normalise_locus(conf$locus), ".", conf$longreads$type, ".", conf$reftype, ".", conf$consensus)
    )))
  }
  conf
}

validate_dr2s_conf <- function(conf) {
  fields <- c("datadir", "outdir", "threshold", "iterations", "mapper", "limits", "haptypes",
              "pipeline", "longreads", "shortreads", "nreads", "opts",
              "sample_id", "locus", "reference", "alternate", "consensus")#"limitA", "limitB",
  if (!all(fields %in% names(conf))) {
    stop("Missing fields <", comma(fields[!fields %in% names(conf)]), "> in config", call. = FALSE)
  }
  conf <- structure(conf[fields], class = c("DR2Sconf", "list"))
  conf$datadir <- normalizePath(conf$datadir, mustWork = TRUE)
  conf$outdir  <- normalizePath(conf$outdir, mustWork = FALSE)
  if (!is.numeric(conf$threshold) || (conf$threshold <= 0 || conf$threshold >= 1)) {
    stop("The polymorphism threshold must be a number between 0 ans 1", call. = FALSE)
  }
  # check iterations
  if (!is.numeric(conf$iterations) || !conf$iterations&&1 ==  0 || conf$iterations > 10) {
    stop("The number of iterations must be integers between 1 and 10", call. = FALSE)
  }
  ## Check mapper
  conf$mapper <- match.arg(conf$mapper, c("bwamem", "graphmap"))
  ## Check pipeline
  pipesteps <- c("clear", "cache", "mapInit", "mapIter", "mapFinal",
                 "partition", "split", "extract", "polish", "report")
  if (!all(conf$pipeline %in% pipesteps)) {
    stop("Invalid pipeline step <", comma(conf$pipeline[!conf$pipeline %in% pipesteps]), ">", call. = FALSE)
  }
  ## Long reads must have a "type" and a "dir" field. "name" and "opts" are optional
  lrdfields <- c("type", "dir")
  if (!all(lrdfields %in% names(conf$longreads))) {
    stop("Missing fields <", comma(lrdfields[!lrdfields %in% names(conf$longreads)]), "> in longreads config", call. = FALSE)
  }
  ## Type must be one of "pacbio" or "nanopore"
  if (!conf$longreads$type %in% c("pacbio", "nanopore")) {
    stop("Unsupported longread type <", conf$longreads$type, ">", call. = FALSE)
  }
  longreaddir <- file.path(conf$datadir, conf$longreads$dir)
  if (!file.exists(longreaddir)) {
    stop("No long read directory <", longreaddir, ">", call. = FALSE)
  }
  if (length(dir(longreaddir, pattern = ".+fastq(.gz)?")) == 0) {
    stop("No FASTQ files in long read directory <", longreaddir, ">", call. = FALSE)
  }
  ## Check short reads (short reads are optional i.e. the field can be NULL)
  if (!is.null(conf$shortreads)) {
    srdfields <- c("type", "dir")
    if (!all(srdfields %in% names(conf$shortreads))) {
      stop("Missing fields <", comma(srdfields[!srdfields %in% names(conf$shortreads)]), "> in shortreads config", call. = FALSE)
    }
    if (!conf$shortreads$type %in% c("illumina")) {
      stop("Unsupported shortread type <", conf$shortreads$type, ">", call. = FALSE)
    }
    shortreaddir <- file.path(conf$datadir, conf$shortreads$dir)
    if (!file.exists(longreaddir)) {
      warning("No short read directory <", shortreaddir, ">", call. = FALSE, immediate. = TRUE)
    }
    if (length(dir(shortreaddir, pattern = ".+fastq(.gz)?")) == 0) {
      warning("No FASTQ files in short read directory <", shortreaddir, ">", call. = FALSE, immediate. = TRUE)
    }
  }
  ## Normalise locus
  ##conf$locus   <- sub("^(HLA-|KIR)", "", normalise_locus(conf$locus))
  conf$locus   <- sub("^HLA-", "", normalise_locus(conf$locus))
  ## Check consensus method
  conf$consensus <- match.arg(conf$consensus, c("multialign", "mapping"))
  if (conf$consensus == "multialign") {
    conf['alternate'] <- list(NULL)
  }
  ## Reference type
  conf$reftype <- if (is.null(conf$reference) && is.null(conf$alternate)) {
    "cons"
  } else if (!is.null(conf$reference) && is.null(conf$alternate)) {
    "ref"
  } else if (!is.null(conf$reference) && !is.null(conf$alternate)) {
    "ref.alt"
  }
  conf
}

#' @export
print.DR2Sconf <- function(x, ...) {
  cat(yaml::as.yaml(x, indent = 4))
  invisible(x)
}

