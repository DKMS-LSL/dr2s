
# reexports ---------------------------------------------------------------

#' @importFrom magrittr %>%
#' @export
magrittr::`%>%`


# HLA helpers -------------------------------------------------------------


.ipdHla <- function() {
  ipdHlaVersion <- getOption("dr2s.ipdHlaVersion")
  #ipdDbVersion = NULL
  if (exists("ipdHlaDb", envir = globalenv())) {
    db <- get("ipdHlaDb", envir = globalenv())
    if (!is.null(ipdHlaVersion)) {
      if (ipdHlaVersion == attr(db, "ipdHlaVersion"))
        return(db)
    } else {
      return(db)
    }
  }
  if (!is.null(ipdHlaVersion)) {
    if (endsWith(ipdHlaVersion, ".sqlite")) {
      assertthat::assert_that(assertthat::is.readable(ipdHlaVersion))
      db <- suppressMessages(ipdDb::loadLocalDataBase(ipdHlaVersion))
      attr(db, "ipdHlaVersion") <- ipdHlaVersion
      assign("ipdHlaDb", db,
             envir = globalenv())
      return(db)
    }
  }
  db <- suppressMessages(ipdDb::loadHlaData(ipdHlaVersion))
  attr(db, "ipdHlaVersion") <- ipdHlaVersion
  assign("ipdHlaDb", db,
         envir = globalenv())
  db
}

.ipdKir <- function() {
  ipdKirVersion <- getOption("dr2s.ipdKirVersion")
  #ipdDbVersion = NULL
  if (exists("ipdKirDb", envir = globalenv())) {
    db <- get("ipdKirDb", envir = globalenv())
    if (!is.null(ipdKirVersion)) {
      if (ipdKirVersion == attr(db, "ipdKirVersion"))
        return(db)
    } else {
      return(db)
    }
  }
  if (!is.null(ipdKirVersion)) {
    if (endsWith(ipdKirVersion, ".sqlite")) {
      assertthat::assert_that(assertthat::is.readable(ipdKirVersion))
      db <- suppressMessages(ipdDb::loadLocalDataBase(ipdKirVersion))
      attr(db, "ipdKirVersion") <- ipdKirVersion
      assign("ipdKirDb", db,
             envir = globalenv())
      return(db)
    }
  }
  db <- suppressMessages(ipdDb::loadKirData(ipdKirVersion))
  attr(db, "ipdKirVersion") <- ipdKirVersion
  assign("ipdKirDb", db,
         envir = globalenv())
  db
}

.normaliseLocus <- function(locus) {
  locus <- gsub(" ", "", locus)
  locus <- sub("(HLA[-_]?|KIR[-_]?)", "", toupper(locus))
  if (locus %in% HLA_LOCI()) {
    if (startsWith(locus, "MIC"))
      return(locus)
    "HLA-" %<<% locus
  } else if (locus %in% KIR_LOCI()) {
    "KIR" %<<% locus
  } else if (locus == "ABO") {
    locus
  } else {
    stop("Unknown locus ", sQuote(locus), call. = FALSE)
  }
}

.expandAllele <- function(allele, locus) {
  locus <- .normaliseLocus(locus)
  locus <- sub("KIR", "", toupper(locus))
  locus <- sub("HLA-", "", toupper(locus))
  if (locus %in% HLA_LOCI()) {
    pattern1 <- "^HLA-" %<<% locus %<<% "[*]\\d\\d\\d?(:.+)*$"
    pattern2 <- "^" %<<% locus %<<% "[*]\\d\\d\\d?(:.+)*$"
    pattern3 <- "^\\d\\d\\d?(:.+)*$"
    if (grepl(pattern1, allele)) {
      allele
    } else if (grepl(pattern2, allele)) {
      "HLA-" %<<% allele
    } else if (grepl(pattern3, allele)) {
      "HLA-" %<<% locus %<<% "*" %<<% allele
    }
  } else if (locus %in% KIR_LOCI()) {
    pattern1 <- "^KIR" %<<% locus %<<% "[*]\\d+$"
    pattern2 <- "^" %<<% locus %<<% "[*]\\d+$"
    pattern3 <- "^\\d+$"
    if (grepl(pattern1, allele)) {
      allele
    } else if (grepl(pattern2, allele)) {
      "KIR" %<<% allele
    } else if (grepl(pattern3, allele)) {
      "KIR" %<<% locus %<<% "*" %<<% allele
    }
  } else allele
}


# utilities ---------------------------------------------------------------


strsplit1 <- function(...) strsplit(...)[[1]]

strsplitN <- function(x, split, n, from = "start", collapse = split, ...) {
  stopifnot(is.vector(x))
  from <- match.arg(from, c("start", "end"))
  xs <- strsplit(x, split, ...)
  end <- vapply(xs, length, 0L)
  if (from == "end") {
    end <- end + 1L
    n <- lapply(end, `-`, n)
    n <- .mapply(`[<-`, list(x = n, i = lapply(n, `<`, 0),
                             value = 0L), NULL)
  }
  else {
    n <- lapply(rep.int(0, length(xs)), `+`, n)
    n <- .mapply(`[<-`, list(x = n, i = Map(`>`, n, end),
                             value = end), NULL)
  }
  n <- lapply(n, sort %.% unique)
  unlist(.mapply(function(x, n) paste0(x[n], collapse = collapse),
                 list(x = xs, n = n), NULL))
}

optstring <- function(opts) {
  opts <- list()
  opts[vapply(opts, isTRUE, FALSE)] <- ""
  paste(sprintf("%s%s", names(opts), opts), collapse = "")
}

compact <- function(x) {
  x[!vapply(x, is.null, FALSE, USE.NAMES = FALSE)]
}

indentation <- function(level = 1, shift = 2) {
  structure(function() {
    paste0(rep(" ", level*shift), collapse = "")
  }, class = c("indenter", "function"))
}

lvl <- function(x) UseMethod("lvl")
lvl.indenter <- function(x) {
  get("level", environment(x))
}

incr <- function(x, i) UseMethod("incr")
incr.indenter <- function(x, i = 1) {
  indentation(lvl(x) + i)
}

# Strip illegal characters from filenames a
strip <- function(x, replace = "") {
  gsub("[\\*:?<>|/]", replace, as.character(x))
}

# Strip insertion and deletion character from sequence
stripIndel <- function(x, replace = "") {
  gsub("[+-]", replace, as.character(x))
}

underscore <- function(x) gsub("\\s+", "_", x)

comma <- function(...) paste0(..., collapse = ", ")

semicolon <- function(...) paste0(..., collapse = ";")

colon <- function(...) paste0(..., collapse = ":")

dot <- function(...) paste0(..., collapse = ".")

litQuote <- function(x) paste0("\"", x, "\"")

litArrows <- function(x) paste0("<", x, ">")

mergeList <- function(x, y, update = FALSE) {
  if (length(x) == 0)
    return(y)
  if (length(y) == 0)
    return(x)
  if (update) {
    x[names(y)] <- y
  } else {
    i <- is.na(match(names(y), names(x)))
    if (any(i)) {
      x[names(y)[which(i)]] <- y[which(i)]
    }
  }
  x
}
## collapse integer vector to ranges
in.seq <- function(x) {
  # returns TRUE for elments within ascending sequences
  (c(diff(x, 1), NA) == 1 & c(NA, diff(x,2), NA) == 2)
}

contractSeqs <-  function(x) {
  # returns string formatted with contracted sequences
  x[in.seq(x)] <- ""
  gsub(",{2,}", "-", paste(x, collapse=","), perl=TRUE)
}

wrap <- function(x, wrap = "\"") {
  stopifnot(is.vector(x))
  sprintf("%s%s%s", wrap, x, wrap)
}

rescale <- function(x, lower, upper) {
  min <- min(x)
  max <- max(x)
  (((upper - lower)*(x - min))/(max - min)) + lower
}


`%||%` <- function(a, b) {
  if (length(a) == 0) b else a
}

`%|na|%` <- function(a, b) {
  ifelse(is.na(a), b, a)
}

`%|ch|%` <- function(a, b) {
  if (!nzchar(a)) b else a
}

`%<<%` <- function(a, b) {
  paste0(as.character(a), as.character(b))
}

`%+%` <- function(a, b) {
  a <- as.character(a)
  b <- as.character(b)
  if (nzchar(a)) paste0(a, b) else a
}

`%.%` <- function(f, g) {
  f <- match.fun(f)
  g <- match.fun(g)
  function(...) f(g(...))
}

.dirCreateIfNotExists <- function(path) {
  assert_that(is.character(path))
  path <- normalizePath(path, mustWork = FALSE)
  vapply(path, function(p) {
    if (!dir.exists(p)) {
      dir.create(p, recursive = TRUE)
    }
    if (startsWith(p, "./")) {
      path <- file.path(getwd(), p)
    }
    normalizePath(p, mustWork = TRUE)
  }, FUN.VALUE = character(1), USE.NAMES = FALSE)
}

.fileDeleteIfExists <- function(path) {
  path <- normalizePath(path, mustWork = FALSE)
  if (any(exP <- file.exists(path))) {
    unlink(path[exP])
  }
  invisible(path)
}

.cropPath <- function(base, path) {
  gsub("^/", "", gsub(base, "", unname(path), fixed = TRUE))
}

.hasCommand <- function(cmd) {
  assert_that(is.string(cmd))
  unname(Sys.which(cmd) != "")
}

## Mode
.getModeValue <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

#' Get the number of idling cores that can be used by calling the linux system
#' program mpstat
#' @param max_cpu_usage Maximum cpu usage (percentage idle time) for a core to
#' still be considered idle.
#'
#' @return An integer giving the number of idling cores
.getIdleCores <- function(max_cpu_usage = 10) {
  # total number of cores
  N_CORES <- parallel::detectCores()

  # if dr2s.max.cores is set, return whichever is smaller
  if (!is.null(MAX_CORES <- getOption("dr2s.max.cores"))) {
    return(min(MAX_CORES, N_CORES))
  }

  # fallback to N_CORES/4 if mpstat not installed
  if (!.hasCommand("mpstat")) {
    flog.warn("Install sysstat to make use of idle core detection", immediate. = TRUE)
    return(max(floor(N_CORES/4), 1))
  }

  # create list for readable lapply output
  cores <- lapply(seq_len(N_CORES), function(x) x - 1)
  names(cores) <- paste0('CPU', seq_len(N_CORES) - 1)
  # use platform specific system commands to get idle time
  proc_idle_time <- vapply(cores, function(x) {
    # assumes linux
    out <- system2(command = 'mpstat', args = c('-P', x), stdout = TRUE)
    as.double(sub(",", ".", unlist(strsplit(out[4], ' {2,}'))[12]))
  }, FUN.VALUE = double(1))

  max(sum(proc_idle_time > (70 - max_cpu_usage)) - 1, 1)
}

.setRunstats <- function(self, name, value) {
  assert_that(
    is.character(name),
    is.list(value)
  )
  value <- mergeList(self$getStats(name), value)
  self$setStats(name, list(value))
  writeDR2SConf(self)
  NULL
}

.coverage <- function(x, probs =  c(0.0, 0.5, 1.0)) {
  assert_that(is(x, "pileup"))
  quantile(
    .rowSums(consmat(x), NROW(consmat(x)), NCOL(consmat(x))),
    probs = probs)
}

.collectPileupParams <- function(...) {
  dots <- list(...)
  pParamList <- list(
    max_depth        = dots$max_depth %||% dots$maxDepth %||% 1e4,
    min_base_quality = dots$min_base_quality %||% dots$minBaseQuality %||% 3,
    min_mapq         = dots$min_mapq %||% dots$minMapq %||% 0,
    min_nucleotide_depth    = dots$min_nucleotide_depth %||% dots$minNucleotideDepth %||% 3,
    min_minor_allele_depth  = dots$min_minor_allele_depth %||% dots$minMinorAlleleDepth %||% 0,
    distinguish_strands     = dots$distinguish_strands %||% dots$distinguishStrands %||% FALSE,
    distinguish_nucleotides = dots$distinguish_nucleotides %||% dots$distinguishNucleotides %||% TRUE,
    ignore_query_Ns   = dots$ignore_query_Ns %||% dots$ignoreQueryNs %||% TRUE,
    include_deletions = dots$include_deletions %||% dots$includeDeletions %||% TRUE,
    include_insertions = dots$include_insertions %||% dots$includeInsertions %||% FALSE,
    left_bins  = dots$left_bins %||% dots$leftBins,
    query_bins = dots$query_bins %||% dots$queryBins,
    cycle_bins = dots$cycle_bins %||% dots$cycleBins
  )
  do.call(Rsamtools::PileupParam, pParamList)
}

editor <- function(x, pos = NULL, useEditor = "xdg-open") {
  useEditor <- match.arg(useEditor, c("xdg-open", "subl", "gvim", "gedit"))
  assert_that(.hasCommand(useEditor))
  if (tryCatch(is.readable(x), assertError = function(e) FALSE)) {
    x <- normalizePath(x, mustWork = TRUE)
    if (!is.null(pos) && useEditor == "subl") {
      x <- paste0(x, ":", pos)
    }
    system(paste(useEditor, x, sep = " "))
  } else {
    tmp <- tempfile()
    write(x, file = tmp)
    system(paste(useEditor, tmp, sep = " "))
  }
}


# Alignment browser -------------------------------------------------------


.browseSeqs <- function(seq,
                        file = tempfile(fileext = ".html"),
                        openURL = TRUE,
                        patterns = CODE_PATTERN(),
                        colors = COL_PATTERN(),
                        colWidth = 120,
                        ...) {
  DECIPHER::BrowseSeqs(
    seq,
    file,
    openURL,
    colorPatterns = TRUE,
    highlight = NA,
    patterns,
    colors,
    colWidth,
    ...
  )
}

.browseAlign <- function(seq,
                         file = tempfile(fileext = ".html"),
                         openURL = TRUE,
                         patterns = CODE_PATTERN(),
                         colors = COL_PATTERN(),
                         colWidth = 120,
                         gapOpening = -14,
                         gapExtension = -2,
                         ...) {
  dnaSeq <- Biostrings::DNAStringSet(gsub("+", "N", seq, fixed = TRUE))
  aln <- DECIPHER::AlignSeqs(
    dnaSeq,
    verbose = FALSE,
    gapOpening = gapOpening,
    gapExtension = gapExtension#,
    # ...
  )

  Biostrings::writeXStringSet(aln, paste0(file, ".fa"))
  DECIPHER::BrowseSeqs(
    aln,
    paste0(file, ".html"),
    openURL,
    colorPatterns = TRUE,
    highlight = NA,
    patterns,
    colors,
    colWidth
  )
}
#' Plot an alignment of all intermediate sequences and reference
#'
#' @param x A \code{\link[=DR2S_]{DR2S}} object.
#' @param onlyFinal restrict to final sequences
#'
#' @return NULL
#' @details Calls \code{\link[DECIPHER]{BrowseSeqs}} to open the alignment of
#' references, final and intermediate (only if \code{onlyFinal} is FALSE) in a
#' browser. The alignment is created with \code{\link[DECIPHER]{AlignSeqs}}.
#' @examples
#' ###
#' @export
plotDiagnosticAlignment <- function(x, onlyFinal = FALSE) {
  assert_that(is(dr2s, "DR2S"))
  # Given Ref
  seqs1 <- x$getRefSeq()
  names(seqs1) <- paste0("0 ", names(seqs1))

  seqs2 <- Biostrings::DNAStringSet(unlist(lapply(x$mapIter, function(y)
    lapply(y, function(a) unlist(a$conseq)))))
  names(seqs2) <- unlist(lapply(names(x$mapIter), function(y)
    paste(x$getHapTypes(), "map", y)))

  # final reference
  seqs3 <- Biostrings::DNAStringSet(lapply(x$mapFinal$seq, function(y)
    unlist(y)))
  names(seqs3) <- paste(names(seqs3), "mapFinal")

  seqs <- c(seqs1, seqs2, seqs3)
  seqs <- seqs[order(names(seqs))]
  if (onlyFinal) {
    DECIPHER::BrowseSeqs(DECIPHER::AlignSeqs(seqs3), colWidth = 120)
  } else {
    DECIPHER::BrowseSeqs(DECIPHER::AlignSeqs(seqs), colWidth = 120)
  }
}

