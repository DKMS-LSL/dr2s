
# reexports ---------------------------------------------------------------

#' @importFrom magrittr %>%
#' @export
magrittr::`%>%`

# Globals -----------------------------------------------------------------

DNA_BASES <- function() {
  c("A", "C", "G", "T", "a", "c", "g", "t")
}

VALID_DNA <- function(include = "del"){
  include <- match.arg(include, c("none", "del", "ins", "indel"))
  if (include == "indel") return(c("G", "A", "T", "C", "-", "+"))
  if (include == "ins") return(c("G", "A", "T", "C", "+"))
  if (include == "del") return(c("G", "A", "T", "C", "-"))
  if (include == "none") return(c("G", "A", "T", "C"))
}

CODE_MAP <- function() {
  c(
    A =  "A",  C = "C",    G = "G",    T = "T",    M = "AC",   R = "AG",   
    W = "AT",  S = "CG",   Y = "CT",   K = "GT",   V = "ACG",  H = "ACT",  
    D = "AGT", B = "CGT",  N = "ACGT", a = "-A",   c = "-C",   g = "-G",  
    t = "-T",  m = "-AC",  r = "-AG",  w = "-AT",  s = "-CG",  y = "-CT", 
    k = "-GT", v = "-ACG", h = "-ACT", d = "-AGT", b = "-CGT", n = "-ACGT"
  )
}

DNA_PROB <- function(include = "indels") {
  if (include == "indels") {
    return(c(A = 0.2, C = 0.2, G = 0.2, T = 0.2, `-` = 0.2, `+` = 0.01))
  } else if (include == "ins") {
    return(c(A = 0.2, C = 0.2, G = 0.2, T = 0.2, `+` = 0.01))
  } else if (include == "del") {
    return(c(A = 0.2, C = 0.2, G = 0.2, T = 0.2, `-` = 0.2))
  }
  return(c(A = 0.25, C = 0.25, G = 0.25, T = 0.25))
}

PARTCOL <- function() {
  # colors from https://rdrr.io/cran/igraph/src/R/palette.R; extendable
  c(
    C = "#B35806",
    A = "#F1A340",
    E = "#FEE0B6",
    `-` = "#F7F7F7",
    F = "#D8DAEB",
    B = "#998EC3",
    D = "#542788",
    N = "#9F9F9F"
  )
  # c(A = "#f1a340", B = "#998ec3", `-` = "#f7f7f7")
}

NUCCOL <- function() {
  c(
    A = "#0087bd", ## BLUE
    T = "#ffd300", ## YELLOW
    G = "#009f6b", ## GREEN
    C = "#c40233", ## RED
    N = "#000000",
    `-` = "purple",
    `+` = "darkblue",
    `=` = "purple",
    " " = "grey80"
  )
}

CODE_PATTERN <- function() {
  c("[-]", "[A]", "[C]", "[G]", "[T]", "[RY]", "[SWKM]")
}

COL_PATTERN <- function() {
  c("#CC007A", "#CC2900", "#CCCC00", "#29CC00", "#00CC7A", "#007ACC", "#2900CC")
}

VALID_LOCI <- function() {
  hlaLoci <- ipdDb::loadHlaData()$getLoci()
  kirLoci <- ipdDb::loadKirData()$getLoci()
  c(hlaLoci, kirLoci)
}

HLA_LOCI <- function() {
  loci <- VALID_LOCI()
  hlaLoci <- loci[startsWith(loci, "HLA")]
  unname(sapply(hlaLoci, function(x) strsplit1(x, "-")[2]))
}

KIR_LOCI <- function() {
  loci <- VALID_LOCI()
  kirLoci <- loci[startsWith(loci, "KIR")]
  gsub(pattern = "KIR", "", kirLoci)
}

# Helpers -----------------------------------------------------------------

.normaliseLocus <- function(locus) {
  locus <- sub("(HLA[-_]?|KIR[-_]?)", "", toupper(locus))
  if (locus %in% HLA_LOCI()) {
    paste0("HLA-", locus)
  } else if (locus %in% KIR_LOCI()) {
    paste0("KIR", locus)
  } else if (locus == "ABO") {
    locus
  } else if (locus %in% c("MICA", "MICB")) {
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
    pattern1 <- paste0("^HLA-", locus, "[*]\\d\\d\\d?:.+$")
    pattern2 <- paste0("^", locus, "[*]\\d\\d\\d?:.+$")
    pattern3 <- "^\\d\\d\\d?:.+$"
    if (grepl(pattern1, allele)) {
      allele
    } else if (grepl(pattern2, allele)) {
      paste0("HLA-", allele)
    } else if (grepl(pattern3, allele)) {
      paste0("HLA-", locus, "*", allele)
    }
  } else if (locus %in% KIR_LOCI()) {
    pattern1 <- paste0("^KIR", locus, "[*]\\d+$")
    pattern2 <- paste0("^", locus, "[*]\\d+$")
    pattern3 <- "^\\d+$"
    if (grepl(pattern1, allele)) {
      allele
    } else if (grepl(pattern2, allele)) {
      paste0("KIR", allele)
    } else if (grepl(pattern3, allele)) {
      paste0("KIR", locus, "*", allele)
    }
  } else allele
}

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

optstring <- function(opts, ...) {
  opts[vapply(opts, isTRUE, FALSE)] <- ""
  paste0(c(sprintf("-%s%s", names(opts), opts), ...), collapse = " ") %|ch|% 
    "default"
}

compact <- function(x) {
  x[!vapply(x, is.null, FALSE, USE.NAMES = FALSE)]
}

usc <- function(x) gsub("[*:?<>|]", "_", x)

underscore <- function(x) gsub("\\s+", "_", x)

comma <- function(...) paste0(..., collapse = ", ")

semicolon <- function(...) paste0(..., collapse = ";")

colon <- function(...) paste0(..., collapse = ":")

dot <- function(...) paste0(..., collapse = ".")

litQuote <- function(x) paste0("\"", x, "\"")

litArrows <- function(x) paste0("<", x, ">")

.mergeList <- function(x, y, update = FALSE) {
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

wrap <- function(x, wrap = "\"") {
  stopifnot(is.vector(x))
  sprintf("%s%s%s", wrap, x, wrap)
}

minimum <- function(n, m) {
  if (n > m) n else m
}

maximum <- function(n, m) {
  if (n < m) n else m
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
  a <- as.character(a)
  b <- as.character(b)
  paste0(a, b)
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
  }, FUN.VALUE = character(1))
}

.fileDeleteIfExists <- function(path) {
  path <- normalizePath(path, mustWork = FALSE)
  if (file.exists(path)) {
    unlink(path)
  }
  invisible(path)
}


.cropPath <- function(base, path) {
  gsub("^/", "", gsub(base,"", path))
}

.hasCommand <- function(cmd) {
  assert_that(is.string(cmd))
  unname(Sys.which(cmd) != "")
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

## Mode
.getModeValue <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
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

  # Given Ref
  seqs1 <- x$getRefSeq()
  names(seqs1) <- paste0("0 ", names(seqs1))
  seqs2 <- Biostrings::DNAStringSet(unlist(lapply(x$mapIter, function(y) 
    sapply(y, function(a) unlist(a$conseq)))))
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

