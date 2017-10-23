
# reexports ---------------------------------------------------------------

#' @importFrom magrittr %>%
#' @export
magrittr::`%>%`

# Globals -----------------------------------------------------------------

DNA_BASES <- function() {
  c("A", "C", "G", "T", "a", "c", "g", "t")
}
DNA_BASES_UPPER <- function() {
  c("A", "C", "G", "T")
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
    A =  "A",  C = "C",   G = "G",   T = "T",    M = "AC",   R = "AG",   W = "AT",
    S = "CG",  Y = "CT",  K = "GT",  V = "ACG",  H = "ACT",  D = "AGT",  B = "CGT",  N = "ACGT",
    a = "-A",  c = "-C",  g = "-G",  t = "-T",   m = "-AC",  r = "-AG",  w = "-AT",
    s = "-CG", y = "-CT", k = "-GT", v = "-ACG", h = "-ACT", d = "-AGT", b = "-CGT", n = "-ACGT"
  )
}

DNA_PROB <- function(include = "indels"){
  if (include == "indels") {
    return(c(A = 0.2, C = 0.2, G = 0.2, T = 0.2, `-` = 0.20, `+` = 0.001))
  } else if (include == "ins"){
    return(c(A = 0.25, C = 0.25, G = 0.25, T = 0.25, `+` = 0.001))
  } else if (include == "del"){
    return(c(A = 0.25, C = 0.25, G = 0.25, T = 0.25, `-` = 0.001))
  }
  return(c(A = 0.25, C = 0.25, G = 0.25, T = 0.25))
}
PARTCOL <- function() {
  # colors from https://rdrr.io/cran/igraph/src/R/palette.R; extendable
  c(C = "#B35806", A = "#F1A340", E = "#FEE0B6", `-` = "#F7F7F7", F = "#D8DAEB", B = "#998EC3",
  D  = "#542788", N = "#9F9F9F" )
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

HLA_LOCI <- function() {
  c("A", "B", "C", "DRB1", "DRB1", "DQB1", "DPB1")
}

KIR_LOCI <- function() {
  c("2DL1", "2DL2", "2DL3", "2DL3", "2DL4", "2DL5", "2DS1", "2DS2", "2DS3",
    "2DS4", "2DS5", "2DP1", "3DL1", "3DL2", "3DL3", "3DS1", "3DP1")
}

# Helpers -----------------------------------------------------------------

normalise_locus <- function(locus) {
  locus <- sub("(HLA[-_]?|KIR[-_]?)", "", toupper(locus))
  if (locus %in% HLA_LOCI()) {
    paste0("HLA-", locus)
  } else if (locus %in% KIR_LOCI()) {
    paste0("KIR", locus)
  } else if (locus == "ABO") {
    locus
  } else {
    stop("Unknown locus ", sQuote(locus), call. = FALSE)
  }
}

expand_hla_allele <- function(allele, locus) {
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
  } else allele
}

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
  unlist(.mapply(function(x, n) paste0(x[n], collapse = collapse), list(x = xs, n = n), NULL))
}

optstring <- function(opts, ...) {
  opts[vapply(opts, isTRUE, FALSE)] <- ""
  paste0(c(sprintf("-%s%s", names(opts), opts), ...), collapse = " ") %|ch|% "default"
}

compact <- function(x) {
  x[!vapply(x, is.null, FALSE, USE.NAMES = FALSE)]
}

usc <- function(x) {
  gsub("[*:?<>|]", "_", x)
}

comma <- function(...) {
  paste0(..., collapse = ", ")
}

merge_list <- function (x, y, update = FALSE)  {
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

colon <- function(...) paste0(..., collapse = ":")

wrap <- function (x, wrap = "\"") {
  stopifnot(is.vector(x))
  sprintf("%s%s%s", wrap, x, wrap)
}

#' @keywords internal
#' @export
minimum <- function(n, m) {
  if (n > m) n else m
}

#' @keywords internal
#' @export
maximum <- function(n, m) {
  if (n < m) n else m
}

det2 <- function(m) {
  m[1, 1] * m[2, 2] - m[1, 2] * m[2, 1]
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

dir_create_if_not_exists <- function(path) {
  path <- normalizePath(path, mustWork = FALSE)
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
  }
  if (startsWith(path, "./")) {
    path <- file.path(getwd(), path)
  }
  normalizePath(path, mustWork = TRUE)
}

file_delete_if_exists <- function(path) {
  path <- normalizePath(path, mustWork = FALSE)
  if (file.exists(path)) {
    unlink(path)
  }
  invisible(path)
}

recode_fastq_header <- function(fqpath) {
  fq <- ShortRead::readFastq(fqpath)
  ids <- as.character(ShortRead::id(fq))
  if (any(grepl(" MD:Z:", ids))) {
    return(TRUE)
  }
  read_id <- sub("(;|\\s+).+$", "", ids)
  mdata <- sub("^[^; ]+[; ]+", "", ids)
  #mdata <- dkms::strsplitN(mdata, " ", 1)
  #mdata <- sub("MD:Z:", "", mdata)
  mdata <- if (all(grepl(";", mdata))) {
    paste0(dkms::strsplitN(mdata, " ", 1), paste0(";BARCODE=", gsub("(\\]|\\[)", "", dkms::strsplitN(mdata, " ", 2))))
  } else if (all(grepl("ONBC", mdata))) {
    paste0("BARCODE=", gsub("(\\]|\\[)", "", mdata))
  } else {
    gsub("(\\s+|:)", ";", mdata)
  }
  ids <- Biostrings::BStringSet(paste0(read_id, " MD:Z:", mdata))
  fq@id <- ids
  fq2 <- paste0(fqpath, "~")
  ShortRead::writeFastq(fq, file = fq2, compress = FALSE)
  if (file.copy(fq2, fqpath, overwrite = TRUE)) {
    unlink(fq2)
    return(TRUE)
  }
  FALSE
}

file_create_if_not_exists <- function(file) {
  path <- file.path(
    dir_create_if_not_exists(dirname(file)),
    basename(file)
  )
  if (!file.exists(path)) {
    file.create(path)
  }
  normalizePath(path, mustWork = TRUE)
}

has_command <- function(cmd) {
  stopifnot(assertthat::is.string(cmd))
  unname(Sys.which(cmd) != "")
}

subl <- function(x, pos = NULL) {
  assertthat::assert_that(has_command('subl'))
  if (tryCatch(assertthat::is.readable(x), assertError = function(e) FALSE)) {
    x <- normalizePath(x, mustWork = TRUE)
    if (!is.null(pos)) {
      x <- paste0(x, ":", pos)
    }
    system(paste0("subl ", x))
  } else {
    tmp <- tempfile()
    write(x, file = tmp)
    system(paste0("subl ", tmp))
  }
}


# Alignment browser -------------------------------------------------------


#' @export
browse_seqs <- function(seq,
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

#' @export
browse_align <- function(seq,
                         file = tempfile(fileext = ".html"),
                         openURL = TRUE,
                         patterns = CODE_PATTERN(),
                         colors = COL_PATTERN(),
                         colWidth = 120,
                         gapOpening = -14,
                         gapExtension = -2,
                         ...) {
  dna_seq <- Biostrings::DNAStringSet(gsub("+", "N", seq, fixed = TRUE))
  aln <- DECIPHER::AlignSeqs(
    dna_seq,
    verbose = FALSE,
    gapOpening = gapOpening,
    gapExtension = gapExtension,
    ...
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

#' @export
plot_diagnostic_alignment <- function(x, onlyFinal = FALSE) {

  # Given Ref
  seqs1 <- x$getRefSeq()
  names(seqs1) <- paste0("0 ", names(seqs1))
  # References from multialign or ref; and mapIter steps
  seqs2 <- Biostrings::DNAStringSet(unlist(lapply(x$mapIter, function(y) sapply(y, function(a) unlist(a$conseq)))))
  names(seqs2) <- unlist(lapply(names(x$mapIter), function(y) paste(x$getHapTypes(), "map", y)))

  # final reference
  seqs3 <- Biostrings::DNAStringSet(lapply(x$mapFinal$seq, function(y) unlist(y)))
  names(seqs3) <- paste(names(seqs3), "mapFinal")

  seqs <- c(seqs1, seqs2, seqs3)
  seqs <- seqs[order(names(seqs))]
  if (onlyFinal){
    DECIPHER::BrowseSeqs(DECIPHER::AlignSeqs(seqs4), colWidth = 120)
  }else {
    DECIPHER::BrowseSeqs(DECIPHER::AlignSeqs(seqs), colWidth = 120)
  }
}

