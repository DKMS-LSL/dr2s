
# reexports ---------------------------------------------------------------

#' @importFrom magrittr %>%
#' @export
magrittr::`%>%`

# Globals -----------------------------------------------------------------

DNA_BASES <- function() {
  c("A", "C", "G", "T", "a", "c", "g", "t")
}

CODE_MAP <- function() {
  c(
    A =  "A",  C = "C",   G = "G",   T = "T",    M = "AC",   R = "AG",   W = "AT",
    S = "CG",  Y = "CT",  K = "GT",  V = "ACG",  H = "ACT",  D = "AGT",  B = "CGT",  N = "ACGT",
    a = "-A",  c = "-C",  g = "-G",  t = "-T",   m = "-AC",  r = "-AG",  w = "-AT",
    s = "-CG", y = "-CT", k = "-GT", v = "-ACG", h = "-ACT", d = "-AGT", b = "-CGT", n = "-ACGT"
  )
}

PARTCOL <- function() {
  c(A = "#f1a340", B = "#998ec3", `-` = "#f7f7f7")
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
  DECIPHER::BrowseSeqs(
    aln,
    file,
    openURL,
    colorPatterns = TRUE,
    highlight = NA,
    patterns,
    colors,
    colWidth
  )
}

#' @export
plot_diagnostic_alignment <- function(x) {
  if (x$hasMultialign()) {
    A1 <- x$map1$A$ref$refseq
    names(A1) <- paste("A", x$getLrdType(), "multialign")
    A2 <- x$map1$A$conseq$multialign
    names(A2) <- paste("A", x$getLrdType(), "map1")
    B1 <- x$map1$B$ref$refseq
    names(B1) <- paste("B", x$getLrdType(), "multialign")
    B2 <- x$map1$B$conseq$multialign
    names(B2) <- paste("B", x$getLrdType(), "map1")
  } else {
    A1 <- B1 <- x$getRefSeq()
    names(A1) <- strsplitN(names(A1), "~", 1, fixed = TRUE)
    names(B1) <- strsplitN(names(B1), "~", 1, fixed = TRUE)
    A2 <- x$map1$A$conseq$reference
    names(A2) <- paste("A", x$getLrdType(), "map1")
    B2 <- x$map1$B$conseq$reference
    names(B2) <- paste("B", x$getLrdType(), "map1")
  }
  A3 <- Biostrings::BStringSet(x$map2$A$conseq[[1]])
  names(A3) <- paste("A", x$getLrdType(), "map2")
  A4 <- Biostrings::BStringSet(x$consensus$seq$HapA)
  names(A4) <- paste("A", x$getSrdType(), "map3")

  B3 <- Biostrings::BStringSet(x$map2$B$conseq[[1]])
  names(B3) <- paste("B", x$getLrdType(), "map2")
  B4 <- Biostrings::BStringSet(x$consensus$seq$HapB)
  names(B4) <- paste("B", x$getSrdType(), "map3")
  seqs <- Biostrings::DNAStringSet(c(A1, A2, A3, A4, B1, B2, B3, B4))
  DECIPHER::BrowseSeqs(DECIPHER::AlignSeqs(seqs), colWidth = 120)
}

