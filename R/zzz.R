## allow setting maximum number of available cores via an environment
## variable
.onLoad <- function(libname, pkgname) {
  op <- options()
  op.dr2s <- list(
    #dr2s.max.cores = max(1L, parallel::detectCores() - 1L)
    dr2s.max.cores = NULL
  )
  toset <- !(names(op.dr2s) %in% names(op))
  if (any(toset)) {
    options(op.dr2s[toset])
  }
  invisible(NULL)
}

## Define global variables for vars in "foreach" and "data.table" to avoid
## R CMD check warnings
globalVariables(c("lrd", "h", "hp", "i", "pos", "sampleId", "nucleotide", "freq",
                  "count", "npoly", "clade", "prob", "haplotype", "read"))

## Define helper functions that emit global configuration variables
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

HLA_LOCI <- function() {
  loci <- .ipdHla()$getLoci()
  unname(vapply(loci, function(x) strsplit1(x, "-")[2],
                FUN.VALUE = character(1)))
}

KIR_LOCI <- function(ipd = NULL) {
  loci <- .ipdKir()$getLoci()
  gsub(pattern = "KIR", "", loci)
}

SR_PIPELINE <- function() {
  c("clear",
    "mapInit",
    "partitionLongreads",
    "mapIter",
    "partitionShortreads",
    "mapFinal",
    "polish",
    "report"
  )
}

LR_PIPELINE <- function() {
  c("clear",
    "mapInit",
    "partitionLongreads",
    "mapIter",
    "report",
    "cache"
  )
}
