

# Consensus Matrix --------------------------------------------------------


#' Construct a consensus matrix from a pile-up
#'
#' @param x A \code{pileup} object.
#' @param freq If \code{TRUE} then frequencies are reported, otherwise counts.
#' @param ... Additional arguments such as \code{n}, \code{offset},
#' \code{insertions}, \code{read_ids}
#' @details
#' \code{consmat}: a \code{matrix} with positions row names and nucleotides as
#' column manes.
#' A \code{consmat} object includes the attributes:
#' \describe{
#'   \item{n}{<integer>; Number of reads per row.}
#'   \item{freq}{<logical>; Is it a frequency matrix or a count matrix.}
#'   \item{offset}{<integer>; Offset}
#'   \item{insertions}{<integer>; Insertions}
#'   \item{read_ids}{<character>; Read IDs}
#' }
#'
#' @return A \code{consmat} object.
#' @export
#' @examples
#' ###
consmat <- function(x, ...) UseMethod("consmat")

# Internal constructor
Consmat_ <- function(x, n, freq, offset = 0L, insertions = NULL, read_ids = NULL) {
  structure(
    x, n = n, freq = freq, offset = offset, insertions = insertions,
    read_ids = read_ids, class = c("consmat", "matrix")
  )
}

#' @export
consmat.matrix <- function(x, freq = TRUE, ...) {
  stopifnot(all(c("A", "C", "G", "T", "-", "+") %in% colnames(x)))
  n <- .rowSums(x, NROW(x), NCOL(x))
  x <- if (freq) sweep(x, 1, n, `/`) else x
  Consmat_(x, n, freq, ...)
}

#' @export
consmat.tbl_df <- function(x, freq = TRUE, drop.unused.levels = TRUE, ...) {
  stopifnot(all(c("pos", "nucleotide", "count") %in% colnames(x)))
  x <- xtabs(formula = count ~ pos + nucleotide, data = x, drop.unused.levels = drop.unused.levels)
  rs <- matrix(x, NROW(x),  NCOL(x))
  dimnames(rs) <- dimnames(x)
  n <- .rowSums(rs, NROW(rs), NCOL(rs))
  rs <- if (freq) sweep(rs, 1, n, `/`) else x
  Consmat_(rs, n, freq, ...)

}

#' @export
consmat.consmat <- function(x, freq = TRUE, ...) {
  if (freq) {
    if (is.freq(x)) {
      x
    } else {
      Consmat_(
        sweep(x, 1, n(x), `/`), n = n(x), freq = freq, offset = offset(x),
        insertions = ins(x), read_ids = ids(x)
      )
    }
  } else if (!freq) {
    ## recalibrate n
    n(x) <- .rowSums(x, NROW(x), NCOL(x))
    if (is.freq(x)) {
      Consmat_(
        sweep(x, 1, n(x), `*`), n = n(x), freq = freq, offset = offset(x),
        insertions = ins(x), read_ids = ids(x)
      )
    } else {
      x
    }
  }
}

#' @export
consmat.pileup <- function(x, freq = TRUE, ...) {
  cm <- if (is(x$consmat, "consmat")) x$consmat else consmat(x$pileup)
  consmat(cm, freq = freq)
}

#' @export
print.consmat <- function(x, n = 25, noHead = FALSE, transpose = FALSE,  ...) {
  if (!noHead) {
    cat("Consensus Matrix: ", NROW(x), " x ", NCOL(x), "\n", sep = "")
  }
  show <- if (transpose) {
    x <- as.matrix(t(x))
    if ((nc_ <- NCOL(x)) > n) {
      cbind(x[, 1:floor(n / 2)], x[, (nc_ - floor(n / 2)):nc_])
    } else x
  } else {
    x <- as.matrix(x)
    if (NROW(x) > n) {
      rbind(head(x, floor(n / 2)), tail(x, floor(n / 2)))
    } else x
  }
  print(show)
}

#' @export
`[.consmat` <- function(x, i, j, ..., drop = FALSE) {
  rs <- NextMethod(drop = drop)
  if (missing(j)) {
    ## recalibrate n
    n <- .rowSums(rs, NROW(rs), NCOL(rs))
    Consmat_(
      rs, n, freq = is.freq(x), offset = offset(x), insertions = ins(x),
      read_ids = compact(ids(x)[as.character(as.numeric(i) - offset(x))])
    )
  } else rs
}

#' @export
`[<-.consmat` <- function(x, i, j, ..., value) {
  rs <- NextMethod()
  ## recalibrate n
  n <- .rowSums(rs, NROW(rs), NCOL(rs))
  Consmat_(
    rs, n, freq = is.freq(x), offset = offset(x), insertions = ins(x),
    read_ids = compact(ids(x))
  )
}

#' @export
as.matrix.consmat <- function(x, ...) {
  rs <- unclass(x)
  attr(rs, "n") <- NULL
  attr(rs, "insertions") <- NULL
  attr(rs, "freq") <- NULL
  attr(rs, "offset") <- NULL
  attr(rs, "read_ids") <- NULL
  rs
}

#' @export
as.data.frame.consmat <- function(x, ...) {
  df <- dplyr::tbl_df(as.data.frame.table(x)) %>%
    dplyr::transmute(
      pos = as.integer(as.character(pos)),
      nucleotide = nucleotide,
      freq = Freq
    ) %>%
    dplyr::arrange(pos, nucleotide)
  df
}

#' @keywords internal
#' @export
n <- function(x, ...) UseMethod("n")
#' @export
n.consmat <- function(x) attr(x, "n")

#' @keywords internal
#' @export
`n<-` <- function(x, value, ...) UseMethod("n<-")
#' @export
`n<-.consmat` <- function(x, value) {
  attr(x, "n") <- value
  x
}

#' @keywords internal
#' @export
ins <- function(x, ...) UseMethod("ins")
#' @export
ins.consmat <- function(x) attr(x, "insertions")

#' @keywords internal
#' @export
`ins<-` <- function(x, value, ...) UseMethod("ins<-")
#' @export
`ins<-.consmat` <- function(x, value) {
  attr(x, "insertions") <- value
  x
}

#' @keywords internal
#' @export
ids <- function(x, ...) UseMethod("ids")
#' @export
ids.consmat <- function(x) attr(x, "read_ids")

#' @keywords internal
#' @export
`ids<-` <- function(x, value, ...) UseMethod("ids<-")
#' @export
`ids<-.consmat` <- function(x, value) {
  attr(x, "read_ids") <- value
  x
}

#' @keywords internal
#' @export
offset <- function(x, ...) UseMethod("offset")
#' @export
offset.consmat <- function(x) attr(x, "offset")

#' @keywords internal
#' @export
`offset<-` <- function(x, value, ...) UseMethod("offset<-")
#' @export
`offset<-.consmat` <- function(x, value) {
  attr(x, "offset") <- value
  x
}

#' @keywords internal
#' @export
is.freq <- function(x, ...) UseMethod("is.freq")
#' @export
is.freq.consmat <- function(x) attr(x, "freq")

#' @keywords internal
#' @export
prune_consensus_matrix <- function(cm, n_look_behind = 36, cutoff = 0.6, verbose = TRUE) {
  rowsum <- rowSums(cm[, 1:4]) ## only consider bases
  m0 <- do.call(cbind, Map(function(n) dplyr::lag(rowsum, n), n_look_behind:1))
  wquant <- suppressWarnings(apply(m0, 1, quantile, probs = 0.75, na.rm = TRUE))
  devi <- (wquant - rowsum)/wquant
  i <- which(devi > cutoff)
  if (length(i) > 0) {
    if (verbose) {
      message("  Pruning positions ", comma(i))
    }
    cm[i, ] <- 0L
    attr(cm, "pruning_positions") <- unname(i)
  }
  cm
}

#' @export
# Create position weight matrix
create_PWM <- function(msa){
  # Need to calc first a count based consensus matrix, while removing "+". Prob is calculated afterwards.
  cmat <- as.matrix(as.matrix(Biostrings::consensusMatrix(msa, as.prob = FALSE))[VALID_DNA(include = "del"),])
  cmat <- sweep(cmat, 2, colSums(cmat), "/")
  ## Add pseudocount
  cmat <- cmat + 1/length(msa)
  ## Divide by DNA_PROBABILITIES
  cmat <- cmat / DNA_PROB(include = "del")
  ## Get log2likelihood ratio
  cmat <- log2(cmat)
  cmat <- rbind(cmat, "+" = 0)
  cmat
}




# Position specific distance matrix
## Better use cpp_PSDM
# PSDM <- function(xseqs, xcons){
#   xseqsmat <- as.matrix(x_sub[1:4])
#   # upper <- unlist(sapply(1:nrow(xseqsmat), function(r) callPSDM(xcons, xseqsmat, r)))
#   upper <- unlist(sapply(1:3, function(r) callPSDM(xcons, xseqsmat, r)))
#   upper <- unlist(sapply(1:nrow(xseqsmat), function(r) callPSDM(xcons, xseqsmat, r)))
#   dist <- matrix(NA, nrow(xseqsmat), nrow(xseqsmat))
#   dist[lower.tri(dist, diag <- TRUE)] <- upper
#   rownames(dist) <- rownames(xseqsmat)
#   colnames(dist) <- rownames(xseqsmat)
#   as.dist(dist)
# }
#
# callPSDM <- function(xcons, xseqsmat, r){
#   l <- r:nrow(xseqsmat)
#   sapply(l, function(a) PSDM_cell(xcons, xseqsmat, r, a))
# }
#
# PSDM_cell <- function(xcons, xseqsmat, r, a){
#   xr <- xseqsmat[r,]
#   xa <- xseqsmat[a,]
#   b <- sapply(1:length(xr), function(x) PSDM_score(x, xcons, xr[x], xa[x]))
#   n <- length(xr) - sum(b == 1)
#   b[b == 1] <- 0
#   sum(b)/n
# }
#
# PSDM_score <- function(x, xcons, r, a){
#    if (r == a){
#      return(0)
#    }
#   else if(r == "-" || a == "-"){
#     return(1)
#   }
#   else {
#     return(xcons[r,x]*xcons[a,x])
#    }
# }
