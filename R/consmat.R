

# Consensus Matrix --------------------------------------------------------


#' Construct a consensus matrix from pileup.
#'
#' Accessor and replacement methods are available: n(consmat), ins(consmat),
#' offsetBases(consmat), is.freq(consmat).
#'
#' @param x A \code{pileup} object. Additional inputs can be a \code{matrix},
#' \code{tbl_df} or \code{consmat} objects.
#' @param freq If \code{TRUE} then frequencies are reported, otherwise counts.
#' @param ... Additional arguments such as \code{n}, \code{offsetBases},
#'   \code{insertions}
#' @param value The value to replace with.
#' @param consmat The \code{\link{consmat}} object that is changed.
#' @details
#' \code{consmat}: a \code{matrix} with positions row names and nucleotides as
#' column manes.
#' A \code{consmat} object includes the attributes:
#' \describe{
#'   \item{n}{<integer>; Number of reads per row.}
#'   \item{freq}{<logical>; Is it a frequency matrix or a count matrix.}
#'   \item{offsetBases}{<integer>; OffsetBases}
#'   \item{insertions}{<integer>; Insertions}
#' }
#' @return A \code{consmat} object.
#' @export
#' @examples
#' print("TODO: ADD EXAMPLES")
#' ###
consmat <- function(x, freq = TRUE, ...) UseMethod("consmat")

# Internal constructor
Consmat_ <- function(x, n, freq, offsetBases = 0L, insertions = NULL) {
  structure(
    x, n = n, freq = freq, offsetBases = offsetBases, insertions = insertions,
    class = c("consmat", "matrix")
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
consmat.tbl_df <- function(x, freq = TRUE, drop.unused.levels = FALSE, ...) {
  stopifnot(all(c("pos", "nucleotide", "count") %in% colnames(x)))
  x <- stats::xtabs(formula = count ~ pos + nucleotide, data = x,
             drop.unused.levels = drop.unused.levels)
  rs <- matrix(x, NROW(x),  NCOL(x))
  dimnames(rs) <- dimnames(x)
  n <- .rowSums(rs, NROW(rs), NCOL(rs))
  rs <- if (freq) sweep(rs, 1, n, `/`) else x
  rs <- rs[,VALID_DNA("indel")]
  Consmat_(rs, n, freq, ...)
}

#' @export
consmat.consmat <- function(x, freq = TRUE, ...) {
  if (freq) {
    if (is.freq(x)) {
      x
    } else {
      Consmat_(
        sweep(x, 1, n(x), `/`), n = n(x), freq = freq,
        offsetBases = offsetBases(x), insertions = ins(x)
      )
    }
  } else if (!freq) {
    ## recalibrate n
    n(x) <- .rowSums(x, NROW(x), NCOL(x))
    if (is.freq(x)) {
      Consmat_(
        sweep(x, 1, n(x), `*`), n = n(x), freq = freq,
        offsetBases = offsetBases(x), insertions = ins(x)
      )
    } else {
      x
    }
  }
}

#' @export
print.consmat <- function(x, n = 25, noHead = FALSE, transpose = FALSE,  ...) {
  if (!noHead) {
    cat("Consensus Matrix: ", NROW(x), " x ", NCOL(x), "\n", sep = "")
  }
  show <- if (transpose) {
    x <- as.matrix(t(x))
    if ((nc_ <- NCOL(x)) > n) {
      cbind(x[, seq_len(floor(n / 2))], x[, (nc_ - floor(n / 2)):nc_])
    } else x
  } else {
    x <- as.matrix(x)
    if (NROW(x) > n) {
      rbind(utils::head(x, floor(n / 2)), utils::tail(x, floor(n / 2)))
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
      rs, n, freq = is.freq(x), offsetBases = offsetBases(x),
      insertions = ins(x))
  } else rs
}

#' @export
`[<-.consmat` <- function(x, i, j, ..., value) {
  rs <- NextMethod()
  ## recalibrate n
  n <- .rowSums(rs, NROW(rs), NCOL(rs))
  Consmat_(
    rs, n, freq = is.freq(x), offsetBases = offsetBases(x), insertions = ins(x)
  )
}

#' @export
as.matrix.consmat <- function(x, ...) {
  rs <- unclass(x)
  attr(rs, "n") <- NULL
  attr(rs, "insertions") <- NULL
  attr(rs, "freq") <- NULL
  attr(rs, "offsetBases") <- NULL
  rs
}

#' @export
as.data.frame.consmat <- function(x, ...) {
  df <- tibble::as_tibble(as.data.frame.table(x)) %>%
    dplyr::transmute(
      pos = as.integer(as.character(.data$pos)),
      nucleotide = .data$nucleotide,
      freq = .data$Freq
    ) %>%
    dplyr::arrange(.data$pos, .data$nucleotide)
  df
}

#' @rdname consmat
#' @export
`consmat<-` <- function(x, value) UseMethod("consmat<-")

#' @rdname consmat
#' @export
n <- function(consmat) UseMethod("n")

#' @export
n.consmat <- function(consmat) attr(consmat, "n")

#' @rdname consmat
#' @export
`n<-` <- function(consmat, value) UseMethod("n<-")
#' @export
`n<-.consmat` <- function(consmat, value) {
  attr(consmat, "n") <- value
  consmat
}

#' @rdname consmat
#' @export
ins <- function(consmat) UseMethod("ins")
#' @export
ins.consmat <- function(consmat) attr(consmat, "insertions")

#' @rdname consmat
#' @export
`ins<-` <- function(consmat, value) UseMethod("ins<-")
#' @export
`ins<-.consmat` <- function(consmat, value) {
  attr(consmat, "insertions") <- value
  consmat
}


#' @rdname consmat
#' @export
offsetBases <- function(consmat) UseMethod("offsetBases")
#' @export
offsetBases.consmat <- function(consmat) attr(consmat, "offsetBases")

#' @rdname consmat
#' @export
`offsetBases<-` <- function(consmat, value) UseMethod("offsetBases<-")
#' @export
`offsetBases<-.consmat` <- function(consmat, value) {
  attr(consmat, "offsetBases") <- value
  consmat
}

#' @rdname consmat
#' @export
is.freq <- function(consmat) UseMethod("is.freq")
#' @export
is.freq.consmat <- function(consmat) attr(consmat, "freq")

## Note for me: flattens the matrix; compare rowsum to rowsums upstream of pos;
## if > t set all to 0
.pruneConsensusMatrix <- function(cm, nLookBehind = 36, cutoff = 0.6,
                                  verbose = TRUE) {
  ## Use names; by helper function
  rowsum <- rowSums(cm[, 1:4]) ## only consider bases
  m0 <- do.call(cbind, Map(function(n) dplyr::lag(rowsum, n), nLookBehind:1))
  wquant <- suppressWarnings(apply(m0, 1,
                                   stats::quantile,
                                   probs = 0.75, na.rm = TRUE))
  devi <- (wquant - rowsum)/wquant
  i <- which(devi > cutoff)
  if (length(i) > 0) {
    if (verbose) {
      message("  Pruning positions ", comma(i))
    }
    cm[i, ] <- 0L
    attr(cm, "pruningPositions") <- unname(i)
  }
  cm
}

#' Create position weight matrix from multiple sequence alignment
#'
#' @param msa A \code{DNAStringSet} object of aligned sequences.
#' @param indelRate <numeric|NULL>; The estimated background rate of Indels.
#' @details
#' \code{PWM}: a \code{matrix} with position as row names and nucleotides as
#' column manes. Values are nucleotide weights at a position
#' A ConsensusMatrix is calculated from the MSA using
#' \code{Biostrings::consensusMatrix} and values are converted to probabilities.
#' Pseudocounts are added and values are divided by DNA probabilities and log2
#' score is reported
#' @return A \code{PWM} matrix.
#'
#' @export
#' @examples
#' print("TODO: Add examples")
#' ###
createPWM <- function(msa, indelRate = NULL) {
  # Need to calc first a count based consensus matrix, while removing "+".
  # Prob is calculated afterwards.
  indent2 <- indentation(2)
  cmat <- as.matrix(Biostrings::consensusMatrix(msa, as.prob = FALSE))
  cmat <- cmat[VALID_DNA("del"), ]
  cmat <- sweep(cmat, 2, colSums(cmat), "/")
  ## Add pseudocount
  cmat <- cmat + 1/length(msa)
  ## Normalise
  cmat <- cmat/colSums(cmat)
  ## Divide by background model
  b <- DNA_PROB(gapfreq = indelRate, include = "del")
  if (!is.null(indelRate)) {
    flog.info("%sUsing background model <%s> to create PWM", indent2(),
              comma(sprintf("%s=%s", dQuote(names(b)), round(b, 4))),
              name = "info")
  }
  cmat <- cmat/b
  ## Get log2likelihood ratio
  cmat <- log2(cmat)
  cmat <- rbind(cmat, "+" = 0)
  cmat
}

# mat <- consmat(pileup)
#removeError = TRUE
.distributeGaps <- function(mat, bamfile, removeError = FALSE, ...) {
  indent <- list(...)$indent %||% indentation()
  seq <- .mat2rle(mat)
  if (removeError) {
    gapError <- .getGapErrorBackground(mat, n = 5)
    flog.info("%sEstimate indel noise <%0.3g> to suppress spurious gaps", indent(),
              gapError, name = "info")
  }
  workers <- min(sum(idx <- seq$length > 5), .getIdleCores())
  bpparam <- BiocParallel::MulticoreParam(workers = workers, log = FALSE)
  flog.info("%sUse %s workers to shift gaps at positions <%s>", indent(),
            workers, comma(which(idx)), name = "info")
  bam <- Rsamtools::BamFile(bamfile)
  if (workers > 1) {
    Rsamtools::open.BamFile(bam)
    on.exit(Rsamtools::close.BamFile(bam))
  }
  ## Collect changed matrix elements
  ## TODO: bplapply throws warning in serialize(data, node$con, xdr = FALSE)
  ## 'package:stats' may not be available when loading
  ## For the time being we suppress this warning
  changeMat <- suppressWarnings(BiocParallel::bplapply(which(idx), function(i, bfl) {
    #  i <- which(seq$length > 5)[2]
    ## Assign new gap numbers to each position starting from left
    #  meanCoverage <- mean(rowSums(mat[(seqStart-2):(seqEnd+2),1:4]))
    seqStart <- sum(seq$lengths[seq_len(i - 1)]) + 1
    seqEnd <- seqStart + seq$lengths[i] - 1
    range <- c(seqStart, seqEnd)
    msa <- .msaFromBam(bfl, range, paddingLetter = "+")
    ## Skip position if empty
    if (length(msa) == 0) {
      NULL
    } else {
      ## Use only sequences spanning the complete region! Every other sequence
      ## gives no Info
      msa <- msa[vapply(msa, function(x) !"+" %in% Biostrings::uniqueLetters(x),
                        FUN.VALUE = logical(1))]
      if (removeError) {
        aFreq <- Biostrings::alphabetFrequency(msa[1])
        nt <- colnames(aFreq)[which.max(aFreq)]
        nucs <- sum(Biostrings::nchar(msa))
        falseGaps <- round(nucs * gapError)
        while (falseGaps > 0) {
          seqsWithGap <- which(vapply(msa, function(x)
            "-" %in% Biostrings::uniqueLetters(x), FUN.VALUE = logical(1)))
          if (length(seqsWithGap) > 1) {
            set.seed(3)
            seqsToChange <- sample(seqsWithGap,
                                   min(falseGaps,
                                       length(seqsWithGap)))
            s1 <- seqsToChange
          } else {
            seqsToChange <- seqsWithGap
          }
          if (length(seqsToChange) == 0)
            break
          msaChanged <- Biostrings::DNAStringSet(
            lapply(msa[seqsToChange], function(x) {
              # x <- msa[seqsToChange][[1]]
              pos <- Biostrings::matchPattern("-", x)[1]
              Biostrings::replaceLetterAt(x,
                                          Biostrings::start(
                                            Biostrings::matchPattern("-", x)[1]),
                                          nt)
            }))
          msa[seqsToChange] <- msaChanged
          falseGaps <- falseGaps - length(seqsToChange)
        }
      }

      msa <- Biostrings::DNAStringSet(vapply(msa, function(x) {
        gapless <- gsub(pattern = "-", "", x)
        nGaps <- Biostrings::nchar(x) - Biostrings::nchar(gapless)
        paste0(rep("-", nGaps), collapse = "") %<<% gapless
      }, FUN.VALUE = character(1)))

      list(
        seqStart = seqStart,
        seqEnd = seqEnd,
        matPart = t(Biostrings::consensusMatrix(msa, as.prob = FALSE))[
          , VALID_DNA(include = "indel")]
      )
    }
    #}, bamfile = bamfile, BPPARAM = bpparam))
  }, bfl = bam, BPPARAM = bpparam))
  ##
  assert_that(!any(vapply(changeMat, is.null, FUN.VALUE = logical(1))))
  ## Change the matrix
  for (i in changeMat) {
    mat[i$seqStart:i$seqEnd, ] <- i$mat
  }
  mat
}

.getGapErrorBackground <- function(mat, n = 5){
  seq <- .mat2rle(mat)
  consecutiveRegions <- rle((seq$lengths > n))
  longestRegion <- which.max(consecutiveRegions$lengths)
  startSeq <- ifelse(longestRegion == 1,
                     1,
                     sum(consecutiveRegions$lengths[
                       seq_len(longestRegion - 1)]))
  endSeq <- startSeq + consecutiveRegions$lengths[longestRegion] - 1
  ## Add an offset of 10 to acknowledge bad quality after homopolymer regions
  startSeq <- startSeq + 10
  background <- mat[startSeq:endSeq,]
  backgroundSums <- colSums(background)
  backgroundSums["-"] / sum(backgroundSums)
}


.extractIdsFromMat <- function(mat, readIds) {
  consmat(t(
      Biostrings::consensusMatrix(mat[readIds], as.prob = FALSE)[
        VALID_DNA(include = "indel"), ]
    ), freq = FALSE)
}
