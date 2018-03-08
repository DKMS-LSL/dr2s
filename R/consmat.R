

# Consensus Matrix --------------------------------------------------------


#' Construct a consensus matrix from pileup
#'
#' @param x A \code{pileup} object.
#' @param freq If \code{TRUE} then frequencies are reported, otherwise counts.
#' @param ... Additional arguments such as \code{n}, \code{offset}, \code{insertions}
#' @details
#' \code{consmat}: a \code{matrix} with positions row names and nucleotides as
#' column manes.
#' A \code{consmat} object includes the attributes:
#' \describe{
#'   \item{n}{<integer>; Number of reads per row.}
#'   \item{freq}{<logical>; Is it a frequency matrix or a count matrix.}
#'   \item{offset}{<integer>; Offset}
#'   \item{insertions}{<integer>; Insertions}
#' }
#' @return A \code{consmat} object.
#' @export
#' @examples
#' ###
consmat <- function(x, ...) UseMethod("consmat")

# Internal constructor
Consmat_ <- function(x, n, freq, offset = 0L, insertions = NULL) {
  structure(
    x, n = n, freq = freq, offset = offset, insertions = insertions,
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
  x <- xtabs(formula = count ~ pos + nucleotide, data = x, drop.unused.levels = drop.unused.levels)
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
        sweep(x, 1, n(x), `/`), n = n(x), freq = freq, offset = offset(x),
        insertions = ins(x)
      )
    }
  } else if (!freq) {
    ## recalibrate n
    n(x) <- .rowSums(x, NROW(x), NCOL(x))
    if (is.freq(x)) {
      Consmat_(
        sweep(x, 1, n(x), `*`), n = n(x), freq = freq, offset = offset(x),
        insertions = ins(x)
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
      rs, n, freq = is.freq(x), offset = offset(x), insertions = ins(x)
    )
  } else rs
}

#' @export
`[<-.consmat` <- function(x, i, j, ..., value) {
  rs <- NextMethod()
  ## recalibrate n
  n <- .rowSums(rs, NROW(rs), NCOL(rs))
  Consmat_(
    rs, n, freq = is.freq(x), offset = offset(x), insertions = ins(x)
  )
}

#' @export
as.matrix.consmat <- function(x, ...) {
  rs <- unclass(x)
  attr(rs, "n") <- NULL
  attr(rs, "insertions") <- NULL
  attr(rs, "freq") <- NULL
  attr(rs, "offset") <- NULL
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

## Note for me: flattens the matrix; compare rowsum to rowsums upstream of pos;
## if > t set all to 0
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

#' Create position weight matrix from multiple sequence alignment
#'
#' @param msa A \code{DNAStringSet} object of aligned sequences.
#' @details
#' \code{PWM}: a \code{matrix} with positions row names and nucleotides as
#' column manes. Values are nucleotide weights at a position
#' A ConsensusMatrix is calculated from the MSA using \code{Biostrings::consensusMatrix} and values are converted to probabilities.
#' Pseudocounts are added and values are divided by DNA probabilities and log2 score is reported
#' @return A \code{PWM} matrix.
#'
#' @export
#' @examples
#' ###
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

# mat <- pileup$consmat
#removeError = TRUE
.distributeGaps <- function(mat, bamfile, reference, removeError = FALSE) {
  seq <- .mat2rle(mat)
  if (removeError) {
    gapError <- .getGapErrorBackground(mat, n = 5)
  }
  for (i in which(seq$length > 5)) {
    ## Assign new gap numbers to each position starting from left
    # meanCoverage <- mean(rowSums(mat[seqStart:seqEnd,1:4]))
    seqStart <- sum(seq$lengths[1:(i - 1)]) + 1
    seqEnd <- seqStart + seq$lengths[i] - 1
    region <- paste0(names(reference), ":", seqStart, "-", seqEnd)

    msa <- msa_from_bam(bamfile, reference, paddingLetter = "+", region = region)
    ## Use only sequences spanning the complete region! Every other sequence
    ## gives no Info
    msa <- msa[sapply(msa, function(x) !"+" %in% Biostrings::uniqueLetters(x))]

    if (removeError) {
      aFreq <- Biostrings::alphabetFrequency(msa[1])
      nt <- colnames(aFreq)[which.max(aFreq)]
      nucs <- sum(Biostrings::nchar(msa))
      falseGaps <- round(nucs * gapError)
      while (falseGaps > 0) {
        seqsWithGap <- which(sapply(msa, function(x)
          "-" %in% Biostrings::uniqueLetters(x)))
        if (length(seqsWithGap) > 1) {
          seqsToChange <- sample(seqsWithGap,
                                 min(falseGaps,
                                     length(seqsWithGap)))
        } else {
          seqsToChange <- seqsWithGap
        }
        if (length(seqsToChange) == 0)
          break
        msaChanged <- Biostrings::DNAStringSet(
          sapply(msa[seqsToChange], function(x) {
            # x <- msa[seqsToChange][[1]]
            pos <- Biostrings::matchPattern("-", x)[1]
            Biostrings::replaceLetterAt(x,
                                     Biostrings::start(Biostrings::matchPattern("-", x)[1]),
                                     nt)
          }))
        msa[seqsToChange] <- msaChanged
        falseGaps <- falseGaps - length(seqsToChange)
      }
    }

    msa <- Biostrings::DNAStringSet(sapply(msa, function(x) {
      gapless <- gsub(pattern = "-", "", x)
      nGaps <- Biostrings::nchar(x) - Biostrings::nchar(gapless)
      paste0(paste0(rep("-", nGaps), collapse = ""), gapless)
    }))
    mat[seqStart:seqEnd,] <- t(Biostrings::consensusMatrix(msa))[,VALID_DNA(include = "indel")]
  }
  mat
}

# .distributeGapsConsMat <- function(mat, bamfile, reference, removeError = FALSE){
#   stopifnot(is(mat,"consmat"))
#   if (! "-" %in% colnames(mat)) {
#     flog.warn("No Gaps to distribute")
#     return(mat)
#   }
#   seq <- .mat2rle(mat)
#
#   if (removeError){
#     for (i in which(seq$length > 1)) {
#       seqStart <- sum(seq$lengths[1:(i-1)])+1
#       seqEnd <- seqStart+seq$lengths[i]-1
#       mat[seqStart:seqEnd,"-"] <- sum(mat[seqStart:seqEnd,"-"])/seq$lengths[i]
#     }
#     backgroundGapError <- .getGapErrorBackground(mat)
#     newGapError <- apply(mat, 1, function(x) sum(x)*backgroundGapError)#
#     mat[,"-"] <- round(mat[,"-"] - newGapError)
#     mat[which(mat[,"-"] < 0),"-"] <- 0
#   }
#   mat <- .moveGapsToLeft(mat, bamfile, reference)
#   mat
# }

.getGapErrorBackground <- function(mat, n = 5){
  seq <- .mat2rle(mat)
  consecutiveRegions <- rle((seq$lengths > n))
  longestRegion <- which.max(consecutiveRegions$lengths)
  startSeq <- ifelse(longestRegion == 1,
                     1,
                     sum(consecutiveRegions$lengths[1:(longestRegion - 1)]))
  endSeq <- startSeq + consecutiveRegions$lengths[longestRegion] - 1
  ## Add an offset of 10 to acknowledge bad quality after homopolymer regions
  startSeq <- startSeq + 10
  background <- mat[startSeq:endSeq,]
  backgroundSums <- colSums(background)
  backgroundSums["-"] / sum(backgroundSums)
}
# .moveGapsToLeft <- function(mat){
#   seq <- .mat2rle(mat)
#   ## debug
#   al
#   i <- 7084
#   ##
#   for (i in which(seq$length > 9)) {
#     ## Assign new gap numbers to each position starting from left
#     # meanCoverage <- mean(rowSums(mat[seqStart:seqEnd,1:4]))
#     ## !! not true now: start is now -1, end +1
#     seqStart <- sum(seq$lengths[1:(i-1)])+1
#     seqEnd <- seqStart+seq$lengths[i]-1
#     length <- seq$length[i]
#
#     msa <-
#
#     gaps <- sum(mat[seqStart:seqEnd, "-"])
#     prevCoverage <- sum(mat[(seqStart-1),1:5])
#     baseSums <- colSums(mat[seqStart:seqEnd,1:4])
#     idx <- 0
#     while(gaps > 0 && !seq$values[i] == "-"){
#       # move the gaps to the left
#       gapsAtPos <- min(gaps, prevCoverage)
#       mat[(seqStart+idx),"-"] <- gapsAtPos
#       # if more gaps than bases move the bases from a column to the right
#       if (gaps >= prevCoverage){
#         length <- length - 1
#         # move the bases to the right; distribute them
#         # set bases at gap column to 0
#         mat[(seqStart+idx),1:4] <- 0
#         ## distribute the bases to all other columns
#         newMeanBases <- baseSums/length
#         for (pos in (seqStart+idx+1):seqEnd){
#           mat[pos,1:4] <- newMeanBases
#         }
#       }
#       if (length == 1)
#         break
#       gaps <- gaps - gapsAtPos
#       idx <- idx+1
#     }
#     newSeqStart <- seqStart + idx
#     if (!newSeqStart > seqEnd)
#       mat[newSeqStart:seqEnd, "-"] <- 0
#   }
#   m <- round(mat)
# }
