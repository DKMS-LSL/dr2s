

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
  assert_that(all(VALID_DNA("indel") %in% colnames(x)))
  n <- .rowSums(x[, VALID_DNA("del")], NROW(x), 5L)
  x <- if (freq) sweep(x, 1, n, `/`) else x
  Consmat_(x, n, freq, ...)
}

#' @export
consmat.tbl_df <- function(x, freq = TRUE, drop.unused.levels = FALSE, ...) {
  assert_that(all(c("pos", "nucleotide", "count") %in% colnames(x)))
  ## Add empty positions to pileup
  # gaps <- which(!1:max(x$pos) %in% x$pos)
  # x <- dplyr::bind_rows(x, lapply(gaps, function(i, seqname) {
  #   tibble::tibble(seqnames = seqname, pos = i, nucleotide = "-", count = 0L)
  # }, seqname = x$seqnames[[1]])) %>%
  #   dplyr::arrange(pos)
  # x$nucleotide <- factor(x$nucleotide)

  ctab <- stats::xtabs(formula = count ~ pos + nucleotide, data = x,
                       drop.unused.levels = drop.unused.levels)
  cmat <- matrix(ctab, NROW(ctab), NCOL(ctab))
  dimnames(cmat) <- dimnames(ctab)
  ## restrict matrix cols to GATC-+ and reorder accordingly
  cmat <- cmat[, VALID_DNA("indel")]
  ## calc number of reads at a position (only using GATC-)
  n <- .rowSums(cmat[, VALID_DNA("del")], NROW(cmat), 5L)
  cmat <- if (freq) sweep(cmat, 1, n, `/`) else cmat
  Consmat_(cmat, n, freq, ...)
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
  }
  else if (!freq) {
    if (is.freq(x)) {
      Consmat_(
        sweep(x, 1, n(x), `*`), n = n(x), freq = freq,
        offsetBases = offsetBases(x), insertions = ins(x)
      )
    } else {
      ## recalibrate n
      n(x) <- .rowSums(x[, VALID_DNA("del")], NROW(x), 5L)
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
    n <- .rowSums(rs[, VALID_DNA("del")], NROW(rs), 5L)
    Consmat_(
      rs, n, freq = is.freq(x), offsetBases = offsetBases(x),
      insertions = ins(x))
  } else rs
}

#' @export
`[<-.consmat` <- function(x, i, j, ..., value) {
  rs <- NextMethod()
  ## recalibrate n
  n <- .rowSums(rs[, VALID_DNA("del")], NROW(rs), 5L)
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
  cmat <- cmat[VALID_DNA("del"), , drop = FALSE]
  cmat <- sweep(cmat, 2, colSums(cmat), "/")
  ## Add pseudocount
  cmat <- cmat + 1/length(msa)
  ## Normalise
  cmat <- scale(cmat, center = FALSE,
                 scale = colSums(cmat))
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
  structure(cmat, class = c("pwm", "matrix"))
}

# mat <- consmat(pileup)
#removeError = TRUE
.distributeGaps <- function(mat, bamfile, removeError = FALSE, ...) {
  indent <- list(...)$indent %||% indentation()
  maxHomopolymerLength <- list(...)$maxHomopolymerLength %||% 5
  seq <- .mat2rle(mat)
  if (removeError) {
    gapError <- .getGapErrorBackground(mat, n = maxHomopolymerLength)
    flog.info("%sEstimate indel noise <%0.3g> to suppress spurious gaps", indent(),
              gapError, name = "info")
  }
  idx <- which(seq$length > maxHomopolymerLength)
  idx <- idx[!seq$values[idx] == "-"]
  workers <- min(length(idx), .getIdleCores())
  bpparam <- BiocParallel::MulticoreParam(workers = workers, log = FALSE)
  flog.info("%sUse %s workers to rightshift gaps at homopolymer positions <%s>", indent(),
            workers, comma(idx), name = "info")
  ## NOTE: Parallelising reading from a bam file fails catastrophically on
  ## the cluster (but not locally) with:
  ##   [E::bgzf_uncompress] Inflate operation failed: invalid distance too far back
  ##   [E::bgzf_read] Read block operation failed with error -1 after 0 of 4 bytes
  ## Also, need to reopen bam file every time before calling stackStringsFromBam()
  ## on a range. Otherwise returns empty DNAStringSets after first range.
  bam <- Rsamtools::BamFile(bamfile)
  refname  <- GenomeInfoDb::seqnames(GenomeInfoDb::seqinfo(bam))
  seqStart <- vapply(idx - 1, function(i) sum(seq$lengths[seq_len(i - 1)]) + 1, FUN.VALUE = numeric(1))
  seqEnd <- seqStart + seq$lengths[idx] - 1
  GRNGS  <- S4Vectors::split(GenomicRanges::GRanges(paste0(refname, ":", seqStart, "-", seqEnd)), seq_along(seqStart))
  MSAs   <- lapply(GRNGS, function(grange) {
    Rsamtools::open.BamFile(bam)
    GenomicAlignments::stackStringsFromBam(
      bam, param = grange, Lpadding.letter = "+", Rpadding.letter = "+", use.names = TRUE)
  })
  Rsamtools::close.BamFile(bam)
  ## Collect changed matrix elements
  changeMat <- BiocParallel::bpmapply(function(msa, start, end) {
    ## Skip position if empty
    if (length(msa) == 0) {
      NULL
    } else {
      ## Use only reads spanning the complete region; every other read
      ## gives no info
      msa <- msa[vapply(msa, function(x) !"+" %in% Biostrings::uniqueLetters(x), FUN.VALUE = logical(1))]
      ## Check again if still sequence there
      if (length(msa) == 0)
        return(NULL)
      if (removeError) {
        aFreq <- colSums(Biostrings::letterFrequency(msa, VALID_DNA("none")))
        nt <- names(aFreq)[which.max(aFreq)]
        nucs <- length(msa)*unique(Biostrings::width(msa))
        falseGaps <- floor(nucs*gapError)
        while (falseGaps > 0) {
          seqsWithGap <- which(Biostrings::letterFrequency(msa, "-") > 0)
          seqsToUpdate <- if (length(seqsWithGap) > 1) {
            set.seed(3)
            sample(seqsWithGap, min(falseGaps, length(seqsWithGap)))
          } else {
            seqsWithGap
          }
          if (length(seqsToUpdate) == 0)
            break
          at  <- matrix(FALSE, nrow = length(msa), ncol = Biostrings::width(msa)[1])
          vm  <- Biostrings::vmatchPattern("-", msa[seqsToUpdate])
          vmi <- vapply(Biostrings::startIndex(vm), `[`, 1, FUN.VALUE = integer(1))
          replace_at <- matrix(nrow = length(vmi), c(seqsToUpdate, vmi))
          at[replace_at] <- TRUE
          nts <- Biostrings::DNAStringSet(ifelse(rowSums(at) == 1, nt, ""))
          msa <- Biostrings::replaceLetterAt(msa, at, nts)
          falseGaps <- falseGaps - length(seqsToUpdate)
        }
      }

      gapI <- Biostrings::startIndex(Biostrings::vmatchPattern("-", msa, fixed = TRUE))
      g0   <- which(!vapply(gapI, function(g) {
        is.null(g) || (g[1L] == 1 && all(diff(g) == 1))
      }, FUN.VALUE = logical(1)))
      msa[g0] <- Biostrings::DNAStringSet(vapply(msa[g0], function(x) {
        gapless <- gsub("-", "", x, fixed = TRUE)
        nGaps <- Biostrings::nchar(x) - Biostrings::nchar(gapless)
        paste0(rep("-", nGaps), collapse = "") %<<% gapless
      }, FUN.VALUE = character(1)))

      list(
        seqStart = start,
        seqEnd = end,
        matPart = t(Biostrings::consensusMatrix(msa, as.prob = FALSE))[, VALID_DNA(include = "indel")]
      )
    }
  }, msa = MSAs, start = seqStart, end = seqEnd, SIMPLIFY = FALSE, BPPARAM = bpparam)
  ##
  ## remove empty entries resulting from too bad coverage
  nullMats <- vapply(changeMat, is.null, FUN.VALUE = logical(1))
  changeMat <- changeMat[!nullMats]
  ## Change the matrix
  for (i in changeMat) {
    mat[i$seqStart:i$seqEnd, ] <- i$mat
  }
  mat
}

.getGapErrorBackground <- function(mat, n = 5) {
  seq <- .mat2rle(mat)
  consecutiveRegions <- rle(seq$lengths > n)
  longestRegion <- which.max(consecutiveRegions$lengths)
  startSeq <- ifelse(longestRegion == 1,
                     1,
                     sum(consecutiveRegions$lengths[seq_len(longestRegion - 1)]))
  endSeq <- startSeq + consecutiveRegions$lengths[longestRegion] - 1
  ## Add an offset of 10 to acknowledge bad quality after homopolymer regions
  startSeq <- startSeq + 10
  background <- mat[startSeq:endSeq, VALID_DNA("del")]
  backgroundSums <- colSums(background)
  backgroundSums["-"]/sum(backgroundSums)
}

.extractIdsFromMat <- function(mat, readIds) {
  consmat(t(
      Biostrings::consensusMatrix(mat[readIds], as.prob = FALSE)[
        VALID_DNA(include = "indel"), ]
    ), freq = FALSE)
}
