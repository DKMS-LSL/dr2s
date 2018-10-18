

# Consensus sequences -----------------------------------------------------


#' Construct a consensus sequence
#'
#' @param x \code{pileup} object.
#' @param name Name for consensus sequence.
#' @param type One of "prob", "ambig" or "simple".
#' @param threshold If \code{type == "ambig"}, threshold to call an ambiguous
#' consensus call.
#' @param excludeGaps Exclude gaps at insertion position from consensus
#' calling.
#' @param gapSuppressionRatio The ratio of base/gap above which gaps at
#' insertion position are excluded from from consensus calling.
#' @param forceExcludeGaps Exclude gaps at from consensus calling irrespective
#' of the frequency of the gap.
#' @param ... Additional arguments.
#' @return A \code{\linkS4class{BStringSet}} object with a metadata list
#' containing the slots:
#' \describe{
#'   \item{zscore}{}
#'   \item{freq}{}
#'   \item{ambigs}{}
#'   \item{insertions}{}
#'   \item{deletions}{}
#'   \item{consmat}{}
#' }
#' @export
#' @examples
#' ###
conseq <- function(x,
                   name                = "conseq",
                   type                = c("prob", "ambig", "simple"),
                   threshold           = NULL,
                   excludeGaps         = TRUE,
                   gapSuppressionRatio = 2/5,
                   forceExcludeGaps    = FALSE,
                   ...)
  UseMethod("conseq")
#' @export
conseq.pileup <- function(x,
                          name                = "conseq",
                          type                = c("prob", "ambig", "simple"),
                          threshold           = NULL,
                          excludeGaps         = TRUE,
                          gapSuppressionRatio = 2/5,
                          forceExcludeGaps    = FALSE, ...) {
  if (is.null(threshold)) {
    threshold <- x$threshold
  }
  x <- consmat(x, freq = FALSE)
  conseq(x, name = name, type = type, threshold = threshold,
         excludeGaps = excludeGaps,
         gapSuppressionRatio = gapSuppressionRatio,
         forceExcludeGaps = forceExcludeGaps, ... )
}
#' @export
conseq.matrix <- function(x,
                          name                = NULL,
                          type                = c("prob", "ambig", "simple"),
                          threshold           = NULL,
                          excludeGaps         = TRUE,
                          gapSuppressionRatio = 2/5,
                          forceExcludeGaps    = FALSE, ...) {
  type <- match.arg(type, c("prob", "ambig", "simple"))
  if (type == "ambig" && is.null(threshold)) {
    stop("Must set threshold for ambiguity consensus calling!")
  }
  x <- consmat(x, freq = FALSE)
  conseq <- switch(type,
    prob  = .makeProbConsensus_(x, excludeGaps = excludeGaps,
                                forceExcludeGaps = forceExcludeGaps,
                                gapSuppressionRatio = gapSuppressionRatio),
    ambig = .makeAmbigConsensus_(x, threshold, excludeGaps = excludeGaps),
    simple = .makeSimpleConsensus_(x)

  )
  names(conseq) <- name
  conseq
}

## strict consensus based on z-scores
## <excludeGaps> affects behaviour at insertion positions
## if <excludeGaps> the gap count at insertion position is set to zero if
## the ratio of base/gap >= gapSuppressionRatio (2/3), which allows calling
## the alternate base even if at lower frequency than the gap.
## if <forceExcludeGaps> all gap counts will be set to zero.
.makeProbConsensus_ <- function(x,
                                excludeGaps = TRUE,
                                forceExcludeGaps = FALSE,
                                gapSuppressionRatio = 2/5,
                                asString = FALSE) {




  #cseq <- conseq(reads, "hap" %<<% hptype, "prob", excludeGaps = TRUE)

  xOri <- x
  if (excludeGaps && length(ins_ <- as.character(ins(x))) > 0) {
    x <- .suppressGaps_(x, ins = ins_,
                        gapSuppressionRatio = gapSuppressionRatio)
  }
  if (forceExcludeGaps) {
    x[, "-"] <- 0

  }
  # don't allow gaps at beginning and end
  maxbases <- names(unlist(unname(apply(x, 1, function(a)
    list(which(a == max(a))[1])))))
  maxbase  <- which(maxbases != "-")
  maxgap   <- which(maxbases == "-")
  # remove insertions for the consensus
  x[, "+"] <- 0

  if (!length(maxgap) == 0) {
    if (min(maxgap) < min(maxbase)) {
      excludeFromStart <- min(
        which(maxbases == "-")):(min(which(maxbases != "-")) - 1)
      x[excludeFromStart,"-"] <- 0
    }
    if (max(maxgap) > max(maxbase)) {
      excludeFromEnd <- (max(
        which(maxbases != "-")) + 1):max(which(maxbases == "-"))
      x[excludeFromEnd,"-"] <- 0
    }
  }

  ## remove rows which have been set to zero
  if (length(i <- which(n(x) == 0L)) > 0) {
    x <- x[-i, ]
  }
  rowsd <- mean(n(x)) / 2
  z <- sweep(sweep(x, 1, .rowMeans(x, NROW(x), NCOL(x)), `-`), 1, rowsd, `/`)
  bases <- colnames(x)[apply(z, 1, which.max)]
  if (asString) {
    return(paste0(bases, collapse = ""))
  }
  dels <- bases == "-"
  seq  <- Biostrings::BStringSet(paste0(bases[!dels], collapse = ""))

  # fix zscore; rm del positions
  z <- z[!dels,]

  metadata(seq) <- list(
    zscore     = unname(apply(z, 1, max)),
    freq       = NULL,
    ambigs     = NULL,
    insertions = ins(xOri),
    deletions  = unname(which(dels)),
    consmat    = xOri
  )
  return(seq)
}

## consensus with ambiguities
## <excludeGaps> affects behaviour at polymorphic positions:
## if <!excludeGaps> a gap ambiguity will be called (small letter)
## if <excludeGaps> the alternate base(s) will be called irrespective
## of the frequency of the gap
.makeAmbigConsensus_ <- function(x,
                                  threshold,
                                  excludeGaps = FALSE,
                                  asString = FALSE) {
  ## Filter all bases with a frequency > threshold
  cmf <- consmat(x, freq = TRUE)
  ## remove rows which have been set to zero
  if (length(i <- which(n(cmf) == 0L)) > 0) {
    cmf <- cmf[-i, ]
  }
  baselist <- apply(cmf, 1, function(m) {
    rs <- m[i <- m > threshold]
    names(rs) <- names(m)[i]
    list(rs)
  })
  s <- lapply(baselist, function(b) {
    b <- unlist(b)
    if (length(b) == 1) {
      list(base = names(b), freq = unname(b))
    }
    else if (length(b) > 1) {
      ## sort by name
      b <- b[order(names(b))]
      ## if we have a gap and excludeGaps == TRUE
      if (any(gap <- names(b) == "-") && excludeGaps) {
        b <- b[!gap]
      }
      NUC <- paste0(names(b), collapse = "")
      list(
        base = names(CODE_MAP())[charmatch(NUC, CODE_MAP())] %|na|% "N",
        freq = b
      )
    }
    else {
      stop("No bases?")
    }
  })
  bases <- vapply(s, `[[`, "base", FUN.VALUE = "")
  if (asString) {
    return(paste0(bases, collapse = ""))
  }
  ambigs <- x[which(!bases %in% DNA_BASES()), ]
  attr(ambigs, "ambiguities") <- unname(bases[which(!bases %in% DNA_BASES())])
  dels <- bases == "-"
  seq  <- Biostrings::BStringSet(paste0(bases[!dels], collapse = ""))
  metadata(seq) <- list(
    zscore     = NULL,
    freq       = unname(vapply(s, function(x) sum(x[["freq"]]), 0)),
    ambigs     = ambigs,
    insertions = ins(x),
    deletions  = unname(which(dels)),
    consmat    = x
  )
  seq
}
.makeSimpleConsensus_ <- function(x) {
  cmf <- consmat(x, freq = TRUE)
  if (any(na_ <- is.na(rowSums(cmf)))) {
    cmf[na_, ] <- 0
    cmf <- cbind(cmf, N = ifelse(na_, 1, 0))
  }
  paste0(colnames(cmf)[apply(cmf, 1, which.max)], collapse = "")
}


## Helper functions
# gapSuppressionRatio = base/gap
.suppressGaps_ <- function(x, ins, gapSuppressionRatio = 2/5) {
  x0 <- x[dimnames(x)$pos %in% ins, ]
  ## if the ratio of the most freqent base to gap is greater/equal to
  ## gapSuppressionRatio set the gap count to zero (i.e. suppress the gap)
  i <- which(apply(x0[, c("A", "C", "G", "T")], 1, max) /
               x0[, "-"] >= gapSuppressionRatio)
  if (length(j <- dimnames(x0)$pos[i]) > 0) {
    x[j, "-"] <- 0
  }

  x
}


.getWriteConseq <- function(pileup, name, type, threshold, forceExcludeGaps,
                           conseqPath) {
  conseq <- conseq(pileup$consmat, name = name, type = type,
                   threshold = threshold, forceExcludeGaps = TRUE)
  Biostrings::writeXStringSet(
    Biostrings::DNAStringSet(gsub("[-+]", "N", conseq)),
    conseqPath)
  conseq
}
