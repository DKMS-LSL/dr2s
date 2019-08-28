

# Consensus sequences -----------------------------------------------------


#' Construct a consensus sequence
#'
#' @param x \code{pileup} object.
#' @param name Name for consensus sequence.
#' @param type One of "prob", "ambig" or "simple".
#' @param threshold If \code{type = "ambig"}, threshold to call an ambiguous
#' consensus call.
#' @param suppressAllGaps If \code{type = "prob"} or \code{type = "ambig"},
#' suppress all gaps from consensus calling irrespective of the frequency of the gap.
#' @param suppressInsGaps If \code{type = "prob"}, suppress gaps at insertion position
#' from consensus calling based on \code{columnOccupancy}.
#' @param columnOccupancy Minimum occupancy (1 - fraction of gap) below which
#' bases at insertion position are excluded from from consensus calling.
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
                   name            = "conseq",
                   type            = c("prob", "ambig", "simple"),
                   threshold       = NULL,
                   suppressAllGaps = FALSE,
                   suppressInsGaps = TRUE,
                   gapThreshold    = NULL,
                   columnOccupancy = 0.4,
                   ...)
  UseMethod("conseq")

#' @export
conseq.pileup <- function(x,
                          name            = "conseq",
                          type            = c("prob", "ambig", "simple"),
                          threshold       = NULL,
                          suppressAllGaps = FALSE,
                          suppressInsGaps = TRUE,
                          gapThreshold    = NULL,
                          columnOccupancy = 0.4,
                          ...) {
  if (is.null(threshold)) {
    threshold <- x$threshold
  }
  x <- consmat(x, freq = FALSE)
  conseq(x, name = name, type = type, threshold = threshold,
         suppressAllGaps = suppressAllGaps, suppressInsGaps = suppressInsGaps,
         gapThreshold = gapThreshold, columnOccupancy = columnOccupancy, ...)
}

#' @export
conseq.matrix <- function(x,
                          name            = NULL,
                          type            = c("prob", "ambig", "simple"),
                          threshold       = NULL,
                          suppressAllGaps = FALSE,
                          suppressInsGaps = TRUE,
                          gapThreshold    = NULL,
                          columnOccupancy = 0.4,
                          ...) {
  type <- match.arg(type, c("prob", "ambig", "simple"))
  if (type == "ambig" && is.null(threshold)) {
    stop("Must set threshold for ambiguity consensus calling!")
  }
  x <- consmat(x, freq = FALSE)
  conseq <- switch(type,
    prob   = .makeProbConsensus_(x,
                                 suppressAllGaps = suppressAllGaps,
                                 suppressInsGaps = suppressInsGaps,
                                 gapThreshold    = gapThreshold,
                                 columnOccupancy = columnOccupancy,
                                 ...),
    ambig  = .makeAmbigConsensus_(x, threshold, suppressAllGaps,
                                  gapThreshold   = gapThreshold,
                                  ...),
    simple = .makeSimpleConsensus_(x)
  )
  names(conseq) <- name
  conseq
}

## strict consensus based on z-scores
## if <suppressAllGaps> all gap counts are set to zero.
## if <suppressInsGaps> only gap counts at insertion position are affected
## and if the maximum base frequency >= <columnOccupancy> at an insertion position
## the gap count is set to zero.
## if gapThreshold is set, all gaps below 1-threshold are set to zero.
.makeProbConsensus_ <- function(x,
                                suppressAllGaps = FALSE,
                                gapThreshold = NULL,
                                suppressInsGaps = TRUE,
                                columnOccupancy = 0.4,
                                asRle = FALSE,
                                ...) {
  xOri <- x
  ## always suppress insertions for consensus call!
  x[, "+"] <- 0
  ## remove rows which have been set to zero
  if (length(i <- which(n(x) == 0L)) > 0) {
    x <- x[-i, ]
  }
  ## suppress all deletions
  if (suppressAllGaps) {
    x[, "-"] <- 0
  }  else { ## or suppress deletions specifically at insertion positions or below threshold
    if (suppressInsGaps && length(ins_ <- as.character(ins(x))) > 0) {
      x <- .suppressGaps_(x, ins = ins_, columnOccupancy = columnOccupancy)
    }
    if (!is.null(gapThreshold)) {
      x[(x[, "-"] / rowSums(x) < (1 - gapThreshold)),"-"] <- 0
      ## Use a lower threshold for gaps inside homopolymers
      seqrle <- rle(VALID_DNA(include = "none")[apply(x[,VALID_DNA(include = "none")], 1, which.max)])
      n <- which(seqrle$lengths > 3)
      end <- cumsum(seqrle$lengths)
      start <- c(0, end)
      hps <- unlist(lapply(n, function(pos) (start[pos]+1):end[pos]))
      hpGaps <- hps[x[hps, "-"] / rowSums(x[hps,]) < (1 - 0.67 * gapThreshold)]
      x[hpGaps,"-"] <- 0
    }
  }
  ## never allow gaps at the beginning and end of a sequence
  if (!suppressAllGaps) {
    maxbases <- colnames(x)[apply(x, 1, which.max)]
    maxbase  <- which(maxbases != "-")
    maxgap   <- which(maxbases == "-")
    if (length(maxgap) > 0) {
      if (min(maxgap) < min(maxbase)) {
        excludeFromStart <- min(
          which(maxbases == "-")):(min(which(maxbases != "-")) - 1)
        x[excludeFromStart, "-"] <- 0
      }
      if (max(maxgap) > max(maxbase)) {
        excludeFromEnd <- (max(
          which(maxbases != "-")) + 1):max(which(maxbases == "-"))
        x[excludeFromEnd, "-"] <- 0
      }
    }
  }
  rowsd <- mean(.rowSums(x, NROW(x), NCOL(x)))/2 ## WHY?
  z <- sweep(sweep(x, 1, .rowMeans(x, NROW(x), NCOL(x)), `-`), 1, rowsd, `/`)
  bases <- colnames(x)[apply(z, 1, function(i) 
    as.numeric(ifelse(max(i) == 0, 5, which.max(i))))]

  if (asRle) {
    return(rle(bases))
  }

  dels <- bases == "-"
  seq  <- Biostrings::BStringSet(paste0(bases[!dels], collapse = ""))
  # fix zscore; rm del positions
  z <- z[!dels, ]

  S4Vectors::metadata(seq) <- list(
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
## <suppressAllGaps> affects behaviour at polymorphic positions:
## if <!suppressAllGaps> a gap ambiguity may be called (small letter)
## if <suppressAllGaps> the alternate base will be called irrespective
## of the frequency of the gap
.makeAmbigConsensus_ <- function(x,
                                 threshold,
                                 suppressAllGaps = FALSE,
                                 suppressInsGaps = TRUE,
                                 gapThreshold = NULL,
                                 columnOccupancy = 0.4,
                                 asString = FALSE,
                                 ...) {
  ## always suppress insertions for consensus call!
  x[, "+"] <- 0
  ## remove rows which have been set to zero
  if (length(i <- which(n(x) == 0L)) > 0) {
    x <- x[-i, ]
  }
  ## suppress all deletions
  if (suppressAllGaps) {
    x[, "-"] <- 0
  } else { ## or suppress deletions specifically at insertion positions or below threshold
    if (suppressInsGaps && length(ins_ <- as.character(ins(x))) > 0) {
      x <- .suppressGaps_(x, ins = ins_, columnOccupancy = columnOccupancy)
    }
    if (!is.null(gapThreshold)) {
      x[(x[, "-"] / rowSums(x) < (1 - gapThreshold)),"-"] <- 0
      ## Use a lower threshold for gaps inside homopolymers
      seqrle <- rle(VALID_DNA(include = "none")[apply(x[,VALID_DNA(include = "none")], 1, which.max)])
      n <- which(seqrle$lengths > 3)
      end <- cumsum(seqrle$lengths)
      start <- c(0, end)
      hps <- unlist(lapply(n, function(pos) (start[pos]+1):end[pos]))
      hpGaps <- hps[x[hps, "-"] / rowSums(x[hps,]) < (1 - 0.67 * gapThreshold)]
      x[hpGaps,"-"] <- 0
    }
  }
  ## never allow gaps at the beginning and end of a sequence
  if (!suppressAllGaps) {
    maxbases <- colnames(x)[apply(x, 1, which.max)]
    maxbase  <- which(maxbases != "-")
    maxgap   <- which(maxbases == "-")
    if (length(maxgap) > 0) {
      if (min(maxgap) < min(maxbase)) {
        excludeFromStart <- min(
          which(maxbases == "-")):(min(which(maxbases != "-")) - 1)
        x[excludeFromStart, "-"] <- 0
      }
      if (max(maxgap) > max(maxbase)) {
        excludeFromEnd <- (max(
          which(maxbases != "-")) + 1):max(which(maxbases == "-"))
        x[excludeFromEnd, "-"] <- 0
      }
    }
  }
  # ## suppress all deletions
  # if (suppressAllGaps) {
  #   x[, "-"] <- 0
  # }
  # ## Filter all bases with a frequency > threshold
  cmf <- consmat(x, freq = TRUE)
  # 
  # ## Take a higher threshold for gaps. Everything with more than 1.5 * threshold
  # ## is probably not present
  # #cmf[11,]
  # #x[11,]
  # 
  # cmf[cmf[,"-"] > 1 - threshold, VALID_DNA("none")] <- 0
  #     
  baselist <- apply(cmf, 1, function(m) {
    rs <- m[i <- m > threshold]
    names(rs) <- names(m)[i]
    list(rs)
  })
  # baselist[469]
  # x[469,]
  # cmf[469,]
  
  s <- lapply(baselist, function(b) {
    b <- unlist(b)
    if (length(b) == 0) {
      ## this arises if all bases are below threshold
      ## we've seen this with promethION data:
      # Consensus Matrix: 3 x 6
      # nucleotide
      # pos            G          A T         C         -          +
      #   3055 0.1166667 0.28333333 0 0.2500000 0.1888889 0.16111111
      list(base = "N", freq = b)
    }
    else if (length(b) == 1) {
      list(base = names(b), freq = unname(b))
    }
    else if (length(b) > 1) {
      ## sort by name
      NUC <- paste0(names(b[order(names(b))]), collapse = "")
      list(
        base = names(CODE_MAP())[charmatch(NUC, CODE_MAP())] %|na|% "N",
        freq = b
      )
    }
  })
  bases <- vapply(s, `[[`, "base", FUN.VALUE = "")
  if (asString) {
    return(paste0(bases, collapse = ""))
  }
  ambigs <- x[which(!bases %in% DNA_BASES()), ]
  attr(ambigs, "ambiguities") <- unname(bases[which(!bases %in% DNA_BASES())])
  dels <- unlist(gregexpr("[a-z]", bases)) > 0
  seq  <- Biostrings::BStringSet(paste0(bases[!dels], collapse = ""))

  S4Vectors::metadata(seq) <- list(
    zscore     = NULL,
    freq       = unname(vapply(s, function(x) sum(x[["freq"]]), 0)),
    ambigs     = ambigs,
    insertions = ins(x),
    deletions  = unname(which(dels)),
    consmat    = x
  )

  return(seq)
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
.suppressGaps_ <- function(x, ins, columnOccupancy = 0.4) {
  x0 <- x[dimnames(x)$pos %in% ins, ]
  ## if frequency of the most frequent base is greater or equal to the
  ## columnOccupancy (1 - gap frequency) set the gap count to zero
  ## (i.e. suppress the gap)
  i <- which(apply(x0[, c("A", "C", "G", "T")], 1, max) /
               x0[, "-"] >= columnOccupancy)
  if (length(j <- dimnames(x0)$pos[i]) > 0) {
    x[j, "-"] <- 0
  }
  x
}

.writeConseq <- function(x, name, type, threshold = NULL, suppressAllGaps = TRUE,
                         replaceIndel = "", conspath, gapThreshold = NULL, ...) {
  conseq0 <- conseq(consmat(x, freq = FALSE), name = name, type = type,
                   threshold = threshold, suppressAllGaps = suppressAllGaps, gapThreshold = gapThreshold, ...)
  conseq1 <- Biostrings::DNAStringSet(stripIndel(conseq0, replace = replaceIndel))
  Biostrings::writeXStringSet(conseq1, conspath, append = FALSE, compress = FALSE)
  return(invisible(conseq1))
}
