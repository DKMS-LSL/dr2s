

# Consensus sequences -----------------------------------------------------


#' Construct a consensus sequence
#'
#' @param x \code{pileup} object.
#' @param name Name for consensus sequence.
#' @param type One of "prob" or "ambig".
#' @param threshold If \code{type == "ambig"}, threshold to call an ambiguous
#' consensus call.
#' @param excludeGaps Exclude gaps at insertion position from consensus
#' calling.
#' @param gapSuppressionRatio The ratio of base/gap above which gaps at
#' insertion position are excluded from from consensus calling.
#' @param forceExcludeGaps Exclude gaps at from consensus calling irrespective
#' of the frequency of the gap.
#' @param ... Additional arguments.
#'
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
conseq <- function(x, ...) UseMethod("conseq")

#' @export
conseq.pileup <- function(x,
                          name = "conseq",
                          type = c("prob", "ambig"),
                          threshold = NULL,
                          excludeGaps = TRUE,
                          gapSuppressionRatio = 2/5,
                          forceExcludeGaps = FALSE, ...) {
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
                          name = NULL,
                          type = c("prob", "ambig"),
                          threshold = NULL,
                          excludeGaps = TRUE,
                          gapSuppressionRatio = 2/5,
                          forceExcludeGaps = FALSE, ...) {
  type <- match.arg(type, c("prob", "ambig"))
  if (type == "ambig" && is.null(threshold)) {
    stop("Must set threshold for ambiguity consensus calling!")
  }
  x <- consmat(x, freq = FALSE)
  conseq <- switch(type,
    prob  = .makeProbConsensus_(x, excludeGaps = excludeGaps,
                                 forceExcludeGaps = forceExcludeGaps,
                                 gapSuppressionRatio = gapSuppressionRatio),
    ambig = .makeAmbigConsensus_(x, threshold, excludeGaps = excludeGaps)
  )
  names(conseq) <- name
  conseq
}

#' @keywords internal
#' @export
simpleConsensus <- function(x) {
  cmf <- consmat(x, freq = TRUE)
  if (any(na_ <- is.na(rowSums(cmf)))) {
    cmf[na_, ] <- 0
    cmf <- cbind(cmf, N = ifelse(na_, 1, 0))
  }
  paste0(colnames(cmf)[apply(cmf, 1, which.max)], collapse = "")
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


# Summarise and plot ------------------------------------------------------


#' Plot probability of consensus bases.
#'
#' @param seqA consensus sequence A.
#' @param seqB consensus sequence B.
#' @param textSize Text size.
#' @param pointSize Point size.
#' @param labelA Label A.
#' @param labelB Label B.
#'
#' @return A \code{ggplot} object.
#' @export
#'
#' @examples
#' ##
# threshold = "auto"
# textSize = 3
# pointSize = 1
plotConseqProbability <- function(cseqs,
                                    threshold = "auto",
                                    textSize = 3,
                                    pointSize = 1) {
  labels <- sapply(cseqs, function(x) x$label)
  seqs   <- sapply(cseqs, function(x) x$cseq)
  # get labels and tags
  if (all(nzchar(labels))) {
    tags   <- lapply(labels, function(x) gsub("[<>]", "", 
                                              strsplit1(x, " ", fixed = TRUE)))
    groups <- lapply(tags, function(x) dot(setdiff(x, Reduce(intersect, tags))))
    label  <- dot(Reduce(intersect, tags))
  } else {
    label  <- ""
    groups <- 1:length(cseqs)
  }

  if (all(unlist(lapply(seqs, function(x) !is.null(metadata(x)$zscore))))) {
    ylabel <- "Z-score probability"
  } else {
    ylabel <- "Frequency"
  }

  df <- dplyr::bind_rows(
    foreach(hp = names(cseqs)) %do% {
      seq <- seqs[[hp]]
      dplyr::data_frame(
        group = groups[[hp]],
        pos   = seq_len(Biostrings::width(seq)),
        base  = strsplit1(as.character(seq), split = ""),
        prob  = if (!is.null(metadata(seq)$zscore)) {
          as.numeric(pnorm(metadata(seq)$zscore))
        } else {
          metadata(seq)$freq
        }
      )
    }
  )

  lower <- if (threshold == "auto") {
    dplyr::summarise(dplyr::group_by(df, group), lower = quantile(prob, 0.0025))
  } else {
    dplyr::data_frame(group = unlist(groups), lower = threshold)
  }

  df2 <- dplyr::filter(dplyr::left_join(df, lower, by = "group"), prob <= lower)
  ggplot(df, aes(pos, prob)) +
    facet_grid(group ~ .) +
    geom_point(aes(colour = base), alpha = 0.5, size = pointSize) +
    scale_color_manual(values = NUCCOL()) +
    geom_label(aes(x = pos, y = prob - prob*0.05, label = pos), data = df2,
               colour = "black", size = textSize) +
    geom_hline(aes(yintercept = lower), linetype = "dashed", 
               colour = "gray20", data = lower) +
    ylim(c(0.25, 1)) +
    xlab("Position [bp]") +
    ylab(ylabel) +
    ggtitle(label) +
    theme_bw()
}
