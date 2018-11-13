
# Constructor -------------------------------------------------------------

#' Calculate pileup for a BAM file.
#'
#' @param bamfile <\code{character}>; path to bamfile.
#' @param reffile <\code{character}>; Path to reference.
#' @param readtype <\code{character}>; one of "illumina", "pacbio", or "nanopore".
#' @param ... Additional parameters passed to \code{\link[Rsamtools]{PileupParam}}.
#' @inheritParams Rsamtools::PileupParam
#' @param pParam A \code{\link[Rsamtools]{PileupParam}} object.
#'
#' @details
#' Returns a \code{pileup} object:
#' A \code{list} with slots:
#' \describe{
#'   \item{bamfile}{<character>; Path to the bamfile used to construct the
#'   pileup}
#'   \item{reffile}{<character>; Path to the reference used to construct the
#'   bamfile}
#'   \item{refname}{<character>; Name of the reference}
#'   \item{readtype}{<character>; One of "illumina", "pacbio", or "nanopore"}
#'   \item{param}{A \code{\link[Rsamtools]{PileupParam}} object}
#'   \item{pilup}{A \code{tbl_df} with colums:
#'     \describe{
#'       \item{seqnames}{<character>; The RNAME field}
#'       \item{pos}{<integer>; Genomic position of base}
#'       \item{nucleotide}{<character>; Base}
#'       \item{count}{<integer>; Count of base}
#'     }
#'   }
#'   \item{consmat}{<consmat>; The consensus matrix}
#'   \item{reads}{<character>; (optional) names of n top-scoring reads}
#'   \item{stats}{<list>; Run statistics}
#' }
#'
#' @return A \code{pileup} object. See \code{Details}.
#' @export
#' @examples
#' ###
pileup <- function(bamfile,
                   reffile,
                   readtype,
                   ...,
                   pParam) {
  bampath <- normalizePath(bamfile, mustWork = TRUE)
  refpath <- normalizePath(reffile, mustWork = TRUE)
  indent <- list(...)$indent %||% indentation()
  if (missing(pParam)) {
    param <- list(...)
    if (is.null(param$distinguish_strands))
      param$distinguish_strands <- FALSE
    pParam <- do.call(Rsamtools::PileupParam, param)
  }
  assert_that(is(pParam, "PileupParam"))

  bam <- Rsamtools::BamFile(bampath)
  Rsamtools::open.BamFile(bam)
  on.exit(Rsamtools::close.BamFile(bam))

  if (Rsamtools::idxstatsBam(bam)$mapped == 0) {
    flog.info("%sNo reads map to reference", indent(), name = "info")
    stop("No reads map to reference")
  }

  sParam <- Rsamtools::ScanBamParam(
    flag = Rsamtools::scanBamFlag(
      isUnmappedQuery = FALSE,
      isSecondaryAlignment = FALSE
    ),
    which = IRanges::IRangesList()
  )
  pileup <- tibble::as_tibble(Rsamtools::pileup(
    file = bam,
    scanBamParam = sParam,
    pileupParam = pParam))
  refname <- GenomeInfoDb::seqnames(Rsamtools::seqinfo(bam))
  Pileup_(refpath = refpath, bampath = bampath, param = pParam,
          pileup = pileup, refname = refname, readtype = readtype)
}

# Helpers -----------------------------------------------------------------

.pileupFindInsertionPositions_ <- function(x, threshold = 0.20) {
  cm  <- consmat(x, freq = TRUE)
  pos <- .ambiguousPositions(cm, threshold = threshold)
  pos[cm[pos, "+"] > threshold]
}

.pileupGetInsertions_ <- function(x, threshold = 0.20, ...) {
  indent <- list(...)$indent %||% indentation()
  res   <- list()
  colnm <- colnames(consmat(x))
  inpos <- .pileupFindInsertionPositions_(x, threshold)
  inpos <- inpos[!inpos %in% 1:5]
  inpos <- inpos[!inpos %in% (NROW(consmat(x)) - 5):NROW(consmat(x))]
  if (length(inpos) > 0) {
    inseqs <- .getInsertions(bamfile = bampath(x), inpos = inpos, readtype = readtype(x), indent = indent)
    inseqs <- inseqs[order(as.integer(names(inseqs)))]
    for (inseq in inseqs) {
      if (length(inseq) < 100) {
        inseq <- inseq[order(Biostrings::width(inseq))]
        maxWidth <- max(Biostrings::width(inseq))
        inseq0 <- inseq[inseq != "-"]
        if (length(inseq0) > 2) {
          gT <- data.frame(seq_along(inseq0))
          inseq1 <- DECIPHER::AdjustAlignment(
            DECIPHER::AlignSeqs(inseq0, verbose = FALSE,
                                iterations = 0, refinements = 0,
                                restrict = c(-500, 2, 10), normPower = 0)
          )
          inseq1 <- c(inseq1, Biostrings::DNAStringSet(
            rep(paste0(rep.int("-", maxWidth), collapse = ""),
                sum(inseq == "-"))
          ))
        } else {
          inseq1 <- inseq
        }
      }  else {
        # get biostringset which fills all positions with gaps where it is
        # not complete, i.e. width is < max width
        dels <- Biostrings::DNAStringSet(vapply(max(Biostrings::width(inseq)) -
                                                  Biostrings::width(inseq),
                                                function(x) {
          paste0(rep.int("-", x), collapse = "")
        }, FUN.VALUE = ""))
        # merge both biostrings so each seq is of same width
        inseq1 <- Biostrings::xscat(inseq, dels)
      }
      cm <- t(Biostrings::consensusMatrix(inseq1))[, colnm, drop = FALSE]
      cmf <- sweep(cm, 1, rowSums(cm), "/")
      cmf[, "-"] <- 0
      cm <- cm[apply(cmf, 1, function(row) any(row > threshold)), ,
               drop = FALSE]
      res <- c( res, list(cm) )
    }
    res
    names(res) <- names(inseqs)
  }
  # Remove all empty positions for now!! ToDo: Better apply this to the initial
  # ins pos calling in the python script
  res <- Filter(length, res)
  res
}

.pileupIncludeInsertions <- function(x, threshold = 0.2, ...) {
  assert_that(is(x, "pileup"))
  indent <- list(...)$indent %||% indentation()
  if (!"+" %in% colnames(consmat(x))) {
    flog.warn("%sNo insertions to call", indent(), name = "info")
    return(x)
  }
  ins_ <- .pileupGetInsertions_(x, threshold, indent = indent)
  if (length(ins_) == 0) {
    flog.info("%sNo insertions found at threshold <%s>", indent(), threshold, name = "info")
    return(x)
  }
  offsetBases <- 0L
  insIdx <- integer()
  insRun <- integer()
  cm <- consmat(x, freq = FALSE)
  cmAttr <- attributes(cm)
  # cm[680:700,]
  # i <- 1
  for (i in seq_along(ins_)) {
    j <- as.integer(names(ins_[i])) + offsetBases
    cm[j, "+"] <- 0L
    cm <- rbind(cm[1:j, ], ins_[[i]], cm[(j + 1):NROW(cm), ])
    ins_[1]
    insLen <- NROW(ins_[[i]])
    insIdx <- c(insIdx, (j + 1):(j + insLen))
    insRun <- c(insRun, rep(i, insLen))
    offsetBases <- offsetBases + insLen
  }
  attr(insIdx, "run") <- insRun
  cmAttr$dim <- dim(cm)
  cmAttr$dimnames$pos <- as.character(seq_len(NROW(cm)))
  cmAttr$n <- .rowSums(cm, NROW(cm), NCOL(cm))
  cmAttr$insertions <- insIdx
  attributes(cm) <- cmAttr
  x$consmat <- cm
  x
}

.msaFromBam <- function(bamfile, range = NULL, paddingLetter = "+") {
  assert_that(
    is(bamfile, "BamFile"),
    is.null(range) || (is.numeric(range) && length(range) == 2 && range[1] <= range[2])
  )
  if (!Rsamtools::isOpen(bamfile)) {
    Rsamtools::open.BamFile(bamfile)
    on.exit(Rsamtools::close.BamFile(bamfile))
  }
  baminfo <- Rsamtools::seqinfo(bamfile)
  param <- if (is.null(range)) {
    GenomicRanges::GRanges(
      seqnames = GenomeInfoDb::seqnames(baminfo),
      ranges = IRanges::IRanges(start = 1L, end = GenomeInfoDb::seqlengths(baminfo)))
  } else {
    GenomicRanges::GRanges(
      seqnames = GenomeInfoDb::seqnames(baminfo),
      ranges = IRanges::IRanges(start = range[1], end = range[2]))
  }
  GenomicAlignments::stackStringsFromBam(
    file = bamfile, index = bamfile, param = param,
    Lpadding.letter = paddingLetter, Rpadding.letter = paddingLetter,
    use.names = TRUE)
}

.topXReads <- function(bamfile, topx = "auto", ...) {
  indent <- list(...)$indent %||% indentation()
  msa <- .msaFromBam(Rsamtools::BamFile(bamfile))
  ## there is no point in sampling if
  if (is.numeric(topx) && length(msa) <= topx) {
    flog.info("%sNothing to extract from %s mapped longreads", indent(), length(msa), name = "info")
    return(names(msa))
  }
  pwm <- createPWM(msa)
  bpparam <- BiocParallel::MulticoreParam(workers = .getIdleCores())
  res <- do.call(dplyr::bind_rows, suppressWarnings(BiocParallel::bplapply(seq_along(msa),
    function(i, msa, pwm) {
      tibble::tibble(read = names(msa[i]), score = .PWMscore(msa[[i]], pwm))
    }, msa = msa, pwm = pwm, BPPARAM = bpparam)))

  # res <- do.call(dplyr::bind_rows, suppressWarnings(BiocParallel::bplapply(seq_along(msa),
  #   function(s, msa, pwm) {
  #     ## TODO: Make this CPP
  #     read     <- strsplit1(as.character(msa[[s]]), split = "")
  #     readname <- names(msa[s])
  #     b <- sum(vapply(seq_along(read), function(i, pwm, read) pwm[read[i], i],
  #                     pwm = pwm, read = read, FUN.VALUE = numeric(1)))/length(read)
  #     tibble::tibble(read = readname, score = b)
  #   }, msa = msa, pwm = pwm, BPPARAM = bpparam)))

  reads <- if (topx == "auto") {
    .pickTopReads(res, f = 3/5)
  } else {
    dplyr::top_n(res, topx, score)$read
  }

  flog.info("%sFrom %s mapped longreads extract the %0.3g%% (%s) top-scoring reads",
            indent(), length(msa), 100*length(reads)/length(msa), length(reads),
            name = "info")

  reads
}

.PWMscore <- function(read, pwm) {
  ## expect read to be a DNAString or DNAStringSet instance
  ## this is 25x faster than the vapply loop
  ## nxm
  n <- NROW(pwm)
  m <- NCOL(pwm)
  read <- strsplit(as.character(read), split = "")[[1L]]
  i <- match(read, rownames(pwm)) + ((0:(m - 1))*n)
  sum(pwm[i])/m
}

## TODO: think about this more carefully
.pickTopReads <- function(x, f = 1) {
  score_rng <- range(x$score)
  scorecut <- seq(score_rng[1], score_rng[2], length.out = 100)
  readsum <- vapply(scorecut, function(cutoff)
    sum(x$score >= cutoff), FUN.VALUE = double(1))
  optcoef <- scorecut[which.max((readsum^f)*scorecut)]
  # df <- tibble::tibble(
  #   x = seq_along(res$score),
  #   score = sort(res$score),
  #   pick = sort(res$score) >= optcoef
  # )
  # p <- ggplot(df, aes(x, score, colour = pick)) +
  #   geom_point()
  #sum(res$score >= optcoef)
  dplyr::filter(x, score >= optcoef)$read
}


# Summarise and plot ------------------------------------------------------


#' Plot pileup coverage
#'
#' @param x A \code{pileup} object.
#' @param threshold At which frequency do we visualise polymorphisms.
#' @param range Which part along the x-axis do we plot.
#' @param thin Plot only a fraction of the (non-polymorphic) positions.
#' @param width Width of polymorphic positions.
#' @param label Optional plot label.
#' @param drop.indels Don't plot indels.
#'
#' @return A \code{ggplot} object.
#' @export
#' @examples
#' ###
plotPileupCoverage <- function(x, threshold = 0.2, range = NULL, thin = 0.1,
                               width = 1, label = "", drop.indels = FALSE) {
  stopifnot(is(x, "pileup"))
  x <- as.data.table(x$pileup)
  if (drop.indels) {
    x <- x[nucleotide != "-" & nucleotide != "+"]
  }
  if (!is.null(range)) {
    x <- x[pos >= range[1] & pos <= range[2]]
  }
  x[, nucleotide := as.character(nucleotide)]
  x[, freq := count/sum(count), by = pos]
  x[nucleotide == "-" & nucleotide == "+", npoly := sum(freq >= max(
    c(threshold, 0.3))), by = pos]
  x[nucleotide != "-" & nucleotide != "+", npoly := sum(freq >= max(
    c(threshold, 0.2))), by = pos]
  #x
  x[npoly == 1 | (npoly > 1 & freq < threshold), nucleotide := " ", by = pos]
  nonpoly <- x[npoly == 1, unique(pos)]
  nonpolyThin <- nonpoly[seq(1, length(nonpoly),
                             floor(length(nonpoly)/(length(nonpoly)*thin)))]
  dt <- x[, list(count = sum(count)), by = list(pos, nucleotide)]
  dtbg <- dt[pos %in% nonpolyThin]
  dtpoly <- dt[!pos %in% nonpoly][order(pos, nucleotide)]
  dtpoly[, nucleotide := factor(nucleotide,
                                levels = c("+", "-", "A", "C", "G", "T", " "),
                                labels = c("+", "-", "A", "C", "G", "T", " "),
                                ordered = TRUE)]
  setkeyv(dtpoly, c("pos", "nucleotide"))

  ## Set the width for neighbouring SNPs to 1
  positions <- dtpoly$pos
  closePositions <- sort(unique(c(which((positions - 1) %in% positions),
                                  which((positions + 1) %in% positions))))
  dtpoly$width <- width
  dtpoly[closePositions]$width <- 1

  ## suppress warning because of potentially overlapping x intervals
  p <- ggplot(dtbg, aes(x = pos, y = count)) +
    geom_bar(stat = "identity", position = position_stack(), fill = "grey80") +
    geom_bar(aes(fill = nucleotide), data = dtpoly, stat = "identity",
             position = position_stack(), width = dtpoly$width) +
    scale_fill_manual(values = NUCCOL(),
                      limits = c("A", "C", "G", "T", "-", "+")) +
    guides(fill = guide_legend(reverse = TRUE, title = "Bases")) +
    labs(x = "Position [bp]", y = "Count", fill = "Nucleotide", title = label) +
    theme_bw(base_size = 12) +
    theme(
      panel.border = element_blank(),
      panel.grid = element_blank(),
      axis.line = element_line(colour = "grey50")
    )
  p
}

#' Plot basecall frequency
#'
#' @param x A \code{pileup} object.
#' @param threshold At which frequency do we visualise polymorphisms.
#' @param label Optional plot label.
#' @param drop.indels Don't plot indels.
#'
#' @return A \code{ggplot} object.
#' @export
#' @examples
#' ###
plotPileupBasecallFrequency <- function(x, threshold = 0.20, label = "",
                                        drop.indels = FALSE) {
  cm <- consmat(x, freq = FALSE)
  if (drop.indels && "-" %in% colnames(cm)) {
    cm[, "-"] <- 0
  }
  if (drop.indels && "+" %in% colnames(cm)) {
    cm[, "+"] <- 0
  }
  cmlong <- dplyr::filter(as.data.frame(consmat(cm, freq = TRUE)), freq > 0)
  ## TOOD check alpha change
  ggplot(cmlong, aes(x =~ pos, y =~ freq, colour =~ nucleotide)) +
    geom_point(aes(alpha =~ ifelse(freq > threshold, 0.4, 0.2)), size = 0.75,
               shape = 15) +
    scale_color_manual(values = NUCCOL()) +
    scale_alpha_continuous(guide = FALSE) +
    guides(colour = guide_legend(title = "Bases")) +
    scale_y_log10() +
    geom_hline(yintercept = threshold, colour = "grey20", size = 0.25) +
    geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey20",
               size = 0.75) +
    labs(x = "Position [bp]", y = "Base frequency", title = label) +
    theme_bw(base_size = 12) +
    theme(
      panel.border = element_blank(),
      axis.line = element_line(colour = "grey50")
    )
}

.getInsertions <- function(bamfile, inpos, readtype = "illumina", reads = NULL, ...) {
  assert_that(is.numeric(inpos))
  inpos <- sort(inpos)
  indent <- list(...)$indent %||% indentation()
  flog.info("%sCall insertions at positions <%s>", indent(), comma(inpos), name = "info")
  ## Get the actual position of the first insertion character, not the last
  ##  matching position, so +1
  inpos   <- inpos + 1

  bam <- Rsamtools::BamFile(bamfile)
  Rsamtools::open.BamFile(bam)
  on.exit(Rsamtools::close.BamFile(bam))

  baminfo <- GenomicRanges::seqinfo(bam)
  inposRanges <- if (readtype == "illumina") {
    GenomicRanges::GRanges(
      seqnames = GenomeInfoDb::seqnames(baminfo),
      ranges = IRanges::IRanges(start = inpos, end = inpos))
  } else {
    GenomicRanges::GRanges(
      seqnames = GenomeInfoDb::seqnames(baminfo),
      ranges =  IRanges::IRanges(start = 1, end = GenomeInfoDb::seqlengths(baminfo)))
  }
  bamParam <- Rsamtools::ScanBamParam(what = "seq", which = inposRanges)
  aln <- GenomicAlignments::readGAlignments(bam, param = bamParam, use.names = TRUE)

  ## Use only reads of interest if specified
  if (!is.null(reads))
    aln <- aln[names(aln) %in% reads]
  ## Use only reads with an insertion
  readsWithIns <- vapply(GenomicAlignments::cigar(aln),
                         function(x) grepl("I",  x),
                         FUN.VALUE = logical(1))
  alnI <- aln[readsWithIns]
  ## Get all insertions from all reads
  ## TODO: bplapply throws warning in serialize(data, node$con, xdr = FALSE)
  ## 'package:stats' may not be available when loading
  ## For the time being we suppress this warning
  bpparam <- BiocParallel::MulticoreParam(workers = .getIdleCores())
  insSeq <- unlist(Biostrings::DNAStringSetList(suppressWarnings(
    BiocParallel::bplapply(seq_along(alnI), function(a, alnI, inpos) {
      .extractInsertion(alnI[a], inpos)
    }, alnI = alnI, inpos = inpos, BPPARAM = bpparam))))
  ## Extract per positions
  insSeqs <- foreach(i = inpos) %do% {
    #message("Position: ", i)
    insSeq[which(names(insSeq) == i)]
  }
  ## decrement to last matching position again to work as expected with
  ## downstream
  names(insSeqs) <- inpos - 1
  insSeqs <- insSeqs[vapply(insSeqs, function(x)
    length(unlist(x)) > 0, FUN.VALUE = logical(1))]
  insSeqs
}

.extractInsertion <- function(read, inpos) {
  cigar <- read@cigar
  ## map insertion position to reference space, get insertions from query space
  insertionQ <- GenomicAlignments::cigarRangesAlongQuerySpace(
    cigar, ops = "I")[[1]]
  insertionR <- GenomicAlignments::cigarRangesAlongReferenceSpace(
    cigar, ops = "I")[[1]]

  qPos <- GenomicAlignments::start(read) + insertionR@start - 1
  foundPositions <- which(qPos %in% inpos)
  fInPos <- qPos[foundPositions]
  insertPos <- insertionQ[foundPositions]

  if (length(insertPos) > 0) {
    seq <- unlist(Biostrings::DNAStringSetList(lapply(seq_along(insertPos),
                                                      function(ip) {
                                                        ipp <- insertPos[ip]
                                                        Biostrings::subseq(read@elementMetadata$seq,
                                                                           start = IRanges::start(ipp),
                                                                           end = IRanges::end(ipp))
                                                      })))
    names(seq) <- fInPos
    return(seq)
  }
  Biostrings::DNAStringSet("-")
}

.checkCoverage <- function(pileup, forceMapping, plotFile, maptag, ...) {
  indent <- list(...)$indent %||% indentation()
  flog.warn("%sShortreads seem corrupted or the reference is bad!", indent(), name = "info")
  maxCov <- max(rowSums(consmat(pileup, freq = FALSE)))
  q75Cov <- quantile(rowSums(consmat(pileup, freq = FALSE)), 0.75)
  flog.warn("%sMaximum of coverage %s / 75%% quantile %s: %s > 5." %<<%
            " No equal distribution of coverage!" %<<%
            " Have a look at the mapInit plot",
            indent(), maxCov, q75Cov, maxCov/q75Cov, name = "info")
  plt <- plotPileupCoverage(
    x = pileup,
    thin = 0.25,
    width = 2,
    label = maptag,
    drop.indels = TRUE
  )
  flog.error("%sAborting. If you want to force processing set" %<<%
             " forceMapping = TRUE in DR2S object initialisation", indent(), name = "info")
  save_plot(plotFile, plt, base_aspect_ratio = 3,
            onefile = TRUE)
  if (!forceMapping) {
    stop("Shortreads probably of bad quality. Bad coverage distribution.
         Run with forceMapping = TRUE to force processing.")
  } else {
    flog.warn("%sContinue. Be aware that resulsts may not be correct!!", indent(), name = "info")
  }
}
