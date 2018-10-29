
# Class: Pileup -----------------------------------------------------------

## Internal onstructor class "pileup"
Pileup_ <- function(...) {
  dots <- list(...)
  structure(
    list(
      bamfile   = dots$bamfile,  # <character>; path to bamfile.
      reffile   = dots$reffile,  # <character>; path to reference.
      refname   = dots$refname,  # <character>; name of reference.
      readtype  = dots$readtype, # <character>; "illumina", "pacbio", "nanopore".
      param     = dots$param,    # <PileupParam>
      pileup    = dots$pileup,   # <tbl_df> with columns: "seqnames", "pos", "nucleotide", "count".
      consmat   = consmat(dots$pileup, freq = FALSE), # <consmat>
      reads     = NULL,          # <character>; names of n top-scoring reads.
      stats     = list()         # <named list>
    ),
    class = c("pileup", "list")
  )
}

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
  bamfile <- normalizePath(bamfile, mustWork = TRUE)
  reffile <- normalizePath(reffile, mustWork = TRUE)

  if (missing(pParam)) {
    param <- list(...)
    if (is.null(param$distinguish_strands))
      param$distinguish_strands <- FALSE
    pParam <- do.call(Rsamtools::PileupParam, param)
  }
  assert_that(is(pParam, "PileupParam"))

  bam <- Rsamtools::BamFile(bamfile)
  Rsamtools::open.BamFile(bam)
  on.exit(Rsamtools::close.BamFile(bam))
  if (Rsamtools::idxstatsBam(bam)$mapped == 0) {
    flog.info("No reads maps to the reference", name = "info")
    stop("No reads maps to the reference")
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
  Pileup_(bamfile = bamfile, reffile = reffile, refname = refname,
          readtype = readtype, param = pParam, pileup = pileup)
}


# Methods: pileup ---------------------------------------------------------


#' @export
print.pileup <- function(x, asString = FALSE, ...) {
  params <- methods::slotNames(x$param)
  names(params) <- params
  values <- lapply(params, methods::slot, object = x$param)
  info <- paste(methods::slotNames(x$param), values,
                sep = ": ", collapse = "; ")
  msg <- if (asString) "" else
    sprintf("An object of class '%s'.\n", class(x)[1])
  msg <- sprintf("%s Bamfile: %s\n Reference: %s\n Readtype: %s\n %s\n",
                 msg, basename(path(x)), basename(refpath(x)), readtype(x),
                 paste0(strwrap(info, initial = "Params: ", exdent = 4),
                        collapse = "\n"))
  if (asString)
    return(msg)
  else {
    cat(msg)
    print(consmat(x), n = 4)
  }
}

#' @export
path.pileup <- function(x, ...) {
  x$bamfile
}

#' @export
refpath.pileup <- function(x, ...) {
  x$reffile
}

#' @export
refname.pileup <- function(x, ...) {
  x$refname
}

#' @export
readtype.pileup <- function(x, ...) {
  x$readtype
}

#' @export
stats.pileup <- function(x, ...) {
  x$stats
}

#' @export
consmat.pileup <- function(x, freq = FALSE, ...) {
  cm <- if (is(x$consmat, "consmat")) x$consmat else consmat(x$pileup)
  consmat(cm, freq = freq)
}

#' @export
`consmat<-.pileup` <- function(x, value) {
  assert_that(is(value, "consmat"))
  x$consmat <- value
  x
}

#' @rdname pileup
#' @export
reads <- function(x, ...) UseMethod("reads")
#' @export
reads.pileup <- function(x, ..) {
  x$reads
}

#' @rdname pileup
#' @export
`reads<-` <- function(x, value) UseMethod("reads<-")
#' @export
`reads<-.pileup` <- function(x, value) {
  assert_that(is.character(value))
  x$reads <- value
  x
}


# Helpers -----------------------------------------------------------------


.pileupFindInsertionPositions_ <- function(x, threshold = 0.20) {
  cm  <- consmat(x, freq = TRUE)
  pos <- .ambiguousPositions(cm, threshold = threshold)
  pos[cm[pos, "+"] > threshold]
}

.pileupGetInsertions_ <- function(x, threshold = 0.20) {
  res   <- list()
  colnm <- colnames(consmat(x))
  inpos <- .pileupFindInsertionPositions_(x, threshold)
  inpos <- inpos[!inpos %in% 1:5]
  inpos <- inpos[!inpos %in% (NROW(consmat(x)) - 5):NROW(consmat(x))]
  if (length(inpos) > 0) {
    inseqs <- .getInsertions(bamfile = path(x), inpos = inpos, readtype = readtype(x))
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

.pileupIncludeInsertions <- function(x, threshold = 0.2) {
  assert_that(is(x, "pileup"))
  if (!"+" %in% colnames(consmat(x))) {
    flog.warn("No insertions to call!", name = "info")
    return(x)
  }
  ins_ <- .pileupGetInsertions_(x, threshold)
  if (length(ins_) == 0) {
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

.msaFromBam <- function(bamfile, region = NULL, paddingLetter = "+") {
  if (!is(bamfile, "BamFile")) {
    bamfile <- Rsamtools::BamFile(bamfile)
  }
  if (!Rsamtools::isOpen(bamfile)) {
    Rsamtools::open.BamFile(bamfile)
    on.exit(Rsamtools::close.BamFile(bamfile))
  }
  if (is.null(region)) {
    baminfo <- Rsamtools::seqinfo(bamfile)
    myParam <- GenomicRanges::GRanges(
      seqnames = GenomeInfoDb::seqnames(baminfo),
      ranges = IRanges::IRanges(start = 1L, end = GenomeInfoDb::seqlengths(baminfo)))
  } else if (is(region, "GRanges")) {
    myParam <- region
  } else {
    assert_that(grepl(pattern = "^[[:alnum:]_\\*#:\\.\\-]+:\\d+-\\d+", region))
    m <- regexpr("(?<seqname>^[[:alnum:]_\\*#:\\.\\-]+):(?<start>\\d+)-(?<end>\\d+)",
                 region, perl = TRUE)
    starts <- attr(m, "capture.start")
    ends   <- starts + attr(m, "capture.length") - 1
    myParam <- GenomicRanges::GRanges(
      seqnames = substr(region, starts[, "seqname"],  ends[, "seqname"]),
      ranges = IRanges::IRanges(
        start = as.integer(substr(region, starts[, "start"],  ends[, "start"])),
        end =  as.integer(substr(region, starts[, "end"],  ends[, "end"]))))
  }
  GenomicAlignments::stackStringsFromBam(
    bamfile,
    param = myParam,
    Lpadding.letter = paddingLetter,
    Rpadding.letter = paddingLetter,
    use.names = TRUE)
}

.topXReads <- function(bamfile, n = 2000) {
  msa <- .msaFromBam(bamfile)
  mat <- createPWM(msa)
  mat["+",] <- 0
  bpparam <- BiocParallel::MulticoreParam(workers = .getIdleCores())
  res <- do.call(dplyr::bind_rows, BiocParallel::bplapply(seq_along(msa),
    function(s, aln, mat) {
      seq <- as.character(aln[[s]])
      seq <- unlist(strsplit(seq, split = ""))
      read <- names(aln[s])
      ## Make this CPP
      b <- sum(vapply(seq_along(seq), function(x, mat, seq) mat[seq[x], x],
                       mat = mat, seq = seq, FUN.VALUE = numeric(1)))
      t <- tibble::tibble(read, b/length(seq))
      names(t) <- c("read", "score")
      t
  }, aln = msa, mat = mat, BPPARAM = bpparam))
  dplyr::arrange(res, dplyr::desc(score))[1:n, ]$read
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

.getInsertions <- function(bamfile, inpos, readtype = "illumina", reads = NULL) {
  assert_that(is.numeric(inpos))
  inpos <- sort(inpos)
  flog.info("  Extracting insertions at positions %s", comma(inpos), name = "info")
  ## Get the actual position of the first insertion character, not the last
  ##  matching position, so +1
  inpos   <- inpos + 1
  bamfile <- Rsamtools::BamFile(bamfile)
  baminfo <- GenomicRanges::seqinfo(bamfile)
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
  bam <- GenomicAlignments::readGAlignments(bamfile, param = bamParam, use.names = TRUE)

  ## Use only reads of interest if specified
  if (!is.null(reads))
    bam <- bam[names(bam) %in% reads]
  ## Use only reads with an insertion
  readsWithIns <- vapply(GenomicAlignments::cigar(bam),
                         function(x) grepl("I",  x),
                         FUN.VALUE = logical(1))
  bamI <- bam[readsWithIns]
  ## Get all insertions from all reads
  ## TODO: bplapply throws warning in serialize(data, node$con, xdr = FALSE)
  ## 'package:stats' may not be available when loading
  ## For the time being we suppress this warning
  bpparam <- BiocParallel::MulticoreParam(workers = .getIdleCores())
  insSeq <- unlist(Biostrings::DNAStringSetList(suppressWarnings(
    BiocParallel::bplapply(seq_along(bamI), function(a, bamI, inpos) {
      .extractInsertion(bamI[a], inpos)
    }, bamI = bamI, inpos = inpos, BPPARAM = bpparam))))
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

.checkCoverage <- function(pileup, forceMapping, plotFile, maptag) {
  flog.warn(" Shortreads seem corrupted or the reference is bad!",
            name = "info")
  maxCov <- max(rowSums(consmat(pileup, freq = FALSE)))
  q75Cov <- quantile(rowSums(consmat(pileup, freq = FALSE)), 0.75)
  flog.warn("   Maximum of coverage %s / 75%% quantile %s: %s > 5." %<<%
              " No equal distribution of coverage!" %<<%
              " Have a look at the mapInit plot",
            maxCov, q75Cov, maxCov/q75Cov, name = "info")
  plt <- plotPileupCoverage(
    x = pileup,
    thin = 0.25,
    width = 2,
    label = maptag,
    drop.indels = TRUE
  )
  flog.error(paste(" Aborting. If you want to force processing set ",
                   "forceMapping = TRUE in DR2S object initialisation",
                   " "),
             name = "info")
  save_plot(plotFile, plt, base_aspect_ratio = 3,
            onefile = TRUE)
  if (!forceMapping) {
    stop("Shortreads probably of bad quality. Bad coverage distribution.
         Run with forceMapping = TRUE to force processing.")
  } else {
    flog.warn(" Continue. Be aware that resulsts may not be correct!!",
              name = "info")
  }
}
