

# Pileup ------------------------------------------------------------------


#' Calculate pile-up for a BAM file.
#'
#' @param bamfile BAM file path.
#' @param threshold Threshold used for SNP calling.
#' @param ... Additional parameters passed to 
#' \code{\link[Rsamtools]{PileupParam}}.
#' @inheritParams Rsamtools::PileupParam
#'
#' @details
#' Returns a \code{pileup} object:
#' A \code{list} with slots:
#' \describe{
#'   \item{bamfile}{<character>; Path to the bam file used to construct the 
#'   pileup}
#'   \item{threshold}{<numeric>; Threshold used for calling SNPs}
#'   \item{param}{A \code{\link[Rsamtools]{PileupParam}} object}
#'   \item{pilup}{A \code{data.frame} with colums:
#'     \describe{
#'       \item{seqnames}{<character>; The RNAME field}
#'       \item{pos}{<integer>; Genomic position of base}
#'       \item{nucleotide}{<character>; Base}
#'       \item{count}{<integer>; Count of base}
#'     }
#'   }
#'   \item{consmat}{Consensus matrix}
#' }
#'
#' @return A \code{pileup} object. See \code{Details}.
#' @export
#' @examples
#' ###
Pileup <- function(bamfile,
                   threshold = 0.20,
                   ...) {
  param <- list(...)
  
  if (is.null(param$distinguish_strands))
      param$distinguish_strands <- FALSE
      
  bfl <- BamFile(bamfile)
  pParam <- do.call(PileupParam, c(param))
  
  sParam <- ScanBamParam(
    flag = scanBamFlag(
      isUnmappedQuery = FALSE,
      isSecondaryAlignment = FALSE
    )
  )
  ## TODO: check if still correct
  # pileup <- Rsamtools:::.pileup(file = bfl, scanBamParam = sParam, 
  #                               pileupParam = pParam)
  pileup <- pileup(file = bfl, scanBamParam = sParam, 
                                pileupParam = pParam)
  pileup <-
    dplyr::mutate(dplyr::tbl_df(pileup),
                  seqnames = strsplitN(as.character(.data$seqnames), "~", 1, 
                                       fixed = TRUE))
  structure(
    list(
      bamfile   = bamfile,
      threshold = threshold,
      param     = pParam,
      pileup    = pileup,
      consmat   = consmat(pileup, freq = FALSE)
    ),
    class = c("pileup", "list")
  )
}

# Methods: pileup ---------------------------------------------------------


#' @export
print.pileup <- function(x, asString = FALSE, ...) {
  values <- sapply(methods::slotNames(x$param), methods::slot, object = x$param)
  info <- paste(methods::slotNames(x$param), values, sep=": ", collapse="; ")
  msg <- if (asString) "" else sprintf("An object of class '%s'.\n", 
                                       class(x)[1])
  msg <- sprintf("%s Bamfile: %s\n Threshold: %s\n %s\n",
                 msg, basename(x$bamfile), x$threshold,
                 paste0(strwrap(info, initial = "Params: ", exdent = 4), 
                        collapse = "\n"))
  if (asString)
    return(msg)
  else {
    cat(msg)
    print(x$consmat, n = 4)
  }
}


# Helpers -----------------------------------------------------------------


.pileupFindInsertionPositions_ <- function(x, threshold) {
  cm  <- consmat(x, freq = TRUE)
  pos <- .ambiguousPositions(cm, threshold = threshold)
  pos[cm[pos, "+"] > threshold]
}

.pileupGetInsertions_ <- function(x, threshold) {
  res <- list()
  colnm <- colnames(x$consmat)
  inpos <- .pileupFindInsertionPositions_(x, threshold)
  inpos <- inpos[!inpos %in% 1:5]
  inpos <- inpos[!inpos %in% (NROW(x$consmat) - 5):NROW(x$consmat)]
  bamfile <- x$bamfile
  if (length(inpos) > 0) {
    inseqs <- .getInsertions(bamfile, inpos)
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
  res <- Filter(length,res)
  res
}

.pileupIncludeInsertions <- function(x, threshold = NULL) {
  stopifnot(is(x, "pileup"))
  if (!"+" %in% colnames(x$consmat)) {
    flog.warn("No insertions to call!", name = "info")
    return(x)
  }
  if (is.null(threshold)) {
    threshold <- x$threshold
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

.msaFromBam <- function(bamfile, refseq = NULL, paddingLetter = "+", 
                        region = NULL) {
  if (is.null(region)) {
    nRefseq <- names(refseq)
    lRefseq <- length(refseq[[1]])
    region <- paste0(nRefseq, ":1-", lRefseq)
  }
  assert_that(grepl(pattern = "^[[:alnum:]_\\*#\\.]+:\\d+-\\d+", region))
  GenomicAlignments::stackStringsFromBam(
    bamfile,
    param = region,
    Lpadding.letter = paddingLetter,
    Rpadding.letter = paddingLetter,
    use.names = TRUE)
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
  
  
  ggplot(dtbg, aes(x = pos, y = count)) +
    geom_bar(stat = "identity", position = position_identity(), 
             fill = "grey80") +
    geom_bar(data = dtpoly, stat = "identity", position = position_stack(),
             aes(fill = nucleotide), width = dtpoly$width) +
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

.getInsertions <- function(bamfile, inpos, reference, reads = NULL) {
  stopifnot(is.numeric(inpos))
  inpos <- sort(inpos)
  flog.info("  Extracting insertions at position %s ...", comma(inpos), 
            name = "info")

  ## Get the actual position of the first insertion character, not the last
  ##  matching position, so +1
  inpos <- inpos + 1

  ## Get the reference
  reference   <- seqinfo(BamFile(bamfile))@seqnames[1]
  inposRanges <- GenomicRanges::GRanges(reference, 
                                        IRanges::IRanges(start = inpos, 
                                                         end = inpos))
  bamParam    <- ScanBamParam(what = "seq", which = inposRanges)
  bam <- GenomicAlignments::readGAlignments(bamfile, param = bamParam, 
                                            use.names = TRUE)
  ## Use only reads of interest if specified
  if (!is.null(reads))
    bam <- bam[names(bam) %in% reads]
  ## Use only reads with an insertion
  bamI <- bam[sapply(GenomicAlignments::cigar(bam), function(x) grepl("I",  x))]
  insSeqs <- foreach(i = inpos) %do% {
    message("Position: ", i)
    bamPos <- bamI[GenomicAlignments::start(bamI) <= i & 
                     GenomicAlignments::end(bamI) >= 1]
    insSeq <- lapply(bamPos, function(a) .extractInsertion(a, i))
    Biostrings::DNAStringSet(insSeq)
  }

  ## decrement to last matching position again to work as expected with 
  ## downstream
  names(insSeqs) <- inpos - 1
  insSeqs
}

.extractInsertion <- function(read, i) {
  cigar <- read@cigar

  seq <- Biostrings::DNAString("-")
  ## map insertion position to reference space, get insertions from query space
  insertionQ <- GenomicAlignments::cigarRangesAlongQuerySpace(
    cigar, ops = "I")[[1]]
  insertionR <- GenomicAlignments::cigarRangesAlongReferenceSpace(
    cigar, ops = "I")[[1]]
  insertPos <- insertionQ[which((
    GenomicAlignments::start(read) + insertionR@start - 1) == i)]
  if (length(insertPos) > 0)
    seq <- unlist(Biostrings::subseq(read@elementMetadata$seq, 
                                     start = IRanges::start(insertPos), 
                                     end = IRanges::end(insertPos)))
  seq
}

