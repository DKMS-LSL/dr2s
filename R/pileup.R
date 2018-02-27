

# Pileup ------------------------------------------------------------------

#' Calculate pile-up for a BAM file.
#'
#' @param bamfile BAM file path.
#' @param threshold Threshold used for SNP calling.
#' @inheritParams Rsamtools::PileupParam
#'
#' @details
#' Returns a \code{pileup} object:
#' A \code{list} with slots:
#' \describe{
#'   \item{bamfile}{<character>; Path to the bam file used to construct the pileup}
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
                   max_depth = 250,
                   ##min_base_quality = 13,
                   min_base_quality = 0,
                   min_mapq = 0,
                   min_nucleotide_depth = 1,
                   min_minor_allele_depth = 0,
                   distinguish_strands = FALSE,
                   distinguish_nucleotides = TRUE,
                   ignore_query_Ns = TRUE,
                   include_deletions = TRUE,
                   include_insertions = FALSE,
                   left_bins = NULL,
                   query_bins = NULL) {
  bfl <- Rsamtools::BamFile(bamfile)
  if (packageVersion("Rsamtools") < "1.24.0") {
    p_param <- Rsamtools::PileupParam(
      max_depth = max_depth,
      min_base_quality = min_base_quality,
      min_mapq = min_mapq,
      min_nucleotide_depth = min_nucleotide_depth,
      min_minor_allele_depth = min_minor_allele_depth,
      distinguish_strands = distinguish_strands,
      distinguish_nucleotides = distinguish_nucleotides,
      ignore_query_Ns = ignore_query_Ns,
      include_deletions = include_deletions,
      include_insertions = include_insertions,
      cycle_bins = NULL %||% numeric()
    )
  } else {
    p_param <- Rsamtools::PileupParam(
      max_depth = max_depth,
      min_base_quality = min_base_quality,
      min_mapq = min_mapq,
      min_nucleotide_depth = min_nucleotide_depth,
      min_minor_allele_depth = min_minor_allele_depth,
      distinguish_strands = distinguish_strands,
      distinguish_nucleotides = distinguish_nucleotides,
      ignore_query_Ns = ignore_query_Ns,
      include_deletions = include_deletions,
      include_insertions = include_insertions,
      left_bins = left_bins,
      query_bins = query_bins
    )
  }
  s_param <- Rsamtools::ScanBamParam(
    flag = Rsamtools::scanBamFlag(
      isUnmappedQuery = FALSE,
      isSecondaryAlignment = FALSE
    )
  )
  pileup <- Rsamtools:::.pileup(file = bfl, scanBamParam = s_param, pileupParam = p_param)
  pileup <-
    dplyr::mutate(dplyr::tbl_df(pileup),
                  seqnames = strsplitN(as.character(seqnames), "~", 1, fixed = TRUE))
  structure(
    list(
      bamfile   = bamfile,
      threshold = threshold,
      param     = p_param,
      pileup    = pileup,
      consmat   = consmat(pileup, freq = FALSE)
    ),
    class = c("pileup", "list")
  )
}


# Methods: pileup ---------------------------------------------------------


#' @export
print.pileup <- function(x, as_string = FALSE, ...) {
  values <- sapply(slotNames(x$param), slot, object = x$param)
  info <- paste(slotNames(x$param), values, sep=": ", collapse="; ")
  msg <- if (as_string) "" else sprintf("An object of class '%s'.\n", class(x)[1])
  msg <- sprintf("%s Bamfile: %s\n Threshold: %s\n %s\n",
                 msg, basename(x$bamfile), x$threshold,
                 paste0(strwrap(info, initial = "Params: ", exdent = 4), collapse = "\n"))
  if (as_string)
    return(msg)
  else {
    cat(msg)
    print(x$consmat, n = 4)
  }
}


# Helpers -----------------------------------------------------------------

pileup_find_insertion_positions_ <- function(x, threshold) {
  cm  <- consmat(x, freq = TRUE)
  pos <- ambiguous_positions(cm, threshold = threshold)
  pos[cm[pos, "+"] > threshold]
}

pileup_get_insertions_ <- function(x, threshold) {
  res <- list()
  colnm <- colnames(x$consmat)
  inpos <- pileup_find_insertion_positions_(x, threshold)
  inpos <- inpos[!inpos %in% 1:5]
  inpos <- inpos[!inpos %in% (NROW(x)-5):NROW(x)]
  bamfile <- x$bamfile
  if (length(inpos) > 0) {
    inseqs <- .getInsertions(bamfile, inpos)
    inseqs <- inseqs[order(as.integer(names(inseqs)))]
    for (inseq in inseqs) {
      if (length(inseq) < 100) {
        inseq <- inseq[order(Biostrings::width(inseq))]
        max_width <- max(Biostrings::width(inseq))
        inseq0 <- inseq[inseq != "-"]
        if (length(inseq0) > 2) {
          gT <- data.frame(seq_along(inseq0))
          inseq1 <- DECIPHER::AdjustAlignment(
            DECIPHER::AlignSeqs(inseq0, verbose = FALSE,
                                iterations = 0, refinements = 0,
                                restrict = c(-500, 2, 10), normPower = 0)
          )
          inseq1 <- c(inseq1, Biostrings::DNAStringSet(
            rep(paste0(rep.int("-", max_width), collapse = ""), sum(inseq == "-"))
          ))
        } else {
          inseq1 <- inseq
        }
      }  else {
        # get biostringset which fills all positions with gaps where it is not complete, i.e. width is < max width
        dels <- Biostrings::DNAStringSet(vapply(max(Biostrings::width(inseq)) - Biostrings::width(inseq), function(x) {
          paste0(rep.int("-", x), collapse = "")
        }, FUN.VALUE = ""))
        # merge both biostrings so each seq is of same width
        inseq1 <- Biostrings::xscat(inseq, dels)
      }
      cm <- t(Biostrings::consensusMatrix(inseq1))[, colnm, drop = FALSE]
      cmf <- sweep(cm, 1, rowSums(cm), "/")
      cmf[, "-"] <- 0
      cm <- cm[apply(cmf, 1, function(row) any(row > threshold)), , drop = FALSE]
      res <- c( res, list(cm) )
    }
    res
    names(res) <- names(inseqs)
  }
  # Remove all empty positions for now!! ToDo: Better apply this to the initial ins pos calling in the python script
  res <- Filter(length,res)
  res
}

#' @keywords internal
#' @export
# debug
# threshold = NULL
# x <- pileup
pileup_include_insertions <- function(x, threshold = NULL) {
  stopifnot(is(x, "pileup"))
  if (! "+" %in% colnames(x$consmat)) {
    flog.warn("No insertions to call!", name = "info")
    return(x)
  }
  if (is.null(threshold)) {
    threshold <- x$threshold
  }
  ins_ <- pileup_get_insertions_(x, threshold)
  if (length(ins_) == 0) {
    return(x)
  }
  offset <- 0L
  ins_idx <- integer()
  ins_run <- integer()
  cm <- consmat(x, freq = FALSE)
  cm_attr <- attributes(cm)
  # i <- 2
  for (i in seq_along(ins_)) {
    j <- as.integer(names(ins_[i])) + offset
    # cm[(j-1):(j+1), ]
    cm[j, "+"] <- 0L
    cm <- rbind(cm[1:j, ], ins_[[i]], cm[(j + 1):NROW(cm), ])
    ins_len <- NROW(ins_[[i]])
    ins_idx <- c(ins_idx, (j + 1):(j + ins_len))
    ins_run <- c(ins_run, rep(i, ins_len))
    offset <- offset + ins_len
  }
  attr(ins_idx, "run") <- ins_run
  cm_attr$dim <- dim(cm)
  cm_attr$dimnames$pos <- as.character(seq_len(NROW(cm)))
  cm_attr$n <- .rowSums(cm, NROW(cm), NCOL(cm))
  cm_attr$insertions <- ins_idx
  attributes(cm) <- cm_attr
  x$consmat <- cm
  x
}
#refseq <- hla.map1$map1$A$ref$refseq
## get msa from bam file
# debug
#bamfile <- self$mapInit$bamfile
#refseq <- self$getRefSeq()
msa_from_bam <- function(bamfile, refseq = NULL, paddingLetter = "+", region = NULL){
  if (is.null(region)) {
    nRefseq <- names(refseq)
    lRefseq <- length(refseq[[1]])
    region <- paste0(nRefseq, ":1-", lRefseq)
  }
  assert_that(grepl(pattern = "^[[:alnum:]_\\*#\\.]+:\\d+-\\d+", region))
  GenomicAlignments::stackStringsFromBam(bamfile, param=region, Lpadding.letter = paddingLetter, Rpadding.letter = paddingLetter, use.names = TRUE)
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
plot_pileup_coverage <- function(x, threshold = 0.2, range = NULL, thin = 0.1, width = 1, label = "", drop.indels = FALSE) {
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
  x[, npoly := sum(freq >= threshold), by = pos]
  x[npoly == 1 | (npoly > 1 & freq < threshold), nucleotide := " ", by = pos]
  nonpoly <- x[npoly == 1, unique(pos)]
  nonpoly_thin <- nonpoly[seq(1, length(nonpoly), floor(length(nonpoly)/(length(nonpoly)*thin)))]
  dt <- x[, list(count = sum(count)), by = list(pos, nucleotide)]
  dtbg <- dt[pos %in% nonpoly_thin]
  dtpoly <- dt[!pos %in% nonpoly][order(pos, nucleotide)]
  dtpoly[, nucleotide := factor(nucleotide, levels = c("+", "-", "A", "C", "G", "T", " "),
                                labels = c("+", "-", "A", "C", "G", "T", " "), ordered = TRUE)]
  setkeyv(dtpoly, c("pos", "nucleotide"))

  ggplot(dtbg, aes(x = pos, y = count)) +
    geom_bar(stat = "identity", position = position_identity(), fill = "grey80") +
    geom_bar(data = dtpoly, stat = "identity", position = position_stack(),
             aes(fill = nucleotide), width = width) +
    scale_fill_manual(values = NUCCOL(), limits = c("A", "C", "G", "T", "-", "+")) +
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
plot_pileup_basecall_frequency <- function(x, threshold = 0.20, label = "", drop.indels = FALSE) {
  cm <- consmat(x, freq = FALSE)
  if (drop.indels && "-" %in% colnames(cm)) {
    cm[, "-"] <- 0
  }
  if (drop.indels && "+" %in% colnames(cm)) {
    cm[, "+"] <- 0
  }
  cmlong <- dplyr::filter(as.data.frame(consmat(cm, freq = TRUE)), freq > 0)
  ggplot(cmlong, aes(x = pos, y = freq, colour = nucleotide)) +
    geom_point(aes(alpha = ifelse(freq > threshold, 0.4, 0.2)), size = 0.75, shape = 15) +
    scale_color_manual(values = NUCCOL()) +
    scale_alpha_continuous(guide = FALSE) +
    guides(colour = guide_legend(title = "Bases")) +
    scale_y_log10() +
    geom_hline(yintercept = threshold, colour = "grey20", size = 0.25) +
    geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey20", size = 0.75) +
    labs(x = "Position [bp]", y = "Base frequency", title = label) +
    theme_bw(base_size = 12) +
    theme(
      panel.border = element_blank(),
      axis.line = element_line(colour = "grey50")
    )
}

# get_reads_at_ppos <- function(bamfile, ppos, nucA = None, nucB = None) {
#
#   ppos
#     if isinstance(ppos, list):
#         ppos = [int(x) for x in ppos]
#     else:
#         ppos = [int(ppos)]
#
#     if nucA and nucB:
#         if isinstance(nucA, list):
#             nucA = [str(x) for x in nucA]
#         else:
#             nucA = [str(nucA)]
#         if isinstance(nucB, list):
#             nucB = [str(x) for x in nucB]
#         else:
#             nucB = [str(nucB)]
#         assert len(nucA) == len(nucB)
#         assert len(ppos) == len(nucA)
#         posdict = {pos: {A if B != "+" else "-": [], B if A != "+" else "-": []} for pos, A, B in zip(ppos, nucA, nucB)}
#     else:
#         posdict = {pos: {"A": [], "C": [], "G": [], "T": [], "-": [], "+": []} for pos in ppos}
#     bam = pysam.AlignmentFile(bamfile, "rb")
#     ref = bam.references[0]
#     pileup = bam.pileup(ref)
#     for pileupCol in pileup:
#         pos = int(pileupCol.reference_pos)
#         if pos not in posdict.keys():
#             pass
#         else:
#             sys.stdout.write("Extracting read IDs for " + " and ".join(posdict[pos].keys()) + " at " + str(pos + 1) + "\n")
#             for read in pileupCol.pileups:
#                 refpos = read.alignment.get_reference_positions(full_length = True)
#                 try:
#                     i = refpos.index(pos)
#                     base = read.alignment.query_sequence[i]
#                     if "+" not in posdict[pos].keys(): # We don't expect an insertion here
#                         try:
#                             posdict[pos][base].append(read.alignment.query_name)
#                         except KeyError:
#                             pass
#                     else: # We do expect an insertion here
#                         j = refpos.index(pos + 1)
#                         ins = read.alignment.query_sequence[i:j][1:]
#                         if ins:
#                             posdict[pos]["+"].append(read.alignment.query_name)
#                         else:
#                             posdict[pos]["-"].append(read.alignment.query_name)
#                 except ValueError:
#                     try:
#                         posdict[pos]["-"].append(read.alignment.query_name)
#                     except KeyError:
#                         pass
#     return dict(posdict)
#
#
#
#
# }
.getInsertions <- function(bamfile, inpos, reference) {
  stopifnot(is.numeric(inpos))
  inpos <- sort(inpos)
  flog.info("Extracting insertions at position %s",
            paste(inpos, collapse = ", "),
            name = "info")

  ## Get the actual position of the first insertion character, not the last
  ##  matching position, so +1
  inpos <- inpos+1

  ## Get the reference
  reference <- Rsamtools::seqinfo(Rsamtools::BamFile(bamfile))@seqnames[1]

  inposRanges <- GenomicRanges::GRanges(reference, IRanges::IRanges(start = inpos, end = inpos))
  bamParam <- Rsamtools::ScanBamParam(what = "seq", which = inposRanges)
  bam <- GenomicAlignments::readGAlignments(bamfile, param = bamParam, use.names = TRUE)
  ## Use only reads with an insertion
  bamI <- bam[sapply(GenomicAlignments::cigar(bam), function(x) grepl("I",  x))]
  insSeqs <- foreach(i = inpos) %do% {
    message("Position: ", i)
    bamPos <- bamI[GenomicAlignments::start(bamI)<=i & GenomicAlignments::end(bamI)>=1]
    insSeq <- lapply(bamPos, function(a) .extractInsertion(a,i))
    Biostrings::DNAStringSet(insSeq)
  }
  names(insSeqs) <- inpos
  insSeqs

}
.extractInsertion <- function(read, i) {
  cigar <- read@cigar
  seq <- Biostrings::DNAString("-")
  insertion <- GenomicAlignments::cigarRangesAlongQuerySpace(cigar, ops = "I")[[1]]
  insertPos <- insertion[which((GenomicAlignments::start(read) + insertion@start - 1) == i)]
  if (length(insertPos) > 0)
    seq <- unlist(Biostrings::subseq(read@elementMetadata$seq, start = IRanges::start(insertPos), end = IRanges::end(insertPos)))
  seq
}
