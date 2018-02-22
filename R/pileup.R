

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

pileup_get_read_ids_ <- function(x,
                                 ppos = NULL,
                                 threshold = NULL,
                                 nucA = NULL,
                                 nucB = NULL,
                                 all_bases = FALSE) {
  stopifnot(is(x, "pileup"))
  if (is.null(threshold)) {
    threshold <- x$threshold
  }
  cm <- consmat(x$pileup, freq = TRUE)
  if (is.null(ppos)) {
    ppos <- cpp_polymorphic_positions(cm, threshold)
  }
  if (length(ppos) == 0) {
    return(NULL)
  }
  if (!all_bases && is.null(nucA) && is.null(nucB)) {
    rs <- cpp_top2_cols(cm[ppos, ])
    nucA <- colnames(cm)[rs$i1]
    nucB <- colnames(cm)[rs$i2]
  }
  bamfile <- x$bamfile
  out <- rPython::python.call("py_get_reads_at_ppos", bamfile = bamfile, ppos = ppos - 1,
                              nucA = nucA, nucB = nucB, simplify = FALSE)
  names(out) <- as.character(as.integer(names(out)) + 1)
  out[order(as.integer(names(out)))]
}

#' @keywords internal
#' @export
pileup_include_read_ids <- function(x, threshold = NULL) {
  rs <- pileup_get_read_ids_(x, threshold = NULL)
  ids(x$consmat) <- rs
  x
}

pileup_find_insertion_positions_ <- function(x, threshold) {
  # polymorphic_positions(consmat(x, freq = TRUE), threshold = threshold) %>%
  #   dplyr::filter(a1 == "+" | a2 == "+") %>%
  #   dplyr::select(position) %>%
  #   unlist() %>%
  #   unname() %>%
  #   as.integer()
  cm  <- consmat(x, freq = TRUE)
  pos <- ambiguous_positions(cm, threshold = threshold)
  pos[cm[pos, "+"] > threshold]
}
# debug
# threshold <- 0.2
pileup_get_insertions_ <- function(x, threshold) {
  res <- list()
  colnm <- colnames(x$consmat)
  inpos <- pileup_find_insertion_positions_(x, threshold)
  inpos <- inpos[!inpos %in% 1:5]
  inpos <- inpos[!inpos %in% (NROW(x)-5):NROW(x)]
  bamfile <- x$bamfile
  if (length(inpos) > 0) {
    inseqs <- rPython::python.call("py_get_insertions", bamfile, inpos, simplify = FALSE)
    inseqs <- inseqs[order(as.integer(names(inseqs)))]
    ## DEBUG ANFANG ##
    # x$consmat[inpos, ]
    #inseq <- inseqs[[3]]
    # table(inseq)
    ## DEBUG ENDE ##
    for (inseq in inseqs) {
      s <- Biostrings::DNAStringSet(inseq)
      if (length(s) < 100) {
        s <- s[order(Biostrings::width(s))]
        max_width <- max(Biostrings::width(s))
        s0 <- s[s != "-"]
        if (length(s0) > 2) {
          gT <- data.frame(seq_along(s0))
          # s1 <- DECIPHER::AdjustAlignment(
          #   DECIPHER::AlignSeqs(s0, guideTree = gT, verbose = FALSE,
          #                       iterations = 0, refinements = 0,
          #                       restrict = c(-500, 2, 10), normPower = 0)
          # )
          s1 <- DECIPHER::AdjustAlignment(
            DECIPHER::AlignSeqs(s0, verbose = FALSE,
                                iterations = 0, refinements = 0,
                                restrict = c(-500, 2, 10), normPower = 0)
          )
          s1 <- c(s1, Biostrings::DNAStringSet(
            rep(paste0(rep.int("-", max_width), collapse = ""), sum(s == "-"))
          ))
        } else {
          s1 <- s
        }
      }  else {
        # get biostringset which fills all positions with gaps where it is not complete, i.e. width is < max width
        dels <- Biostrings::DNAStringSet(vapply(max(Biostrings::width(s)) - Biostrings::width(s), function(x) {
          paste0(rep.int("-", x), collapse = "")
        }, FUN.VALUE = ""))
        # merge both biostrings so each seq is of same width
        s1 <- Biostrings::xscat(s, dels)
      }
      cm <- t(Biostrings::consensusMatrix(s1))[, colnm, drop = FALSE]
      cmf <- sweep(cm, 1, rowSums(cm), "/")
      cmf[, "-"] <- 0
      cm <- cm[apply(cmf, 1, function(row) any(row > threshold)), , drop = FALSE]
      res <- c( res, list(cm) )
    }
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
  if (! "X" %in% colnames(x$consmat)) {
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
  if (length(ids_ <- as.integer(names(cm_attr$read_ids))) > 0) {
    ## shift polymorphic positions with insertions
    run <- rle(ins_run)
    vals <- ins_idx[match(run$values, ins_run)]
    for (i in seq_along(vals)) {
      ids_[ids_ >= vals[i] - 1] <- ids_[ids_ >= vals[i] - 1] + run$lengths[i]
    }
    names(cm_attr$read_ids) <- as.character(ids_)
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
msa_from_bam <- function(bamfile, refseq, paddingLetter = "+", region = NULL){
  nRefseq <- names(refseq)
  lRefseq <- length(refseq[[1]])
  if (is.null(region))
    region <- paste0(nRefseq, ":1-", lRefseq)
  else
    region <- paste(nRefseq, region, sep = ":")
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
