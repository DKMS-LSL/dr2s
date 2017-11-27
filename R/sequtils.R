#' Extract (a subset of) reads from a bamfile
#'
#' @param x Path to bamfile.
#' @param qnames Sample a fixed subset of reads.
#' @param n Sample a fixed number of n reads.
#' @param replace Sample with or without replacemant?
#'
#' @return A \code{\linkS4class{ShortReadQ}} object
#' @export
#'
#' @examples
#' ###
extract_fastq.default <- function(x, qnames = NULL, n = 100, replace = TRUE) {
  stopifnot(is.character(x) && length(x) == 1)
  bam <- Rsamtools::scanBam(x)[[1]]
  qn <-
    if (!is.null(qnames) && is.null(n)) {
      qnames
    } else if (!is.null(qnames) && !is.null(n)) {
      qn <- sample(qnames, size = n, replace = replace)
      q  <- attr(qnames, "q")[match(qn, qnames)]
      o  <- order(q, decreasing = TRUE)
      qn <- qn[o]
      q  <- q[o]
      attr(qn, "q") <- q
      qn
    } else if (is.null(qnames) && !is.null(n)) {
      sample(bam$qname, size = n, replace = replace)
    } else if (is.null(qnames) && is.null(n)) {
      bam$qname
    }
  trim_softclipped_ends(bam, qn)
}

filter_reads <- function(bam, qnames = NULL, preserve_ref_ends = FALSE) {
  #message("Entering trim_softclipped_ends")
  #on.exit(message("Exiting trim_softclipped_ends"))
  i <- if (is.null(qnames)) {
    seq_along(bam$qname)
  } else {
    match(qnames, bam$qname)
  }
  if (preserve_ref_ends) {
    not_trim_starts <- which(bam$pos[i] == 1)
    not_trim_ends   <- which(bam$pos[i] > max(bam$pos[i]) - 250)
  }
  cigar <- bam$cigar[i]
  id    <- if (!is.null(attr(qnames, "q"))) {
    paste0(bam$qname[i], " q=", round(attr(qnames, "q"), 4))
  } else bam$qname[i]
  sc <- map_softclip(cigar)
  if (preserve_ref_ends) {
    sc[not_trim_starts, "start"] <- NA_integer_
    sc[not_trim_ends, "end"] <- NA_integer_
  }

  trim <- sc %>%
    dplyr::transmute(keep = is.na(start) & is.na(end))

  readlength <- 250
  alignmentScoreThreshold = 0.2*readlength # Change to dynamic
  badScore <- bam$qname[bam$tag$AS < alignmentScoreThreshold]
  trim <- id[trim==FALSE]
  flog.info(" Filter %s softclipping reads of %s total ...", length(trim), length(id), name = "info")
  flog.info(" Filter %s reads with alignment score < %s ...", length(badScore), alignmentScoreThreshold, name = "info")
  # try only discrding by quality
  # unique(c(trim, badScore))
  unique(badScore)

}

trim_softclipped_ends <- function(bam, qnames = NULL, preserve_ref_ends = FALSE) {
  #message("Entering trim_softclipped_ends")
  #on.exit(message("Exiting trim_softclipped_ends"))
  i <- if (is.null(qnames)) {
    seq_along(bam$qname)
  } else {
    match(qnames, bam$qname)
  }
  if (preserve_ref_ends) {
    not_trim_starts <- which(bam$pos[i] == 1)
    not_trim_ends   <- which(bam$pos[i] > max(bam$pos[i]) - 250)
  }
  cigar <- bam$cigar[i]
  sread <- bam$seq[i]
  qual  <- bam$qual[i]
  id    <- if (!is.null(attr(qnames, "q"))) {
    paste0(bam$qname[i], " q=", round(attr(qnames, "q"), 4))
  } else bam$qname[i]
  sc <- map_softclip(cigar)
  if (preserve_ref_ends) {
    sc[not_trim_starts, "start"] <- NA_integer_
    sc[not_trim_ends, "end"] <- NA_integer_
  }
  s  <- sc$start + 1L
  e  <- Biostrings::width(sread) - sc$end
  rs <- ShortRead::ShortReadQ(
    sread = XVector::subseq(sread, s, e),
    quality = XVector::subseq(qual, s, e),
    id = Biostrings::BStringSet(id)
  )
  if (preserve_ref_ends) {
    rs@sread@metadata <- list(
      not_trim_starts = not_trim_starts,
      not_trim_ends   = not_trim_ends
    )
  }
  rs
}

map_softclip <- function(cigar) {
  m  <- gregexpr("^\\d+S", cigar)
  sc <- regmatches(cigar, m)
  sc[!vapply(sc, function(x) length(x) > 0, FUN.VALUE = FALSE)] <- NA_character_
  sc_start <- as.integer(gsub("S", "", unlist(sc)))
  m  <- gregexpr("\\d+S$", cigar)
  sc <- regmatches(cigar, m)
  sc[!vapply(sc, function(x) length(x) > 0, FUN.VALUE = FALSE)] <- NA_character_
  sc_end <- as.integer(gsub("S", "", unlist(sc)))
  tibble::tibble(
    start = sc_start,
    end   = sc_end
  )
}

#' Filter a bamfile by qnames
#'
#' @param x Path to bamfile.
#' @param qnames Sample a fixed subset of reads.
#' @param dest Location where the filtered output file will be created.
#'
#' @return A character vector
#' @export
#'
#' @examples
#' ###
filter_bam <- function(x, qnames, dest) {
  stopifnot(is.character(x) && length(x) == 1)
  filt <- S4Vectors::FilterRules(list(function(x) x$qname %in% qnames))
  Rsamtools::filterBam(x, dest, filter = filt, indexDestination = TRUE)
}

inject_deletions <- function(seq) {
  mdata <- metadata(seq)
  if (!is.null(mdata$deletions) || length(mdata$deletions) > 0) {
    dels <- itertools::ihasNext(iter(mdata$deletions))
    while(itertools::hasNext(dels)) {
      d <- nextElem(dels)
      XVector::subseq(seq, d, d - 1) <- "-"
    }
    metadata(seq) <- mdata
  }
  seq
}

#' Generate reference sequences in a FASTA file
#'
#' @param HLA A \code{HLAGene} object
#' @param allele Allele name.
#' @param outdir Output directory
#' @param fullname Truncate allele names
#'
#' @return Path to refseq fasta file.
#' @keywords internal
#' @examples
#' ###
generate_reference_sequence <- function(hla, allele, outdir, fullname = TRUE) {
  if (is.null(allele)) {
    return(NULL)
  }
  assertthat::assert_that(
    file.exists(outdir),
    assertthat::is.dir(outdir),
    assertthat::is.writeable(outdir)
  )
  sref <- foreach(al = allele, .combine = "c") %do% {
    sref <- if (al == "consensus") {
      hla$cons
    } else {
      hla$get_reference_sequence(al)
    }
    if (fullname) {
      names(sref) <- gsub(" +", "_", names(sref))
    } else {
      names(sref) <- gsub("[*:]", "", al)
    }
    sref
  }
  assertthat::assert_that(is(sref, "BStringSet"))
  # workaround for these damn windows filename conventions
  allele_nm <- gsub("[*]", "#", gsub("[:]", "_", paste0(allele, collapse = "~")))
  #allele_nm <- paste0(allele, collapse = "~")
  outpath <- normalizePath(file.path(outdir, paste0(allele_nm, ".ref.fa")), mustWork = FALSE)

  Biostrings::writeXStringSet(sref, outpath)
  outpath
}

stouffers_zscore <- function(z, w = rep(1, length(z))) {
  sum(z*w, na.rm = TRUE)/sqrt(sum(w^2, na.rm = TRUE))
}

merge_conseqs_ <- function(seqs, scores, verbose = TRUE) {
  aln <- DECIPHER::AdjustAlignment(DECIPHER::AlignSeqs(
    seqs, verbose = FALSE, gapOpening = -12, gapExtension = -2,
    iterations = 1, refinements = 1, restrict = c(-500, 2, 10)
  ))
  #browse_align(seqs)
  split1 <- strsplit(toString(aln[[1]]), "")[[1]]
  split2 <- strsplit(toString(aln[[2]]), "")[[1]]
  SEQ1 <- hlatools::ihasNext(iter(split1))
  SEQ2 <- hlatools::ihasNext(iter(split2))
  SCR1 <- hlatools::ihasNext(iter(scores[[1]]))
  SCR2 <- hlatools::ihasNext(iter(scores[[2]]))
  refseq <- vector("character", unique(Biostrings::width(aln)))
  refscr <- vector("numeric", unique(Biostrings::width(aln)))
  i <- 0L
  while(hlatools::hasNext(SEQ1) && hlatools::hasNext(SEQ2)) {
    s1 <- nextElem(SEQ1)
    s2 <- nextElem(SEQ2)
    i <- i + 1L
    if (s1 == s2) {
      z <- c(nextElem(SCR1), nextElem(SCR2))
      refseq[i] <- s1
      refscr[i] <- stouffers_zscore(z)
    } else if (s1 == "-") {
      refseq[i] <- s2
      refscr[i] <- nextElem(SCR2)
      message("   Deletion at seq1 at ", i)
    } else if (s2 == "-") {
      refseq[i] <- s1
      refscr[i] <- nextElem(SCR1)
      message("   Deletion at seq2 at ", i)
    } else {
      scrs <- c(nextElem(SCR1), nextElem(SCR2))
      max.scr <- which.max(scrs)
      refseq[i] <- c(s1, s2)[max.scr]
      refscr[i] <- scrs[max.scr]
      message("   Mismatch ", s1, " (", round(pnorm(scrs[1]), 3), ") <=> ",
              s2, " (", round(pnorm(scrs[2]), 3),") at ", i)
    }
  }
  conseq <- Biostrings::BStringSet(paste0(refseq, collapse = ""))
  names(conseq) <- paste0(names(seqs), collapse = "_")
  metadata(conseq) <- list(zscores = refscr)
  conseq
}
#x <- self
multialign <- function(x, hptype, n, align = list(
  iterations = 0,
  refinements = 0,
  gapOpening = -2,
  gapExtension = -1,
  perfectMatch = 2,
  misMatch = -5,
  terminalGap = -5,
  gapPower = 0,
  normPower = 1
)) {
  stopifnot(is(x$mapIter[["0"]][[hptype]], "mapIter"))
  hptype <- match.arg(hptype, x$getHapTypes())
  readfile <- x$mapIter[["0"]][[hptype]]$reads
  fq <- ShortRead::readFastq(readfile)
  ids  <- as.character(ShortRead::id(fq))
  q    <- as.numeric(DR2S:::strsplitN(DR2S:::strsplitN(ids, " ", 2, fixed = TRUE), "=", 2, fixed = TRUE))
  w    <- Biostrings::width(fq)
  wq   <- (w/max(w)) * q^2
  n    <- if (length(wq) > n) n else length(wq)
  i    <- which(wq > quantile(wq, 1 - n/length(wq)))
  seqs <- ShortRead::sread(fq[i])
  DECIPHER::AlignSeqs(seqs,
                      iterations   = align$iterations,
                      refinements  = align$refinements,
                      gapOpening   = align$gapOpening,
                      gapExtension = align$gapExtension,
                      perfectMatch = align$perfectMatch,
                      misMatch     = align$misMatch,
                      terminalGap  = align$terminalGap,
                      gapPower     = align$gapPower,
                      normPower    = align$normPower,
                      processors   = 8
                      )

}

multialign_consensus <- function(aln) {
  aln
  n      <- length(aln)
  mat    <- t(Biostrings::consensusMatrix(aln))[, c(1:4, 16)]
  ## always accept the first and last base irrespective of how many gaps there are
  mat[1, 5] <- 0
  mat[NROW(mat), 5] <- 0
  ##
  tbl    <- table(mat[, 5]/n)
  cutoff <- as.numeric(names(tbl)[which.min(tbl)])
  mat0   <- mat[mat[, 5]/n <= cutoff, 1:4]/n
  score  <- apply(mat0, 1, max)
  cons   <- Biostrings::BStringSet(paste0(colnames(mat0)[apply(mat0, 1, which.max)], collapse = ""))
  cons@metadata <- list(score = score, cutoff = cutoff, gaptable = tbl)
  cons
}

trim_polymorphic_ends <- function(fq, min_len = 50) {
  #message("Entering trim_polymorphic_ends")
  sr <- fq@sread
  if (length(sr@metadata) > 0) {
    not_trim_starts <- sr@metadata$not_trim_starts
    not_trim_ends   <- sr@metadata$not_trim_ends
    preserve_ref_ends <- TRUE
  } else {
    preserve_ref_ends <- FALSE
  }
  start <- rep(1, length(sr))
  i <- if (preserve_ref_ends) {
    seq_along(start)[-not_trim_starts]
  } else seq_along(start)
  pos   <- 1
  while (length(i) > 0) {
    message("Trimming ", length(i), " reads from start")
    i <- i[Biostrings::width(sr[i]) > pos]
    if (length(i) == 0)
      next
    i <- i[Biostrings::subseq(sr[i], pos, pos) == Biostrings::subseq(sr[i], pos + 1, pos + 1)]
    start[i] <- start[i] + 1
    pos <- pos + 1

  }
  start[start > 1] <- start[start > 1] + 1

  end <- Biostrings::width(sr)
  i   <- if (preserve_ref_ends) {
    seq_along(end)[-not_trim_ends]
  } else seq_along(end)
  pos <- end
  while (length(i) > 0) {
    message("Trimming ", length(i), " reads from end")
    ##i <- i[Biostrings::width(sr[i]) > pos]
    i <- i[pos[i] > 1]
    if (length(i) == 0)
      next
    i <- i[Biostrings::subseq(sr[i], pos[i], pos[i]) == Biostrings::subseq(sr[i], pos[i] - 1, pos[i] - 1)]
    end[i] <- end[i] - 1
    pos[i] <- pos[i] - 1
  }
  end[end > 1] <- end[end > 1] - 1

  ## filter low complexity reads
  if (any(j <- (end - start) <= 0)) {
    remove   <- which(j)
    fq    <- fq[-remove]
    start <- start[-remove]
    end   <- end[-remove]
  }

  fq@sread <- Biostrings::subseq(fq@sread, start, end)
  fq@quality@quality <- Biostrings::subseq(fq@quality@quality, start, end)
  fq[Biostrings::width(fq) >= min_len]
}
get_seqs_from_mat <- function(mat){
  seqs <- apply(mat, 1, function(t) c(unlist(paste(t, collapse = ""))))
  seqs <- seqs[nchar(gsub("-", "",seqs))>0]
  Biostrings::DNAStringSet(seqs)
}
