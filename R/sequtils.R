## Extract (a subset of) reads from a bamfile
.extractFastq <- function(x, qnames=NULL, n=100, replace=TRUE) {
  stopifnot(is.character(x) && length(x) == 1)
  bam <- Rsamtools::scanBam(x)[[1]]
  qn <-
    if (!is.null(qnames) && is.null(n)) {
      qnames
    } else if (!is.null(qnames) && !is.null(n)) {
      qn <- sample(qnames, size=n, replace=replace)
      q  <- attr(qnames, "q")[match(qn, qnames)]
      o  <- order(q, decreasing=TRUE)
      qn <- qn[o]
      q  <- q[o]
      attr(qn, "q") <- q
      qn
    } else if (is.null(qnames) && !is.null(n)) {
      sample(bam$qname, size=n, replace=replace)
    } else if (is.null(qnames) && is.null(n)) {
      bam$qname
    }
  .trimSoftclippedEnds(bam, qn)
}

# todo add length and cut option
.filterReads <- function(bam, qnames=NULL, preserve_ref_ends=TRUE) {
  #message("Entering trim_softclipped_ends")
  #on.exit(message("Exiting trim_softclipped_ends"))
  readlength <- 250
  i <- if (is.null(qnames)) {
    seq_along(bam$qname)
  } else {
    match(qnames, bam$qname)
  }
  if (preserve_ref_ends) {
    not_trim_starts <- which(bam$pos[i] == 1)
    not_trim_ends   <- which(bam$pos[i] > max(bam$pos[i]) - readlength)
  }
  cigar <- bam$cigar[i]
  id    <- if (!is.null(attr(qnames, "q"))) {
    paste0(bam$qname[i], " q=", round(attr(qnames, "q"), 4))
  } else bam$qname[i]
  sc <- .mapSoftclip(cigar)
  if (preserve_ref_ends) {
    sc[not_trim_starts, "start"] <- NA_integer_
    sc[not_trim_ends, "end"] <- NA_integer_
  }
  trim <- sc %>%
    dplyr::transmute(keep = ifelse(is.na(start) & is.na(end),
                                   readlength,
                                   readlength -
                                     abs(ifelse(is.na(end), 0, end) -
                                           ifelse(is.na(start), 0, start))))

  alignmentScoreThreshold = 0.3*readlength # Change to dynamic
  badScore <- bam$qname[bam$tag$AS < alignmentScoreThreshold]
  trim <- id[trim < 0.3 * readlength]
  flog.info(paste0("  Filtering %s softclipping reads from %s reads in total;",
                   "removing reads < %s bp ..."),
            length(trim),
            length(id),
            0.3*readlength,
            name = "info")
  flog.info("  Filtering %s reads with an alignment score < %s ...",
            length(badScore), alignmentScoreThreshold, name = "info")
  unique(c(trim, badScore))
}

.trimSoftclippedEnds <- function(bam, qnames=NULL, preserve_ref_ends=FALSE) {
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
  sc <- .mapSoftclip(cigar)
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

.mapSoftclip <- function(cigar) {
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
generateReferenceSequence <- function(allele, locus, outdir, dirtag=NULL,
                                        fullname=TRUE) {
  if (is.null(allele)) {
    return(NULL)
  }
  locus <- normalise_locus(locus)
  if (!allele %in% ipd.Hsapiens.db::getAlleles(locus)) 
    stop(sprintf("Allele %s not found in database", allele))
    
  dir_create_if_not_exists(normalizePath(
    file.path(outdir, dirtag), mustWork=FALSE))
  assertthat::assert_that(
    file.exists(outdir),
    assertthat::is.dir(outdir),
    assertthat::is.writeable(outdir)
  )
  sref <- foreach(al=allele, .combine="c") %do% {
    sref <- ipd.Hsapiens.db::getClosestComplete(al, locus)
    if (fullname) {
      names(sref) <- gsub(" +", "_", names(sref))
    } else {
      names(sref) <- gsub("[*:]", "", al)
    }
    sref
  }
  assertthat::assert_that(is(sref, "DNAStringSet"))
  # workaround for these damn windows filename conventions
  allele_nm <- gsub("[*]", "#", gsub("[:]", "_", paste0(allele, collapse="~")))
  filename <- ifelse(is.null(dirtag), paste0(allele_nm, ".ref.fa"),
                     file.path(dirtag, paste0(allele_nm, ".ref.fa")))
  outpath <- normalizePath(file.path(outdir, filename), mustWork=FALSE)

  Biostrings::writeXStringSet(sref, outpath)
  filename
}

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

.collapseHomopolymers <- function(path, n = 5) {

  path <- normalizePath(path, mustWork = TRUE)
  assertthat::is.count(n)

  seq <- Biostrings::readBStringSet(path)
  seqname <- names(seq)
  seq <- strsplit1(as.character(seq), "")
  seq <- rle(seq)
  seq$lengths[seq$lengths > n] <- n
  seq <- paste(inverse.rle(seq), collapse = "")
  seq <- Biostrings::DNAStringSet(gsub("[-+]", "", seq))
  names(seq) <- seqname
  Biostrings::writeXStringSet(seq, path)
  path
}

.mat2rle <- function(mat) {
  seq <- conseq(mat, type = "prob", force_exclude_gaps = TRUE)
  .seq2rle(seq)
}

.seq2rle <- function(seq) {
  rle(strsplit1(as.character(seq), ""))
}


#' Get a distribution of homopolymer count in alleles;
#'
#' @param x a DR2S object
#' @param count the minimal length of a homopolymer to be found (10)
#' @param map Which result to use. Either "mapFinal" or "refine"
#'
#' @return plot a pdf with length histogram and write mode to log
#' @examples
#' ###

#' @export
checkHomoPolymerCount <- function(x, count = 10, map = "mapFinal") {
  map <- match.arg(map, c("mapFinal", "refine"))
  stopifnot(is(x, "DR2S"))
  hptypes <- x$getHapTypes()
  if (map == "mapFinal") {
    if (length(hptypes) > 1) {
      lastIter <- x$mapIter[[max(names(x$mapIter))]]
      refseqs <- sapply(hptypes, function(x) lastIter[[x]]$conseq)
      bambase <- x$mapFinal$bamfile
    } else {
      refseqs <- list(A = x$mapInit$SR1$conseq)
      bambase <- x$mapInit$SR2$bamfile
    }
    x$mapFinal$homopolymers <- list()
    plotname <- "plot.homopolymers.pdf"
  } else if ( map == "refine" ) {
    refseqs <- sapply(x$consensus$refine$ref, function(hp) {
      seq <- Biostrings::readDNAStringSet(file.path(x$getOutdir(), hp))
      names(seq) <- strsplit1(names(seq), " ")[1]
      seq
    })
    x$consensus$homopolymers <- list()
    bambase <- x$consensus$refine$bamfile
    plotname <- "plot.homopolymers.refine.pdf"
  }
  p <- foreach(hp = hptypes) %do% {
    seq <- refseqs[[hp]]
    seqrle <- .seq2rle(seq)
    n <- which(seqrle$lengths > count)
    if (length(n) == 0) {
      return(NULL)
    }
    bamfile <- ifelse(length(hptypes) == 1, bambase, bambase[paste0("SR", hp)])
    bamfile <- file.path(x$getOutdir(), bamfile)
    homopolymersHP <- foreach(pos = n, .combine = rbind) %do% {
      positionHP <- sum(seqrle$length[1:(pos-1)])+1
      lenHP <- seqrle$lengths[pos]
      msa <- msa_from_bam(bamfile,
                          refseq = seq,
                          region = sprintf("%s:%s-%s",
                                           names(seq),
                                           positionHP-10,
                                           positionHP + lenHP + 10))


      covering <- sapply(msa, function(a) nchar(gsub(pattern = "\\+",
                                                     "",
                                                     toString(a))) == lenHP + 21)
      msa <- msa[covering]
      readsOI <- names(msa)
      ins <- .getInsertions(bamfile, inpos = positionHP-1, reads = readsOI)[[1]]
      msa <- unlist(Biostrings::DNAStringSetList(sapply(names(msa), function(a) {
        if (a %in% names(ins)) {
          firstSeq <- unlist(Biostrings::subseq(msa[a], start = 1, width = 10))
          insert <- unlist(ins[a])
          lastSeq <- unlist(Biostrings::subseq(msa[a], start = 11, width = lenHP + 10))
          Biostrings::DNAStringSet(c(firstSeq, insert, lastSeq))
        } else {
          msa[a]
        }
      })))
      msarle <- lapply(msa, .seq2rle)
      lens <- sapply(msarle, function(a) max(a$lengths))
      dplyr::data_frame(haptype = hp, position = positionHP, length = lens)
    }

    modes <- homopolymersHP %>%
      dplyr::group_by(position) %>%
      dplyr::summarize(mode = .getModeValue(length))
    for (i in 1:NROW(modes)) {
      modeHP <- modes[i,]
      flog.info("%s: Mode for homopolymer at position %s: %s",
                hp, modeHP$position, modeHP$mode)
      if (map == "mapFinal") {
        x$mapFinal$homopolymers[[hp]][[as.character(modeHP$position)]] <- modeHP$mode
      } else if (map == "refine") {
        x$consensus$homopolymers[[hp]][[as.character(modeHP$position)]] <- modeHP$mode
      }
    }
    ggplot(data = homopolymersHP) +
      geom_histogram(aes(length), binwidth = 1) +
      scale_x_continuous(breaks = seq(min(homopolymersHP$length), max(homopolymersHP$length), 1)) +
      facet_grid(position~., scales = "free_y") +
      ggtitle(paste("Homopolymer length", hp)) +
      theme_minimal()
  }
  if (!all(is.null(unlist(p)))) {
    p1 <- cowplot::plot_grid(plotlist = p, nrow = length(n))
    cowplot::save_plot(p1, filename = x$absPath("plot.Homopolymers.pdf"),
                base_width = 7*length(hptypes),
                title     = paste(x$getLocus(), x$getSampleId(), sep = "." ),
                base_height = 7*length(n))
    cowplot::save_plot(p1, filename = x$absPath(".plots/plot.Homopolymers.svg"),
              base_width = 7*length(hptypes),
              base_height = 7*length(n))
  }
  return(invisible(x))
}
