## Extract (a subset of) reads from a bamfile
.extractFastq <- function(x, qnames = NULL) {
  assert_that(is.string(x))
  bam <- Rsamtools::scanBam(x)[[1]]
  qn <- if (!is.null(qnames)) {
      qnames
    } else if (is.null(qnames)) {
      bam$qname
    }
  .trimSoftclippedEnds(bam, qn)
}

# todo add length and cut option
.filterReads <- function(bam, qnames = NULL, preserveRefEnds = TRUE, ...) {
  indent <- list(...)$indent %||% indentation()
  readlength <- 250
  i <- if (is.null(qnames)) {
    seq_along(bam$qname)
  } else {
    match(qnames, bam$qname)
  }
  if (preserveRefEnds) {
    notTrimStarts <- which(bam$pos[i] == 1)
    notTrimEnds   <- which(bam$pos[i] > max(bam$pos[i]) - readlength)
  }
  cigar <- bam$cigar[i]
  id    <- if (!is.null(attr(qnames, "q"))) {
    bam$qname[i] %<<% " q=" %<<% round(attr(qnames, "q"), 4)
  } else bam$qname[i]
  sc <- .mapSoftclip(cigar)
  if (preserveRefEnds) {
    sc[notTrimStarts, "start"] <- NA_integer_
    sc[notTrimEnds, "end"] <- NA_integer_
  }
  trim <- sc %>%
    dplyr::transmute(keep = ifelse(is.na(.data$start) & is.na(.data$end),
                                   readlength,
                                   readlength -
                                     abs(ifelse(is.na(.data$end), 0,
                                                .data$end) -
                                           ifelse(is.na(.data$start), 0,
                                                  .data$start))))

  alignmentScoreThreshold = 0.3*readlength # Change to dynamic
  badScore <- bam$qname[bam$tag$AS < alignmentScoreThreshold]
  trim <- id[trim < 0.3 * readlength]
  flog.info("%sFiltering %s softclipping reads from %s reads in total;" %<<%
            "removing reads < %s bp", indent(), length(trim), length(id),
            0.3*readlength, name = "info")
  flog.info("%sFiltering %s reads with an alignment score < %s ...",
            indent(), length(badScore), alignmentScoreThreshold, name = "info")
  unique(c(trim, badScore))
}

.trimToRefLen <- function(bamfile, fqFileWrite) {
  bam <- Rsamtools::scanBam(bamfile, param = Rsamtools::ScanBamParam(what = c("qname", "pos", "seq", "qual", "cigar")))[[1]]
  sc <- .mapSoftclip(bam$cigar, which = "start")
  clipPos <- sc - bam$pos
  clipPos[clipPos < 0 | is.na(clipPos)] <- 0
  clipPos[clipPos > 0] <- clipPos[clipPos > 0] + 2
  clipPos[clipPos == 0] <- 1
  fqrClip <- ShortRead::ShortReadQ(sread = bam$seq, quality = Biostrings::BStringSet(bam$qual), id = Biostrings::BStringSet(bam$qname))
  fqrClipped <- ShortRead::narrow(fqrClip, start = clipPos, end = Biostrings::width(bam$seq))
  fqFileWrite <- .fileDeleteIfExists(fqFileWrite)
  ShortRead::writeFastq(fqrClipped, fqFileWrite, compress = TRUE)
  TRUE
}


.trimSoftclippedEnds <- function(bam, qnames = NULL, preserveRefEnds = FALSE) {
  i <- if (is.null(qnames)) {
    seq_along(bam$qname)
  } else {
    match(qnames, bam$qname)
  }
  if (preserveRefEnds) {
    notTrimStarts <- which(bam$pos[i] == 1)
    notTrimEnds   <- which(bam$pos[i] > max(bam$pos[i]) - 250)
  }
  cigar <- bam$cigar[i]
  sread <- bam$seq[i]
  qual  <- bam$qual[i]
  id    <- if (!is.null(attr(qnames, "q"))) {
    bam$qname[i] %<<% " q=" %<<% round(attr(qnames, "q"), 4)
  } else bam$qname[i]
  sc <- .mapSoftclip(cigar)
  if (preserveRefEnds) {
    sc[notTrimStarts, "start"] <- NA_integer_
    sc[notTrimEnds, "end"] <- NA_integer_
  }
  s  <- sc$start + 1L
  e  <- Biostrings::width(sread) - sc$end
  rs <- ShortRead::ShortReadQ(
    sread = XVector::subseq(sread, s, e),
    quality = XVector::subseq(qual, s, e),
    id = Biostrings::BStringSet(id)
  )
  if (preserveRefEnds) {
    rs@sread@metadata <- list(
      notTrimStarts = notTrimStarts,
      notTrimEnds   = notTrimEnds
    )
  }
  rs
}

.mapSoftclip <- function(cigar, which = "both") {
  m  <- gregexpr("^\\d+S", cigar)
  sc <- regmatches(cigar, m)
  sc[!vapply(sc, function(x) length(x) > 0, FUN.VALUE = FALSE)] <- NA_character_
  scStart <- as.integer(gsub("S", "", unlist(sc)))
  if (which == "start")
    return(scStart)
  m  <- gregexpr("\\d+S$", cigar)
  sc <- regmatches(cigar, m)
  sc[!vapply(sc, function(x) length(x) > 0, FUN.VALUE = FALSE)] <- NA_character_
  scEnd <- as.integer(gsub("S", "", unlist(sc)))
  if (which == "end")
    return(scStart)
  tibble::tibble(
    start = scStart,
    end   = scEnd
  )
}

#' Generate reference sequences in a FASTA file
#'
#' @param locus Locus name (HLA or KIR).
#' @param allele Allele name or NULL.
#' @param outdir Output directory.
#' @param dirtag TODO
#' @return Path to refseq fasta file.
#' @keywords internal
#' @examples
#' ###
.generateReferenceSequence <- function(locus, allele = NULL, outdir, 
                                       utrLength = NULL, dirtag = NULL) {
  if (is.null(allele)) {
    ## fetch the generic reference for <locus>
    allele <- .normaliseLocus(locus)
    sref <- findRef(allele)
  } else {
    ## fetch the allele-specific reference for <locus>
    locus <- .normaliseLocus(locus)
    if (startsWith(locus, "KIR")) {
      ipd <- .ipdKir()
    } else {
      ipd <- .ipdHla()
    }
    if (startsWith(locus, "MIC"))
      locus <- paste("HLA", locus, sep = "-")
    if (!allele %in% ipd$getAlleles(locus))
      stop(sprintf("Allele %s not found in database", allele))
    sref <- foreach(i = allele, .combine = "c") %do% {
      sref <- ipd$getClosestComplete(i, locus)
      srefFeatures <- ipd$getStructure(names(sref))
      names(sref) <- gsub(" +", "_", names(sref))
      if (!is.null(utrLength)) {
        assert_that(all(c("start", "end") %in% names(utrLength)))
        if (utrLength$start == -1) {
          start <- 1
        } else {
          start <- srefFeatures[srefFeatures$name == "5' UTR"]
          start <- IRanges::end(start) - utrLength$start
        }
        if (utrLength$end == -1) {
          end <- Biostrings::width(sref)
        } else {
          end <- srefFeatures[srefFeatures$name == "3' UTR"]
          end <- IRanges::start(end) + utrLength$end
        }
        sref <- Biostrings::subseq(sref, start = start, 
                                   end = end)
      }
      sref
    }
  }

  assert_that(is(sref, "DNAStringSet"))
  .dirCreateIfNotExists(normalizePath(
    file.path(outdir, dirtag), mustWork = FALSE))
  assert_that(
    file.exists(outdir),
    is.dir(outdir),
    is.writeable(outdir)
  )

  # Strip illegal characters from filenames
  alleleNm <- strip(allele, "_")
  filename <- ifelse(is.null(dirtag), alleleNm %<<% ".ref.fa",
                     file.path(dirtag, alleleNm %<<% ".ref.fa"))
  outpath <- normalizePath(file.path(outdir, filename), mustWork = FALSE)
  Biostrings::writeXStringSet(sref, outpath)
  filename
}

.trimPolymorphicEnds <- function(fq, minLen = 50) {
  #message("Entering .trimPolymorphicEnds")
  sr <- fq@sread
  if (length(sr@metadata) > 0) {
    notTrimStarts <- sr@metadata$notTrimStarts
    notTrimEnds   <- sr@metadata$notTrimEnds
    preserveRefEnds <- TRUE
  } else {
    preserveRefEnds <- FALSE
  }
  start <- rep(1, length(sr))
  i <- if (preserveRefEnds) {
    seq_along(start)[-notTrimStarts]
  } else seq_along(start)
  pos   <- 1
  while (length(i) > 0) {
    #message("Trimming ", length(i), " reads from start")
    i <- i[Biostrings::width(sr[i]) > pos]
    if (length(i) == 0)
      next
    i <- i[Biostrings::subseq(sr[i], pos, pos) == Biostrings::subseq(sr[i],
                                                                     pos + 1,
                                                                     pos + 1)]
    start[i] <- start[i] + 1
    pos <- pos + 1

  }
  start[start > 1] <- start[start > 1] + 1

  end <- Biostrings::width(sr)
  i   <- if (preserveRefEnds) {
    seq_along(end)[-notTrimEnds]
  } else seq_along(end)
  pos <- end
  while (length(i) > 0) {
    message("Trimming ", length(i), " reads from end")
    ##i <- i[Biostrings::width(sr[i]) > pos]
    i <- i[pos[i] > 1]
    if (length(i) == 0)
      next
    i <- i[Biostrings::subseq(sr[i], pos[i],
                              pos[i]) == Biostrings::subseq(sr[i],
                                                            pos[i] - 1,
                                                            pos[i] - 1)]
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
  fq[Biostrings::width(fq) >= minLen]
}

.getSeqsFromMat <- function(mat){
  seqs <- apply(mat, 1, function(t) c(unlist(paste(t, collapse = ""))))
  seqs <- seqs[nchar(stripIndel(seqs, replace = "")) > 0]
  Biostrings::DNAStringSet(seqs)
}

.collapseHomopolymers <- function(path, n = 5) {
  path <- normalizePath(path, mustWork = TRUE)
  assert_that(is.count(n))

  seq <- Biostrings::readBStringSet(path)
  seqname <- names(seq)
  seq <- strsplit1(as.character(seq), "")
  seq <- rle(seq)
  seq$lengths[seq$lengths > n] <- n
  seq <- paste(inverse.rle(seq), collapse = "")
  seq <- Biostrings::DNAStringSet(stripIndel(seq))
  names(seq) <- seqname
  Biostrings::writeXStringSet(seq, path)
  path
}

.mat2rle <- function(mat) {
  assert_that(is(mat, "consmat"))
  .makeProbConsensus_(x = mat, suppressAllGaps = TRUE, asRle = TRUE)
}

.seq2rle <- function(seq) {
  rle(strsplit1(as.character(seq), ""))
}


#' Get a distribution of homopolymer counts in alleles
#'
#' @param x a DR2S object.
#' @param hpCount the minimal length of a homopolymer to be checked (10).
#' @param map Which result to use. Either "mapFinal" or "remap".
#' @return plot a pdf with length histogram and write mode to log
#' @examples
#' ###
#' @export
checkHomopolymerCount <- function(x, hpCount = 10) {
  assert_that(is(x, "DR2S"))
  hptypes <- x$getHapTypes()
  x$consensus$homopolymers <- list()
  plotname <- "plot.homopolymers.pdf"
  #hp <- "A"
  p <- foreach(hp = hptypes) %do% {
    bamfile <- ifelse(x$hasShortreads(),
        x$absPath(bampath(x$remap[["SR"]][[hp]])),
        x$absPath(bampath(x$remap[["LR"]][[hp]])))
    ref <- ifelse(x$hasShortreads(),
      x$absPath(refpath(x$remap$SR[[hp]])),
      x$absPath(refpath(x$remap$LR[[hp]])))
    seq <- Biostrings::readDNAStringSet(ref)
    seqrle <- .seq2rle(seq)
    n <- which(seqrle$lengths >= hpCount)
    if (length(n) == 0) {
      return(NULL)
    }
    bam <- Rsamtools::BamFile(bamfile)
    #pos = n[1]
    homopolymersHP <- foreach(pos = n, .combine = rbind) %do% {
      positionHP <- sum(seqrle$length[seq_len(pos - 1)])
      lenHP <- seqrle$lengths[pos]
      range <- c(positionHP - 10, positionHP + lenHP + 10)
      Rsamtools::open.BamFile(bam)
      msa <- .msaFromBam(bam, range)
      Rsamtools::close.BamFile(bam)
      covering <- vapply(msa, function(a, lenHP) {
        nchar(gsub(pattern = "\\+", "", toString(a))) == lenHP + 21
      }, lenHP = lenHP, FUN.VALUE = logical(1), USE.NAMES = TRUE)
      msa <- msa[covering]
      readsOI <- names(msa)
      ins <- .getInsertions(bamfile,
                            inpos = positionHP,
                            reads = readsOI,
                            readtype = "illumina")
      if (length(ins) > 0) {
        ins <- ins[[1]]
        s <- msa[names(msa) %in% names(ins)]
        pre <- Biostrings::subseq(s, start = 1, width = 11)
        post <- Biostrings::subseq(s, start = 12, end = length(s[[1]]))
        s <- Biostrings::xscat(pre, ins, post)
        msa[names(s)] <- s
      }
      msarle <- lapply(msa, .seq2rle)
      lens <- vapply(msarle, function(a) max(a$lengths), integer(1))
      tibble::tibble(haptype = hp, position = positionHP, length = lens,
                     read = names(lens), cons = lenHP)
    }

    modes <- homopolymersHP %>%
      dplyr::group_by(.data$position, .data$cons) %>%
      dplyr::summarise(mode = .getModeValue(length)) %>% 
      dplyr::ungroup()
    x$consensus$homopolymers[[hp]] <- modes
    cl <- dplyr::pull(x$srpartition$A$srpartition$haplotypes, read)
    homopolymersHP$clustered <- homopolymersHP$read %in% cl
    max(homopolymersHP$length)
    plots <- ggplot(data = homopolymersHP) +
      geom_histogram(aes(length, fill = clustered), binwidth = 1) +
      scale_x_continuous(breaks = seq(min(homopolymersHP$length),
                                      max(homopolymersHP$length), 2)) +
      facet_wrap(. ~ position, scales = "free") +
      ggtitle(paste("Homopolymer length", hp)) +
      cowplot::theme_cowplot() + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
            panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank(),
            panel.grid.major.x = element_line(size = .1, color = "black"),
            panel.grid.minor.x = element_line(size = .005, color = "black") )
    list(n = n, plot = plots)
  }
  
  if (!all(is.null(unlist(p)))) {
    plots <- lapply(p, function(x) x$plot)
    n <- max(vapply(p, function(x) length(x$n), FUN.VALUE = integer(1)))
    p1 <- cowplot::plot_grid(plotlist = plots, ncol = length(n))
    cowplot::save_plot(p1, dpi = 150, filename = x$absPath("plot.homopolymers.png"),
                       base_width = 7*length(hptypes),
                       base_height = 7*length(n),
                       title = dot(c(x$getLocus(), x$getSampleId())))
    cowplot::save_plot(p1, filename = x$absPath(".plots/plot.homopolymers.svg"),
                       base_width = 7*length(hptypes),
                       base_height = 7*length(n))
  }

  return(invisible(x))
}

.getAmbigPositions <- function(seqs) {
  ## Check for ambiguous bases
  seqLetters <- Biostrings::uniqueLetters(seqs)
  ambigLetters <- seqLetters[which(!seqLetters %in% VALID_DNA(include = "del"))]
  ambigPositions <- list()
  if (!length(ambigLetters) == 0) {
    ambigPositions <- stats::setNames(lapply(ambigLetters, function(l, seqs)
      unlist(Biostrings::vmatchPattern(l, seqs)), seqs = seqs), ambigLetters)
  }
  ambigPositions
}
