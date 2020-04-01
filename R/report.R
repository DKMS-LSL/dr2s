#' @export
report.DR2S <- function(x, whichMap = NULL, threshold = NULL, ...) {
  ## Collect start time for report runstats
  start.time <- Sys.time()

  flog.info("# report", name = "info")
  if (is.null(threshold))
    threshold <- x$getThreshold()

  ## Export report config to function environment
  args <- x$getOpts("report")
  list2env(args, envir = environment())
  assert_that(
    exists("blockWidth") && is.numeric(blockWidth),
    exists("remap") && is.logical(remap),
    exists("createIgv") && is.logical(createIgv),
    is.double(threshold),
    exists("checkHpCount") && is.logical(checkHpCount),
    exists("hpCount") && is.numeric(hpCount)
  )

  outdir <- .dirCreateIfNotExists(x$absPath("report"))
  .reportMap_(x, map = whichMap, outdir = outdir, blockWidth = blockWidth,
              remap = remap, createIgv = createIgv, ...)

  ## set report runstats
  .setRunstats(x, "report",
               list(Runtime = format(Sys.time() - start.time)))
  flog.info("Done", name = "info")

  return(invisible(x))
}

# debug
#map = "mapFinal"
.reportMap_ <- function(x, map, outdir, blockWidth, remap = TRUE, createIgv = TRUE, ...) {
  if (is.null(map)) {
    if (x$hasMapFinal()) {
      map <- "mapFinal"
   } else if (x$hasMapIter()) {
      map <- "mapIter"
   } else stop("Nothing to report")
  }
  map <- match.arg(tolower(map), c("mapfinal", "mapiter"))

  ## Initiate indenter
  indent <- indentation(1)
  ## write latest consensus sequence to consensus accessor
  seqs <- x$getLatestRef()
  names(seqs) <- "hap" %<<% x$getHapTypes()
  x$consensus$seq <- seqs

  threshold <- max(x$getThreshold(), 0.3)
  if (remap) {
    x <- remapAndReport(x, report = TRUE, threshold = threshold)
  }

  writeReportFiles(x = x, map = map, blockWidth = blockWidth, outdir = outdir)
  cache(x)
  invisible(x)
}

#' Report haplotype consensus files as verified
#'
#' @param x  A \code{\link[=DR2S_]{DR2S}} object.
#' @param map Which mapping do we want to update.
#' @return Side effect
#' @details
#' \code{reportCheckedConsensus} will create a \code{checked} directory
#' containing the consensus sequences
#' \code{\{A,B\}.\{readtype\}.\{mapper\}.checked.fa}
#' and a html alignment file \code{aln.\{readtype\}.\{mapper\}.checked.html}.
#' @family DR2S mapper functions
#' @export
reportCheckedConsensus <- function(x, map = "mapFinal") {
  map <- match.arg(tolower(map), c("mapfinal", "mapiter"))
  ending <- ifelse(length(x$getHapTypes()) == 2,
                   "psa",
                    "msa")
  if (x$hasShortreads()) {
    haplotypes <- "SR" %<<% x$getHapTypes()
    readtype <- x$getSrdType()
    mapper <-  x$getSrdMapper()
  } else {
    haplotypes <- "LR" %<<% x$getHapTypes()
    readtype <- x$getLrdType()
    mapper <-  x$getLrdMapper()
  }
  pairfileUnchecked <- dot(c(map, "aln", readtype, mapper, "unchecked", ending))
  pairfileUnchecked <- normalizePath(x$absPath(file.path("report", pairfileUnchecked)),
                                     mustWork = TRUE)
  assert_that(is.readable(pairfileUnchecked))
  pairfileChecked <- dot(c(map, "aln", readtype, mapper, "checked", ending))
  pairfileChecked <- normalizePath(x$absPath(file.path("report", pairfileChecked)),
                                   mustWork = FALSE)

  if (!file.exists(pairfileChecked)) {
    msg <- sprintf("'%s' must be saved as '%s' before invoking this function",
                   pairfileUnchecked, pairfileChecked)
    stop(msg, call. = FALSE)
  }

  outdir <- .dirCreateIfNotExists(x$absPath("checked"))
  rs <- readPairFile(pairfileChecked)
  seqs <- Biostrings::DNAStringSet(
    lapply(rs, function(s) Biostrings::DNAString(stripIndel(s))))

  ## Alignment
  ref <- x$getRefSeq()
  names(ref) <- strsplitN(names(ref), "~", 1, fixed = TRUE)

  seqsAll <- c(ref, seqs)
  alnFile <-  dot(c("aln", x$getLrdType(), x$getLrdMapper(), "html"))
  .browseAlign(seqsAll, file = file.path(outdir,alnFile), openURL = FALSE)

  ## Export FASTA
  files <- vapply(seq_along(seqs), function(sq, seqs, x) {
    file <- dot(c(x$getSampleId(), x$getLocus(), names(seqs[sq]),
                  x$getLrdType(), x$getLrdMapper(), "fa"))
    seq <- seqs[sq]
    seqname <- paste(x$getSampleId(), sub("^hap", "", names(seq)), sep = "_")
    sampleDetails <- x$getSampleDetails()
    names(seq) <-  paste(seqname,
                         semicolon(c("haplotype=" %<<% litQuote(sub("^hap", "", names(seq))),
                                     sampleDetails,
                                     "date=" %<<% litQuote(Sys.Date()),
                                     "status=" %<<% litQuote("checked"))))
    Biostrings::writeXStringSet(
      seq,
      filepath = file.path(outdir,file),
      format = "fasta"
    )
    file
    }, seqs = seqs, x = x, FUN.VALUE = character(1)
  )
}

#' Manually check the alignment of consensus sequences using an editor.
#'
#'
#' @param x  A \code{\link[=DR2S_]{DR2S}} object.
#' @param map Which stage of mapping to use. Either "mapFinal" (default) or
#' "mapIter".
#' @param where Where to focus with the editor.
#' @param editor Which editor to use. Either "xdg-open" (default) for standard
#' system editor, "subl" for sublime, "gvim" or "gedit"
#' @param openEditor should the editor be opened or just create the files?
#' @return The path to the newliy created alignment file to check the polished
#'   consensus sequences.
#' @export
checkAlignmentFile <- function(x, map = "mapFinal", where = 0,
                                 editor = "xdg-open", openEditor = TRUE) {
  map <- match.arg(tolower(map), c("mapfinal", "mapiter"))
  ending <- ifelse(length(x$getHapTypes()) == 2,
                   "psa",
                   "msa")
  if (x$hasShortreads()) {
    haplotypes <- "SR" %<<% x$getHapTypes()
    readtype <- x$getSrdType()
    mapper <-  x$getSrdMapper()
  } else {
    haplotypes <- "LR" %<<% x$getHapTypes()
    readtype <- x$getLrdType()
    mapper <-  x$getLrdMapper()
  }
  pairfileUnchecked <- dot(c(map, "aln", readtype, mapper, "unchecked", ending))
  pairfileUnchecked <- normalizePath(x$absPath(file.path("report", pairfileUnchecked)),
                                     mustWork = TRUE)
  assert_that(is.readable(pairfileUnchecked))
  pairfileChecked <- dot(c(map, "aln", readtype, mapper, "checked", ending))
  pairfileChecked <- normalizePath(x$absPath(file.path("report", pairfileChecked)),
                                   mustWork = FALSE)
  if (!file.exists(pairfileChecked)) {
    file.copy(pairfileUnchecked, pairfileChecked, overwrite = FALSE)
  }
  where <- 23 + 4*floor((where - 1)/80)
  if (openEditor) {
    editor(pairfileChecked, where, useEditor = editor)
  }
  pairfileChecked
}

#' Remap an allele to the checked consensus
#'
#' @param x  A \code{\link[=DR2S_]{DR2S}} object.
#' @param hptype The allele to remap
#' @param report Should the function read a reference from the report folder?
#' Only internal parameter.
#' @return Returns the updated \code{\link[=DR2S_]{DR2S}} object.
#' @return Creates an executable bash file for inspecting the mapping with IGV
#' @family DR2S mapper functions
#' @export
remapAlignment <- function(x, hptype, report = FALSE, createIgv = TRUE,
                           threshold = NULL, ...) {
  # hptype <- "A"
  indent <- list(...)$indent %||% indentation()
  opts <- list(...)$opts
  ## Overide default arguments
  args <- x$getOpts("remap")
  if (!is.null(args)) {
    env  <- environment()
    list2env(args, envir = env)
  }
  mappings <- list(LR = NA, SR = NA)
  reftag <- "remap"
  outdir <- .dirCreateIfNotExists(x$absPath(reftag))
  if (length(x$getHapTypes()) == 1) {
    readpathLR  <- x$getLongreads()
    if (x$hasShortreads())
      readpathSR <- x$getShortreads()
  } else {
    readpathLR <- x$absPath(readpath(x$mapFinal$LR[[hptype]]))
    if (x$hasShortreads())
      readpathSR <- x$absPath(readpath(x$mapFinal$SR[[hptype]]))
  }
  refpath <- if (report) {
    if (x$hasShortreads()) {
      cmat <- x$mapFinal$SR[[hptype]]$pileup$consmat
      gapThreshold <- 1.5 * threshold
    } else {
      cmat <- x$mapFinal$LR[[hptype]]$pileup$consmat
      gapThreshold <- 1.5 * threshold
    }
    seq <- conseq(cmat, "hap" %<<% hptype, type = "prob", suppressAllGaps = FALSE,
                  gapThreshold = gapThreshold)
    file <- dot(c(names(seq), x$getLrdType(), x$getLrdMapper(), "remap", "fa"))
    names(seq) <- names(seq) %<<% " LOCUS=" %<<%
      x$getLocus() %<<% ";REF=" %<<%  x$getReference()
    seqPath <- x$absPath(file.path("remap", file))
    Biostrings::writeXStringSet(seq, filepath = seqPath, format = "fasta" )
    seqPath
  } else {
    .getUpdatedSeqs(x, hptype)
  }
  if (!is.null(attr(refpath, "updated")) & !is.null(x$remap)) {
    flog.info("%sReference of %s is not changed. Keep old mapping", indent(),
              hptype, name = "info")
    return(list(LR = x$remap$LR[[hptype]],
                SR = x$remap$SR[[hptype]]))
  }
  ## Remap long reads to the same reference sequences as short reads
  flog.info("%sRemapping haplotype <%s>", indent(), hptype, name = "info")
  mapgroupLR <- "LR" %<<% hptype
  pileup <- mapReads(
    mapfun = x$getLrdMapFun(), maplabel = reftag, reffile = refpath,
    refname = mapgroupLR, readfile = readpathLR, readtype = x$getLrdType(),
    opts = opts, outdir = outdir, clean = TRUE, force = TRUE, includeDeletions = TRUE,
    includeInsertions = TRUE, callInsertions = FALSE, indent = incr(indent))

  if (createIgv) {
    igv <- createIgvJsFiles(
      refpath(pileup), bampath(pileup), x$getOutdir(), sampleSize = 100,
      fragmentReads = TRUE)
  }

  mappings$LR  <- MapList_(
        ## mapdata
        readpath  = x$relPath(readpathLR),
        refpath   = x$relPath(refpath),
        bampath   = x$relPath(bampath(pileup)),
        # conspath  = self$relPath(conspath),
        pileup    = pileup,
        stats     = list(coverage = .coverage(pileup)),
        ## required metadata
        maplabel  = reftag,
        refname   = mapgroupLR,
        mapper    = x$getLrdMapper(),
        opts      = opts,
        ## additional metadata
        igv       = igv
  )
  ## Map short reads
  if (x$hasShortreads()) {
    mapgroupSR <- "SR" %<<% hptype
    ## Mapper
    pileup <- mapReads(
      mapfun = x$getSrdMapFun(), maplabel = reftag, reffile = refpath,
      refname = mapgroupSR, readfile = readpathSR, readtype = x$getSrdType(),
      opts = opts, outdir = outdir, clean = TRUE, includeDeletions = FALSE,
      includeInsertions = TRUE, callInsertions = FALSE, clip = TRUE, indent = incr(indent))

    if (createIgv) {
      if (!x$getHomozygous()) {
        clusteredReads <- dplyr::pull(tibble::as_tibble(x$srpartition[[hptype]]$srpartition$haplotypes), read)
      } else {
        clusteredReads <- NULL
      }
      igv <- createIgvJsFiles(
        reference = refpath(pileup), bamfile = bampath(pileup),
        outdir = x$getOutdir(), paired = FALSE,
        sampleSize = 100, clusteredReads = clusteredReads)
    } 
    mappings$SR <- MapList_(
      ## mapdata
      readpath  = x$relPath(readpathSR),
      refpath   = x$relPath(refpath),
      bampath   = x$relPath(bampath(pileup)),
      # conspath  = self$relPath(conspath),
      pileup    = pileup,
      stats     = list(coverage = .coverage(pileup)),
      ## required metadata
      maplabel  = reftag,
      refname   = mapgroupSR,
      mapper    = x$getSrdMapper(),
      opts      = opts,
      ## additional metadata
      igv       = igv
    )
  }


  invisible(mappings)
}


# Helpers -----------------------------------------------------------------

.getUpdatedSeqs <- function(x, hptype){
  ending <- ifelse(length(x$getHapTypes()) == 2,
                   "psa", "msa")
  if (x$hasShortreads()) {
    haplotypes <- "SR" %<<% x$getHapTypes()
    readtype <- x$getSrdType()
    mapper <-  x$getSrdMapper()
  } else {
    haplotypes <- "LR" %<<% x$getHapTypes()
    readtype <- x$getLrdType()
    mapper <-  x$getLrdMapper()
  }
  pairfileChecked <- dot(c("mapFinal", "aln", readtype, mapper, "checked", ending))
  pairfileChecked <- normalizePath(x$absPath(file.path("report", pairfileChecked)),
                                   mustWork = FALSE)
  if (!file.exists(pairfileChecked))
    stop("check out the alignment first!")
  rs <- readPairFile(pairfileChecked)


  seqs <- Biostrings::DNAStringSet(
    lapply(rs, function(s) Biostrings::DNAString(stripIndel(s))))

  ## Export FASTA
  seq <- seqs["hap" %<<% hptype]
  file <- dot(c(names(seq), x$getLrdType(), x$getLrdMapper(), "remap", "fa"))
  names(seq) <- names(seq) %<<% " LOCUS=" %<<% x$getLocus() %<<% ";REF=" %<<% x$getReference()
  seqPath <- x$absPath(file.path("remap", file))
  if (file.exists(seqPath)) {
    oldSeq <- Biostrings::readDNAStringSet(seqPath)
    if (identical(as.character(oldSeq), as.character(seq)))
      attr(seqPath, "updated") <- FALSE
  }
  Biostrings::writeXStringSet(seq, filepath = seqPath, format = "fasta")
  seqPath
}

#' @export
readPairFile <- function(pairfile) {
  assert_that(is.readable(pairfile))
  if (endsWith(pairfile, "psa")) {
    rs <- readLines(pairfile)
    ## This assignment relies on the premise that hapA is always used!
    ## This should be true, bcs it is the cluster with the most reads
    rsA <- rs[grepl("^hapA", rs)]
    rsB <- rs[grepl("^hap[B-Z]", rs)]
    hap <- c(
      Biostrings::DNAStringSet(.collapsePairLines_(rsA)),
      Biostrings::DNAStringSet(.collapsePairLines_(rsB))
    )

    names(hap) <- c("hapA", strsplit1(rsB, "\\s")[1])
  } else if (endsWith(pairfile, "msa")) {
    hap <- readMSA(pairfile)
  }

  ## Check for ambiguous positions in sequence
  ambigPositions <- .getAmbigPositions(hap)
  if (length(ambigPositions) > 0) {
    msg <- vapply(seq_along(ambigPositions), function(l, ambigPositions)
      extractAmbigLetters(ambigPositions, names(ambigPositions)[l]),
      ambigPositions = ambigPositions, FUN.VALUE = character(1))
    flog.info(msg, name = "info")
    stop(paste("Check reported reference! Ambiguous positions were found",
               msg, sep = "\n"))
  }
  hap
}

.collapsePairLines_ <- function(x) {
  paste0(vapply(strsplit(x, split = "\\s+"), `[[`, 3, FUN.VALUE = ""),
         collapse = "")
}

#' Read a multiple sequence alignment from file.
#' @param file The absolute path to the file storing the MSA.
#' @details The format  a \code{phylip} like format written by
#' \code{\link{writeMSA}}.
#' @export
readMSA <- function(file) {
  seqs <- c()
  rows <- scan(file, what = "", sep = "\n", strip.white = TRUE,
               quiet = TRUE, blank.lines.skip = TRUE)

  for (i in seq_len(length(rows))) {
    if (!grepl("^[A-Za-z]", rows[i]))
      next
    line <- unlist(strsplit(x = rows[i], "\\s+"))
    name <- line[1]
    if (is.null(seqs) || is.na(seqs[name]))
      seqs[name] <- ""
    seq <- paste0(line[3:(length(line) - 1)], collapse = "")
    seqs[name] <- paste0(seqs[name], seq, collapse = "")
  }
  seqs <- Biostrings::DNAStringSet(seqs)
  seqs
}

## HELPER ##
extractAmbigLetters <- function(irange, ambigLetter){
  msg <- sprintf("Found %s; decide for %s at positions: \n", ambigLetter,
                 paste(strsplit1(CODE_MAP()[ambigLetter], ""),
                       collapse = " or "))
  ambigPositionLetter <- irange[[ambigLetter]]

  msg %<<% paste0( vapply(seq_along(ambigPositionLetter), function(x)
    sprintf("%s: %s",
            names(ambigPositionLetter[x]),
            ambigPositionLetter[x]), character(1)), collapse = "\n") %<<% "\n"
}

#' Write a multiple sequence alignment in a phylip like format.
#' Formatted for easily accessing positions and differences
#' @param aln The alignment as a \code{\link[Biostrings]{DNAStringSet}}.
#' @param file The filename and path.
#' @param block.width The block width of the resulting alignment.
#' @examples
#' \dontrun{
#' library(Biostrings)
#' msa <- DNAStringSet(c("AAA", "AGA", ))
#' print("TODO: write rep seq example which runs")
#' }
#' @export
writeMSA <- function(aln, file="", block.width = 50){
  # ToDo add check for length
  if (!is(aln, "DNAStringSet"))
    stop("'aln' must be a DNAStringSet object and each sequence of same length")
  if (!is.numeric(block.width))
    stop("'block.width' must be a single number")

  .getConsPosition <- function(pos) {
   if (length(unique(pos)) == 1) {
     return(".")
   } else if ("-" %in% pos) {
     return("+")
   } else {
     return("|")
   }
  }

  alnMat <- as.matrix(aln)
  alnStatLine <- Biostrings::BStringSet(paste0(apply(alnMat, 2, function(x)
    .getConsPosition(x)), collapse = ""))
  alignment <- c(Biostrings::BStringSet(aln), alnStatLine)

  ## Write Header
  cat("#=======================================\n", file = file)
  cat("#\n", file = file, append = TRUE)
  cat("# Aligned_sequences: ", length(aln)," \n", file = file, append = TRUE)
  vapply(seq_along(aln), function(x, file) {
    cat(sprintf("# %s: %s\n",
                x, names(aln[x])),
                file = file, append = TRUE)
    TRUE
  }, file = file, FUN.VALUE = logical(1) )
  cat("#\n#\n", file = file, append = TRUE)
  cat("#=======================================\n", file = file,append = TRUE)

  # Write sequences
  lstart <- 1L
  lend <- block.width
  alignmentLength <- Biostrings::width(alignment[1])
  startWidth <- nchar(as.character(1 + alignmentLength))
  nameWidth <- max(20L, nchar(names(alignment)))
  nblock <- ceiling(alignmentLength / block.width)
  for (i in seq_len(nblock)) {
    to <- i * block.width
    from <- to - block.width + 1L
    if (to > alignmentLength)
      to <- alignmentLength
    names <- names(alignment)
    strings <- Biostrings::subseq(alignment, from, to)

    ## Split the seq every 10 chars
    sp <- "(.{10})"
    addSp <- "\\1 "
    a <- vapply(seq_len(length(alignment) - 1), function(x, nameWidth, lstart,
                                                         startWidth, sp, addSp,
                                                         lend, file) {
      cat(format(names[x], width = nameWidth), " ",
          format(lstart, justify = "right", width = startWidth), " ",
          gsub(sp, addSp, Biostrings::toString(strings[x])), " ",
          format(lend, justify = "right"), "\n",
          sep = "", file = file, append = TRUE)
      TRUE
    }, nameWidth = nameWidth, lstart = lstart, startWidth = startWidth, sp = sp,
    addSp = addSp, lend = lend, file = file, FUN.VALUE = logical(1))
    cat(format(" ", width = nameWidth), " ",
          format(" ", justify = "right", width = startWidth), " ",
          gsub(sp, addSp, Biostrings::toString(strings[length(alignment)])),
        "\n\n", sep = "", file = file, append = TRUE)

    lstart <- lend + 1
    lend <- lstart + block.width - 1
  }
}

#' @export
remapAndReport <- function(x, report = FALSE, threshold = NULL, plot = TRUE, ...) {
  if (is.null(threshold))
    threshold <- max(0.3, x$getThreshold())
  dots   <- list(...)
  indent <- dots$indent %||% indentation()
  ## Export report config to function environment
  args <- x$getOpts("report")
  args$createIgv <- dots$createIgv %||% args$createIgv
  args$checkHpCount <- dots$checkHpCount %||% args$checkHpCount
  args$hpCount <- dots$hpCount %||% args$hpCount
  list2env(args, envir = environment())
  assert_that(
    exists("createIgv") && is.logical(createIgv),
    is.double(threshold),
    exists("checkHpCount") && is.logical(checkHpCount),
    exists("hpCount") && is.numeric(hpCount),
    is.logical(report)
  )
  flog.info("%sRemapping final sequences", indent(), name = "info")

  # bpparam <- BiocParallel::MulticoreParam(workers = .getIdleCores())
  # bpparam <- BiocParallel::MulticoreParam(workers = .getIdleCores())
  mappings <- parallel::mclapply(x$getHapTypes(), function(h, x, report,
                                                           createIgv, threshold) {
    remapAlignment(x, h, report = report, createIgv = createIgv, threshold)
  }, x = x, report = report, createIgv = createIgv, threshold = threshold, mc.cores = .getIdleCores())
  names(mappings) <- x$getHapTypes()
  x$remap <- purrr::transpose(mappings)
  hptypes   <- x$getHapTypes()
  if (plot && !is.null(x$plotRemapSummary)) {
    flog.info("%sPlot remap summary", indent(), name = "info")

    ## Coverage and base frequency
    readtypes <- if (x$hasShortreads()) c("LR", "SR") else "LR"
    plotRows  <- if (x$hasShortreads()) 2 else 1
    ## readtype = "LR"
    plotlist <- foreach(readtype = readtypes) %do% {
      suppressWarnings(x$plotRemapSummary(readtype = readtype, thin = 0.25, width = 20))
    }
    p <- cowplot::plot_grid(plotlist = plotlist, nrow = plotRows, labels = readtypes, hjust = -0.25)
    cowplot::save_plot(p, dpi = 150, filename = x$absPath("plot.remap.png"),
                       base_width = 6*plotRows*length(hptypes),
                       base_height = 6/plotRows*length(readtypes),
                       title = dot(c(x$getLocus(), x$getSampleId())))
    cowplot::save_plot(p, filename = x$absPath(".plots/plot.remap.svg"),
                       base_width = 6*plotRows*length(hptypes),
                       base_height = 6/plotRows*length(readtypes))


  }

  ## Do what polish did
  flog.info("%sReport variants", indent(), name = "info")
  vars <- .callVariants(x, threshold)
  ## Polish the reference. This is necessary for Intron 2 of KIR genes and the
  ## beginning of MICB, because
  ## shortreads are not well mapped. Thus set each base of the reference at
  ## positions 400 to 1,500 to the major base, if the longreads show no
  ## variation and the base is supported with at least 90%.
  if (x$getLocus() %in% c(paste0("KIR", KIR_LOCI()), "MICB")) {
    polished <- lapply(hptypes, function(hptype, seq, vars, report) {
      # hap <- paste0("hap", hptype)
      seq <- if (report) {
        x$consensus$seq[[paste0("hap", hptype)]]
      } else {
        unlist(Biostrings::readDNAStringSet(.getUpdatedSeqs(x, hptype)))
      }
      polishRange <- POLISH_RANGE(x$getLocus())
      ## filter variants of the allele of interest
      ## Use only haplotypes in the polishRange range, no gap variants and LR
      ## support > 80%
      hapVar <- dplyr::filter(vars,
                              .data$haplotype == hptype,
                              grepl("Ambiguous position in short reads", .data$warning),
                              .data$pos %in% polishRange,
                              .data$refSR != "-",
                              .data$altSR != "-",
                              .data$refLR != "-",
                              .data$altLR != "-",
                              .data$supportLR > 0.80)  %>%
        dplyr::mutate(pos = as.numeric(.data$pos))
      if (NROW(hapVar) > 0) {
        bases <- apply(hapVar, 1, function(var) {
          # var <- hapVar[1,]
          seqBase <- as.character(seq[as.numeric(var[["pos"]])])
          # DECIPHER::BrowseSeqs(Biostrings::DNAStringSet(seq))
          if (grepl(var[["refLR"]], CODE_MAP()[seqBase])) {
            return(var[["refLR"]])
          }
          flog.warn("Polishing base at %s is not possible. The base should include %s,
               but is %s", var[["pos"]], var[["refLR"]], seqBase, name = "info")
          warning(sprintf("Polishing base at %s is not possible. The base should include %s,
               but is %s", var[["pos"]], var[["refLR"]], seqBase))
          return(NA_character_)
        })

        positions <- hapVar$pos[!is.na(bases)]
        if (report) {
          seq <- Biostrings::replaceLetterAt(seq, at = positions, na.omit(bases))
        }
        return(list(seq = seq, rmVariants = dplyr::filter(hapVar, pos %in% positions)))
      } else {
        return(list(seq = seq, rmVariants = hapVar))
      }
    }, vars = vars, seq = seq, report = report)
    polished <- purrr::transpose(polished)
    polishedSeqs <- Biostrings::DNAStringSet(polished$seq)
    names(polishedSeqs) <- paste0("hap", hptypes)
    x$consensus$seq <- polishedSeqs
    vars <- dplyr::setdiff(vars, dplyr::bind_rows(polished$rmVariants))
  }
  if (report) {
    seqs <- x$consensus$seq
    names(seqs) <- x$getHapTypes()
    ambigPos <- .getAmbigPositions(seqs)
    if (length(ambigPos) > 0) {
      vars <- dplyr::bind_rows(vars, lapply(seq_along(ambigPos), function(iPos, ambigPos) {
        pos <- ambigPos[[iPos]]
        ambigLetter <- names(ambigPos)[iPos]
        tibble::tibble(
          haplotype = names(pos),
          pos = IRanges::start(pos),
          ref = "",
          alt = "",
          warning = paste0("Ambiguous position ",  ambigLetter,
                           " in consensus | Decide for ",
                           paste(strsplit1(CODE_MAP()[ambigLetter], ""), collapse = " or ")),
          refSR = "",
          altSR = "",
          refLR = "",
          altLR = "")
      }, ambigPos = ambigPos))
    }
  }
  ## Check long isertions
  insPosVars <- dplyr::bind_rows(lapply(hptypes, function(hp, x) {
    bamRel <- x$remap$LR[[hp]]$bampath
    if (is.null(bamRel)) {
      return(tibble::tibble(
        haplotype = NA_character_,
        pos = NA_character_,
        ref = NA_character_,
        alt = NA_character_,
        warning = NA_character_,
        refSR = NA_character_,
        altSR = NA_character_,
        refLR = NA_character_,
        altLR = NA_character_,
        supportSR = NA_character_,
        supportLR = NA_character_
      ))
    }
    bamfile <- x$absPath(x$remap$LR[[hp]]$bampath)
    len <- Biostrings::width(Biostrings::readDNAStringSet(
      x$absPath(x$remap$LR[[hp]]$refpath)))
    winLen = 50
    starts <- seq(1, len - winLen, winLen)
    param <- Rsamtools::ScanBamParam(what=c("rname", "pos", "cigar"))
    bam <- Rsamtools::scanBam(bamfile, param=param)[[1]]
    qr <- GenomicAlignments::cigarRangesAlongQuerySpace(
      bam$cigar, ops = "I")
    rr <- GenomicAlignments::cigarRangesAlongReferenceSpace(
      bam$cigar, ops = "I")
    IRanges::width(rr) <- IRanges::width(qr)
    slidingIfrac <- vapply(starts, function(start, rr, winLen) {
      median(sum(IRanges::width(rr[Biostrings::match(IRanges::start(rr), (start:(start + winLen)), nomatch = 0) > 0 ]))/winLen)
    }, numeric(1), rr = rr, winLen = winLen)
    # plot(slidingIfrac)
    # any(slidingIfrac > 0.10)
    tibble::tibble(
      haplotype = hp,
      pos = starts[which(slidingIfrac > 0.1)] + winLen,
      ref = NA_character_,
      alt = NA_character_,
      warning = "Potential long insertion in long reads",
      refSR = NA_character_,
      altSR = NA_character_,
      refLR = NA_character_,
      altLR = NA_character_,
      supportSR = NA_character_,
      supportLR = NA_character_
      )
  }, x = x))
  vars <- dplyr::bind_rows(vars, insPosVars)

  ## Check homopolymer count; Only check if the count is found in both
  if (checkHpCount & x$hasShortreads()) {
    x <- checkHomopolymerCount(x, hpCount = hpCount)
    for (hp in hptypes) {
      homopolymers <- x$consensus$homopolymers[[hp]]
      diffHp <- which(homopolymers$cons != homopolymers$mode)
      adjHp <- lapply(diffHp, function(i, homopolymers, hp) {
        row <- homopolymers[i, ]
        tibble::tibble(haplotype = hp,
                       pos       = row$position,
                       warning   = sprintf(paste(
                         "Homopolymer is of",
                                   "length %s, but should be %s"),
                                   row$cons, row$mode),
                       ref = "",
                       alt = "",
                       refSR     = "",
                       altSR     = "",
                       refLR     = "",
                       altLR     = "")
      }, homopolymers = homopolymers, hp = hp)
      if (!NROW(adjHp) == 0)
        vars <- dplyr::bind_rows(vars, adjHp)
    }
  }

  ## Report problematic Variants
  probvarFile <- dot(c("problems", "tsv"))
  outdir <- .dirCreateIfNotExists(x$absPath("report"))
  vars <- dplyr::arrange(vars, as.numeric(.data$pos), .data$haplotype)
  readr::write_tsv(vars, path = file.path(outdir, probvarFile),
                   append = FALSE, col_names = TRUE)
  x$consensus$variants <- vars
  ## Write the igv config files
  createIgvConfigs(x, map = "remap", open = FALSE)
  return(invisible(x))
}

writeReportFiles <- function(x, map, blockWidth, outdir) {
  ref <- Biostrings::BStringSet(x$getRefSeq())
  if (x$hasShortreads()) {
    haplotypes <- "SR" %<<% x$getHapTypes()
    readtype <- x$getSrdType()
    mapper <-  x$getSrdMapper()
  } else {
    haplotypes <- "LR" %<<% x$getHapTypes()
    readtype <- x$getLrdType()
    mapper <-  x$getLrdMapper()
  }

  ## Write html alignment file
  alnFile <- dot(c(map, "aln", readtype, mapper, "unchecked"))
  # get all seqs as stringset
  seqs <- x$consensus$seq
  refseqs <- c(ref, seqs)
  .browseAlign(refseqs, file = file.path(outdir, alnFile), openURL = FALSE)

  ## Write consensus FASTA files
  # hp = "SRA"
  # hp = "SRB
  for (hp in haplotypes) {
    hp0 <- substr(hp, 3, 3)
    hpFile <- dot(c(map, hp, readtype, mapper, "unchecked.fa"))
    seq <- seqs[endsWith(names(seqs), hp0)]
    seqname <- paste(x$getSampleId(), hp0, sep = "_")
    sampleDetails <- x$getSampleDetails()
    names(seq) <- paste(seqname,
                        semicolon(c(
                          "haplotype=" %<<% litQuote(hp0),
                          sampleDetails,
                          "date=" %<<% litQuote(Sys.Date()),
                          "status=" %<<% litQuote("unchecked"))))
    Biostrings::writeXStringSet(
      seq,
      filepath = file.path(outdir, hpFile),
      format = "fasta"
    )
  }

  ## Write Pairwise or Multiple Alignment
  if (length(x$getHapTypes()) == 2) {
    alnFile <- dot(c(map, "aln", readtype, mapper, "unchecked", "psa"))
    aln <- Biostrings::pairwiseAlignment(pattern = seqs[1], subject = seqs[2], type = "global")
    Biostrings::writePairwiseAlignments(aln, file.path(outdir, alnFile), block.width = blockWidth)
  } else {
    alnFile <- dot(c(map, "aln", readtype, mapper, "unchecked", "msa"))
    if (length(x$getHapTypes()) > 2) {
      aln <- DECIPHER::AlignSeqs(Biostrings::DNAStringSet(seqs), verbose = FALSE)
    } else {
      aln <- Biostrings::DNAStringSet(seqs[1])
    }
    writeMSA(aln, file = file.path(outdir, alnFile))
  }
}

# Helpers -----------------------------------------------------------------


.getVariants <- function(variants) {
  if (all(vapply(variants, function(x) length(x) == 0, logical(1)))) {
    return(
      tibble::tibble(
        haplotype = character(0), pos = integer(0), ref = character(0),
        alt = character(0), warning = character(0), refSR = character(0),
        altSR = character(0), refLR = character(0), altLR = character(0)
      )
    )
  }

  dplyr::bind_rows(lapply(names(variants), function(h)
    dplyr::bind_rows(lapply(variants[[h]],
                            .extractVariant_, h = h))))
}


.extractVariant_ <- function(v, h) {
  data.frame(
    haplotype = h,
    pos = attr(v, "position") %||% NA,
    ref = v[[1]] %||% NA,
    alt = v[[2]] %||% NA,
    warning = attr(v, "warning") %||% NA,
    refSR = attr(v, "srBases")[[1]] %||% NA,
    altSR = attr(v, "srBases")[[2]] %||% NA,
    refLR = attr(v, "lrBases")[[1]] %||% NA,
    altLR = attr(v, "lrBases")[[2]] %||% NA,
    supportSR = attr(v, "srSupport") %||% NA,
    supportLR = attr(v, "lrSupport") %||% NA,
    stringsAsFactors = FALSE
  )
}
