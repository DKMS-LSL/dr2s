#' @export
report.DR2S <- function(x, whichMap, ...) {
  ## If reporting is already done exit safely
  if (.checkReportStatus(x)) return(invisible(x))

  ## Collect start time for report runstats
  start.time <- Sys.time()

  flog.info("# report", name = "info")

  ## Export report config to function environment
  args <- x$getOpts("report")
  list2env(args, envir = environment())
  assert_that(
    exists("blockWidth") && is.numeric(blockWidth),
    exists("noRemap") && is.logical(noRemap),
    exists("createIgv") && is.logical(createIgv)
  )

  outdir <- .dirCreateIfNotExists(x$absPath("report"))
  if (missing(whichMap)) {
    ## if `which` is unspecified choose `mapFinal` if available,
    ## otherwise try `mapIter`
    if (x$hasMapFinal()) {
      .reportMap_(x, map = "mapFinal", outdir = outdir, blockWidth = blockWidth,
                  noRemap = noRemap, createIgv = createIgv, ...)
    }
    else if (x$hasMapIter()) {
      .reportMap_(x, map = "mapIter", outdir = outdir, blockWidth = blockWidth,
                  noRemap = noRemap, createIgv = createIgv, ...)
    }
    else stop("Nothing to report")
  } else {
    whichMap <- match.arg(tolower(whichMap), c("mapFinal", "mapIter"))
    .reportMap_(x, map = whichMap, outdir = outdir, blockWidth = blockWidth,
                noRemap = noRemap, createIgv = createIgv, ...)
  }

  ## set report runstats
  .setRunstats(x, "report",
               list(Runtime = format(Sys.time() - start.time)))
  flog.info("Done", name = "info")

  return(invisible(x))
}

# debug
#map = "mapFinal"
.reportMap_ <- function(x, map, outdir, blockWidth, noRemap = FALSE, createIgv = TRUE, ...) {
  map <- match.arg(tolower(map), c("mapfinal", "mapiter"))
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

  ## Initiate indenter
  indent <- indentation(1)
  ## Write html alignment file
  alnFile <- dot(c(map, "aln", readtype, mapper, "unchecked"))
  # get all seqs as stringset

  haps <- x$getLatestRef()
  seqs <- Biostrings::BStringSet(haps[grepl(paste0(haplotypes, collapse = "|"), names(haps))])
  refseqs <- c(ref, seqs)
  .browseAlign(refseqs, file = file.path(outdir, alnFile), openURL = FALSE)

  ## Write consensus FASTA files
  # hp = "SRA"
  # hp = "SRB
  for (hp in haplotypes) {
    hp0 <- substr(hp, 3, 3)
    hpFile <- dot(c(map, hp, readtype, mapper, "unchecked.fa"))
    seq <- seqs[endsWith(names(seqs), hp)]
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
  if (length(x$getHapTypes()) == 1) {
    alnFile <- dot(c(map, "aln", readtype, mapper, "unchecked", "fa"))
    Biostrings::writeXStringSet(seqs[1], filepath = file.path(outdir, alnFile), format = "fasta")
  }
  else if (length(x$getHapTypes()) == 2) {
    alnFile <- dot(c(map, "aln", readtype, mapper, "unchecked", "psa"))
    aln <- Biostrings::pairwiseAlignment(pattern = seqs[1], subject = seqs[2], type = "global")
    Biostrings::writePairwiseAlignments(aln, file.path(outdir, alnFile), block.width = blockWidth)
  }
  else if (length(x$getHapTypes()) > 2) {
    alnFile <- dot(c(map, "aln", readtype, mapper, "unchecked", "msa"))
    aln <- DECIPHER::AlignSeqs(Biostrings::DNAStringSet(seqs), verbose = FALSE)
    writeMSA(aln, file = file.path(outdir, alnFile))
  }

  if (map == "mapfinal") {
    ## Report problematic Variants
    probvarFile <- dot(c("problems", readtype, mapper, "tsv"))
    vars <- x$consensus$variants %>%
      dplyr::arrange(as.numeric(.data$pos), .data$haplotype)
    readr::write_tsv(vars, path = file.path(outdir, probvarFile),
                     append = FALSE, col_names = TRUE)

    ## Write postprocessing scripts
    .writeCheckConsensus(path = x$getOutdir())
    .writeReportCheckedConsensus(path = x$getOutdir())
    .writeRefineAlignments(path = x$getOutdir(), haptypes = x$getHapTypes())

    if (!noRemap) {
      flog.info("%sRemapping final sequences", indent(), name = "info")
      lapply(x$getHapTypes(), function(h)
        refineAlignment(x, h, report = TRUE, createIgv = createIgv, indent = indent))
    }
  }

  x$setReportStatus(TRUE)
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
  map <- match.arg(tolower(map), c("mapfinal", "mapiter", noRemap = FALSE))
  ending <- ifelse(length(x$getHapTypes()) == 2,
                   "psa",
                   ifelse(length(x$getHapTypes()) == 1,
                          "fa",
                          "msa"))
  pairfileUnchecked <- dot(c(map, "aln", x$getLrdType(), x$getLrdMapper(), "unchecked", ending))
  pairfileChecked   <- dot(c(map, "aln", x$getLrdType(), x$getLrdMapper(), "checked", ending))
  pairfileChecked   <- normalizePath(x$absPath(file.path("report",
                                                pairfileChecked)),
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
    file <- dot(c(names(seqs[sq]), x$getLrdType(), x$getLrdMapper(), "fa"))
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
  # flog.info(x$getSampleDetails(), name = "info")
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
                   ifelse(length(x$getHapTypes()) == 1,
                          "fa",
                          "msa"))
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

.checkReportStatus <- function(x){
  if (x$getReportStatus()) {
    currentCall <- strsplit1(strsplit1(deparse(sys.call()), "\\$")[2], "\\(")[1]
    flog.info("%s: Reporting already done! Nothing to do." %<<%
                " Exit safely for downstream analysis",
              currentCall, name = "info")
  }
  return(x$getReportStatus())
}

#' Refine the mapping of an allele by remapping to the checked consensus
#'
#' @param x  A \code{\link[=DR2S_]{DR2S}} object.
#' @param hptype The allele to refine
#' @param report Should the function read a reference from the report folder?
#' Only internal parameter.
#' @return Returns the updated \code{\link[=DR2S_]{DR2S}} object.
#' @return Creates an executable bash file for inspecting the mapping with IGV
#' @family DR2S mapper functions
#' @export
refineAlignment <- function(x, hptype, report = FALSE, createIgv = TRUE, ...) {

  indent <- list(...)$indent %||% indentation()
  opts <- list(...)$opts

  ## Overide default arguments
  args <- x$getOpts("refineMapping")
  if (!is.null(args)) {
    env  <- environment()
    list2env(args, envir = env)
  }

  ## stop if no shortreads provided
  if (!x$hasShortreads()) {
    flog.warn("%sCannot refine mapping. No shortreads provided", indent(), name = "info")
    return(invisible(x))
  }

  if (is.null(x$consensus$refine)) {
    x$consensus$refine <- list(
      mapgroup = list(),
      bamfile = list(),
      consensus = list(),
      ref = list(),
      igv = list()
    )
  }
  reftag <- "refine"
  outdir <- .dirCreateIfNotExists(x$absPath(reftag))
  if (length(x$getHapTypes()) == 1) {
    readpathLR  <- x$absPath(eadpath(x$mapInit))
    if (x$getFilterScores()) {
      readpathSR <- x$absPath(x$mapInit$SR1$reads)
    } else {
      readpathSR <- x$getShortreads()
    }
  } else {
    readpathLR <- x$absPath(readpath(x$mapFinal$LR[[hptype]]))
    readpathSR <- x$absPath(readpath(x$mapFinal$SR[[hptype]]))
  }
  refpath <- ifelse(report, {
    seq <- x$consensus$noAmbig[[hptype]]
    file <- dot(c(names(seq), x$getLrdType(), x$getLrdMapper(), "polished", "fa"))
    names(seq) <- names(seq) %<<% " LOCUS=" %<<%
      x$getLocus() %<<% ";REF=" %<<%  x$getReference()
    seqPath <- x$absPath(file.path("refine", file))
    Biostrings::writeXStringSet(seq, filepath = seqPath, format = "fasta" )
    stats::setNames(seqPath, hptype)
  }, .getUpdatedSeqs(x, hptype))

  x$consensus$refine$ref[[hptype]] <- x$relPath(refpath)
  names(refpath) <- x$relPath(refpath)

  ## Remap long reads to the same reference sequences as short reads
  flog.info("%sRefine mapping for haplotype <%s>", indent(), hptype, name = "info")
  mapgroupLR <- "LR" %<<% hptype
  maptagLR <- dot(c("refine", mapgroupLR, x$getLrdType(), x$getLrdMapper()))
  pileup <- mapReads(
    mapfun = x$getLrdMapFun(), maplabel = reftag, reffile = refpath,
    refname = mapgroupLR, readfile = readpathLR, readtype = x$getLrdType(),
    opts = opts, outdir = outdir, clean = TRUE, includeDeletions = FALSE,
    includeInsertions = FALSE, indent = incr(indent))

  x$consensus$refine$bamfile[[mapgroupLR]] = x$relPath(bampath(pileup))

  if (createIgv)
    x$consensus$refine$igv[[mapgroupLR]] <- createIgvJsFiles(
      refpath(pileup), bampath(pileup), x$getOutdir(), sampleSize = 100,
      fragmentReads = TRUE)

  ## Map short reads
  if (!is.null(unlist(readpathSR))) {
    mapgroupSR <- "SR" %<<% hptype
    maptagSR <- dot(c("refine", mapgroupSR, x$getSrdType(), x$getSrdMapper()))
    readfiles <- readpathSR
    ## Mapper
    pileup <- mapReads(
      mapfun = x$getSrdMapFun(), maplabel = reftag, reffile = refpath,
      refname = mapgroupSR, readfile = readpathSR, readtype = x$getSrdType(),
      opts = opts, outdir = outdir, clean = TRUE, includeDeletions = FALSE,
      includeInsertions = FALSE, clip = TRUE, indent = incr(indent))

    x$consensus$refine$bamfile[[mapgroupSR]] = x$relPath(bampath(pileup))
    x$consensus$refine$igv[[mapgroupSR]] <- createIgvJsFiles(
      refpath(pileup), bampath(pileup), x$getOutdir(), sampleSize = 100)
  }

  # calc new consensus
  consName <- "refine" %<<% hptype
  cseq <- conseq(consmat(pileup), name = consName, type = "ambig",
                 threshold = x$getThreshold(), suppressAllGaps = FALSE)
  x$consensus$refine$consensus[[hptype]] <- cseq
  createIgvConfigs(x, map = "refine", open = FALSE)
  invisible(x)
}


# Helpers -----------------------------------------------------------------

.getUpdatedSeqs <- function(x, hptype){
  ending <- ifelse(length(x$getHapTypes()) == 1, "fa",
                   ifelse(length(x$getHapTypes()) == 2, "psa",
                          "msa"))
  pairfileChecked <- dot(c("mapfinal.aln", x$getLrdType(), x$getLrdMapper(), "checked", ending))
  pairfileChecked <- normalizePath(x$absPath(file.path("report", pairfileChecked)),
                                   mustWork = FALSE)
  rs <- readPairFile(pairfileChecked)
  seqs <- Biostrings::DNAStringSet(
    lapply(rs, function(s) Biostrings::DNAString(stripIndel(s))))

  ## Export FASTA
  seq <- seqs["hap" %<<% hptype]
  file <- dot(c(names(seq), x$getLrdType(), x$getLrdMapper(), "refined", "fa"))
  names(seq) <- names(seq) %<<% " LOCUS=" %<<% x$getLocus() %<<% ";REF=" %<<% x$getReference()
  seqPath <- x$absPath(file.path("refine", file))
  Biostrings::writeXStringSet(seq, filepath = seqPath, format = "fasta")
  seqPath
}

readPairFile <- function(pairfile) {
  if (endsWith(pairfile, "psa")) {
    rs <- readLines(pairfile)
    ## This assignment relies on the premise that hapA is always used!
    ## Usually this should be true, bcs it is the cluster with the most reads
    rsA <- rs[grepl("^hapA", rs)]
    rsB <- rs[grepl("^hap[B-Z]", rs)]
    hap <- c(
      Biostrings::DNAStringSet(.collapsePairLines_(rsA)),
      Biostrings::DNAStringSet(.collapsePairLines_(rsB))
    )

    names(hap) <- c("hapA", strsplit1(rsB, "\\s")[1])
  } else if (endsWith(pairfile, "msa")) {
    hap <- readMSA(pairfile)
  } else if (endsWith(pairfile, "fa")) {
    hap <- Biostrings::readDNAStringSet(pairfile)
    names(hap) <- "hapA"
  }

  ## Check for ambiguous bases
  seqLetters <- Biostrings::uniqueLetters(hap)
  ambigLetters <- seqLetters[which(!seqLetters %in% VALID_DNA(include = "del"))]
  if (!length(ambigLetters) == 0) {
    ambigPositions <- stats::setNames(lapply(ambigLetters, function(x, hap)
      unlist(Biostrings::vmatchPattern(x, hap)), hap = hap), ambigLetters)
    msg <- vapply(seq_along(ambigPositions), function(x, ambigPositions)
      extractAmbigLetters(ambigPositions, names(ambigPositions)[x]),
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
