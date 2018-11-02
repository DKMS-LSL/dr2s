#' @export
report.DR2S <- function(x, which, blockWidth = 80, noRemap = FALSE, createIgv = TRUE, ...) {
  flog.info("# report", name = "info")

  ## Collect start time for report runstats
  start.time <- Sys.time()

  ## Check if reporting is finished and exit safely for downstream analysis
  if (x$getReportStatus()) {
    currentCall <- strsplit1(deparse(sys.call()), "\\.")[1]
    flog.info(strwrap("%s: Reporting already done! Nothing to do.
                      Exit safely for downstream analysis."),
              currentCall, name = "info")
    return(invisible(x))
  }

  args <- x$getOpts("report")
  if (!is.null(args)) {
    env  <- environment()
    list2env(args, envir = env)
  }

  outdir <- .dirCreateIfNotExists(x$absPath("report"))
  if (missing(which)) {
    ## if `which` is unspecified choose `mapFinal` if available,
    ## otherwise try `mapIter`
    if (is(x$mapFinal, "mapFinal")) {
      .reportMap_(x, map = "mapFinal", outdir = outdir, blockWidth = blockWidth,
                  noRemap = noRemap, createIgv = createIgv, ...)
    } else if (is(x$mapIter$`0`$A, "mapIter")) {
      .reportMap_(x, map = "mapIter", outdir = outdir, blockWidth = blockWidth, ...)
    } else
      stop("Nothing to report")
  } else {
    which <- match.arg(tolower(which), c("mapFinal", "mapIter"))
    .reportMap_(x, which, outdir, blockWidth = blockWidth, ...)
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
  haps <- x$getLatestRef()
  ## Initiate indenter
  indent <- indentation(1)
  ## Write html alignment file
  alnFile <- paste(map, "aln", x$getLrdType(), x$getLrdMapper(), "unchecked", sep = ".")
  # get all seqs as stringset
  seqs <- unlist(Biostrings::BStringSetList(haps))
  names(seqs) <- names(haps)
  seqs <- c(ref, seqs)
  .browseAlign(seqs, file = file.path(outdir, alnFile), openURL = FALSE)

  ## Write consensus FASTA files
  # hptype = "A"
  for (hptype in x$getHapTypes()) {
    hapFile <- paste(map, hptype, x$getLrdType(), x$getLrdMapper(), "unchecked.fa", sep = ".")
    seq <- haps[[hptype]]
    seqname <- paste(x$getSampleId(), hptype, sep = "_")
    sampleDetails <- x$getSampleDetails()
    names(seq) <-  paste(seqname,
                         paste("haplotype=" %<<% litQuote(hptype),
                               sampleDetails,
                               "date=" %<<% litQuote(Sys.Date()),
                               "status=" %<<% litQuote("unchecked"),
                               sep = ";"))

    Biostrings::writeXStringSet(
      seq,
      filepath = file.path(outdir, hapFile),
      format = "fasta"
    )
  }

  ## Write Pairwise or Multiple Alignment
  if (length(x$getHapTypes()) == 1) {
    alnFile <- paste(map, "aln", x$getLrdType(), x$getLrdMapper(), "unchecked", "fa", sep = ".")
    Biostrings::writeXStringSet(haps[[1]], filepath = file.path(outdir,alnFile), format = "fasta")
  } else if (length(x$getHapTypes()) == 2) {
    pairFile <- paste(map, "aln", x$getLrdType(), x$getLrdMapper(), "unchecked", "psa", sep = ".")
    aln <- Biostrings::pairwiseAlignment(pattern = haps[[x$getHapTypes()[[1]]]],
                                         subject = haps[[x$getHapTypes()[[2]]]],
                                         type = "global")
    Biostrings::writePairwiseAlignments(aln, file.path(outdir, pairFile), block.width = blockWidth)
  } else if (length(x$getHapTypes()) > 2) {
    alnFile <- paste(map, "aln", x$getLrdType(), x$getLrdMapper(), "unchecked", "msa", sep = ".")
    aln <- DECIPHER::AlignSeqs(Biostrings::DNAStringSet(
      foreach(h = haps, .combine = c) %do% h
    ), verbose = FALSE)
    writeMSA(aln, file = file.path(outdir, alnFile))
    # Biostrings::write.phylip(Biostrings::DNAMultipleAlignment(aln),
    # file.path(outdir, alnFile))
  }

  if (map == "mapfinal") {
    ## Report problematic Variants
    probvarFile <- paste("problems", x$getLrdType(), x$getLrdMapper(), "tsv", sep = ".")
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
#' @param which Which mapping do we want to update.
#' @return Side effect
#' @details
#' \code{reportCheckedConsensus} will create a \code{checked} directory
#' containing the consensus sequences
#' \code{\{A,B\}.\{readtype\}.\{mapper\}.checked.fa}
#' and a html alignment file \code{aln.\{readtype\}.\{mapper\}.checked.html}.
#' @family DR2S mapper functions
#' @export
reportCheckedConsensus <- function(x, which = "mapFinal") {
  map <- match.arg(tolower(which), c("mapfinal", "mapiter", noRemap = FALSE))
  ending <- ifelse(length(x$getHapTypes()) == 2,
                   "psa",
                   ifelse(length(x$getHapTypes()) == 1,
                          "fa",
                          "msa"))
  pairfileUnchecked <- paste(map, "aln", x$getLrdType(), x$getLrdMapper(),
                              "unchecked", ending, sep = ".")
  pairfileChecked   <- paste(map, "aln", x$getLrdType(), x$getLrdMapper(),
                              "checked", ending, sep = ".")
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
    lapply(rs, function(s) Biostrings::DNAString(gsub("-", "", s))))

  ## Alignment
  ref <- x$getRefSeq()
  names(ref) <- strsplitN(names(ref), "~", 1, fixed = TRUE)

  seqsAll <- c(ref, seqs)
  alnFile <-  paste("aln", x$getLrdType(), x$getLrdMapper(), "html", sep = ".")
  .browseAlign(seqsAll, file = file.path(outdir,alnFile), openURL = FALSE)

  ## Export FASTA
  files <- vapply(seq_along(seqs), function(sq, seqs, x) {
    file <- paste(names(seqs[sq]), x$getLrdType(), x$getLrdMapper(), "fa",
                  sep = ".")
    seq <- seqs[sq]
    seqname <- paste(x$getSampleId(), sub("^hap", "", names(seq)), sep = "_")
    sampleDetails <- x$getSampleDetails()
    names(seq) <-  paste(seqname,
                         paste("haplotype=" %<<%
                                 litQuote(sub("^hap", "", names(seq))),
                               sampleDetails,
                               "date=" %<<%
                                 litQuote(Sys.Date()),
                               "status=" %<<%
                                 litQuote("checked"),
                               sep = ";"))

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
#' @param which Which stage of mapping to use. Either "mapFinal" (default) or
#' "mapIter".
#' @param where Where to focus with the editor.
#' @param editor Which editor to use. Either "xdg-open" (default) for standard
#' system editor, "subl" for sublime, "gvim" or "gedit"
#' @param openEditor should the editor be opened or just create the files?
#' @return The path to the newliy created alignment file to check the polished
#'   consensus sequences.
#' @export
checkAlignmentFile <- function(x, which = "mapFinal", where = 0,
                                 editor = "xdg-open", openEditor = TRUE) {
  which <- match.arg(tolower(which), c("mapfinal", "mapiter"))
  ending <- ifelse(length(x$getHapTypes()) == 2,
                   "psa",
                   ifelse(length(x$getHapTypes()) == 1,
                          "fa",
                          "msa"))
  pairfileUnchecked <- paste(which, "aln", x$getLrdType(), x$getLrdMapper(),
                              "unchecked", ending, sep = ".")
  pairfileUnchecked <- normalizePath(x$absPath(file.path("report",
                                                pairfileUnchecked)),
                                      mustWork = FALSE)
  assert_that(
    file.exists(pairfileUnchecked),
    is.readable(pairfileUnchecked)
  )
  pairfileChecked <- paste(which, "aln", x$getLrdType(), x$getLrdMapper(),
                            "checked", ending, sep = ".")
  pairfileChecked <- normalizePath(x$absPath(file.path("report",
                                              pairfileChecked)),
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
                "Exit safely for downstream analysis ...",
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
    readpathLR  <- "/" %<<% x$mapInit$reads
    if (x$getFilterScores()) {
      readpathSR <- x$absPath(x$mapInit$SR1$reads)
    } else {
      readpathSR <- x$getShortreads()
    }
  } else {
    readpathLR <- x$absPath(x$mapFinal$lreads[hptype])
    readpathSR <- x$absPath(x$mapFinal$sreads[[hptype]])
  }
  refpath <- ifelse(report, {
    seq <- x$consensus$noAmbig[[hptype]]
    file <- paste(names(seq), x$getLrdType(), x$getLrdMapper(), "polished", "fa", sep = ".")
    names(seq) <- names(seq) %<<% " LOCUS=" %<<%
      x$getLocus() %<<% ";REF=" %<<%  x$getReference()
    seqPath <- x$absPath(file.path("refine", file))
    Biostrings::writeXStringSet(seq, filepath = seqPath, format = "fasta" )
    set_names(seqPath, hptype)
  },
  set_names(.getUpdatedSeqs(x, hptype)))

  x$consensus$refine$ref[[hptype]] <- x$relPath(refpath)
  names(refpath) <- hptype

  ## Remap long reads to the same reference sequences as short reads
  flog.info("%sRefine mapping for haplotype <%s>", indent(), hptype, name = "info" )
  mapgroupLR <- "LR" %<<% hptype
  maptagLR <- paste("refine", mapgroupLR, x$getLrdType(), x$getLrdMapper(), sep = ".")
  pileup <- mapReads(
    mapFun = x$getLrMapFun(), maptag = maptagLR, reffile = refpath,
    readfile = readpathLR, allele = mapgroupLR, readtype = x$getLrdType(),
    outdir = outdir, force = TRUE, clean = TRUE, includeDeletions = FALSE,
    includeInsertions = FALSE, refname = hptype)

  x$consensus$refine$bamfile[[mapgroupLR]] = x$relPath(path(pileup))

  if (createIgv)
    x$consensus$refine$igv[[mapgroupLR]] <- createIgvJsFiles(
      refpath(pileup), path(pileup), x$getOutdir(), sampleSize = 100,
      fragmentReads = TRUE)

  ## Map short reads
  if (!is.null(unlist(readpathSR))) {
    mapgroupSR <- "SR" %<<% hptype
    maptagSR <- paste("refine", mapgroupSR, x$getSrdType(), x$getSrdMapper(), sep = ".")
    readfiles <- unlist(readpathSR)

    ## Mapper
    pileup <- mapReads(
      mapFun = x$getSrMapFun(), maptag = maptagSR, reffile = refpath,
      readfile = readfiles, allele = mapgroupSR, readtype = x$getSrdType(),
      outdir = outdir, force = TRUE, clean = TRUE, includeDeletions = FALSE,
      includeInsertions = FALSE, refname = hptype, clip = TRUE)

    x$consensus$refine$bamfile[[mapgroupSR]] = x$relPath(path(pileup))
    x$consensus$refine$igv[[mapgroupSR]] <- createIgvJsFiles(
      refpath(pileup), path(pileup), x$getOutdir(), sampleSize = 100)
  }

  # calc new consensus
  consName <- "refine" %<<% hptype
  cseq <- conseq(consmat(pileup), name = consName, type = "ambig",
                 threshold = x$getThreshold(), suppressAllGaps = FALSE)
  x$consensus$refine$consensus[[hptype]] <- cseq
  x$cache()
  createIgvConfigs(x, map = "refine", open = FALSE)
  invisible(x)
}


# Helpers -----------------------------------------------------------------

.getUpdatedSeqs <- function(x, hptype){
  ending <- ifelse(length(x$getHapTypes()) == 1, "fa",
                   ifelse(length(x$getHapTypes()) == 2, "psa",
                          "msa"))
  pairfileChecked <- paste("mapfinal.aln", x$getLrdType(), x$getLrdMapper(),
                              "checked", ending, sep = ".")
  pairfileChecked <- normalizePath(x$absPath(file.path("report",
                                                pairfileChecked)),
                                      mustWork = FALSE)
  rs <- readPairFile(pairfileChecked)
  seqs <- Biostrings::DNAStringSet(
    lapply(rs, function(s) Biostrings::DNAString(gsub("-", "", s))))

  ## Export FASTA
  seq <- seqs["hap" %<<% hptype]
  file <- paste(names(seq),
                x$getLrdType(),
                x$getLrdMapper(),
                "refined",
                "fa",
                sep = ".")
  names(seq) <- names(seq) %<<% " LOCUS=" %<<% x$getLocus() %<<% ";REF=" %<<%
                          x$getReference()
  seqPath <- x$absPath(file.path("refine", file))
  Biostrings::writeXStringSet(
    seq,
    filepath = seqPath,
    format = "fasta"
  )
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
    ambigPositions <- set_names(lapply(ambigLetters, function(x, hap)
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
