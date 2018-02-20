#' @export
#report(i
# x <-hla.mapFinal
# debug
#x <- dedk.bla
#block_width = 80
report.DR2S <- function(x, which, block_width = 80, ...) {

  ## Check if reporting is already finished and exit safely for downstream analysis
  if (x$getReportStatus()) {
    currentCall <- strsplit(deparse(sys.call()), "\\.")[[1]][1]
    flog.info(strwrap("%s: Reporting already done! Nothing to do.
                      Exit safely for downstream analysis ..."),
              currentCall, name = "info")
    return(invisible(x))
  }

  args <- x$getOpts("report")
  if (!is.null(args)) {
    env  <- environment()
    list2env(args, envir = env)
  }

  outdir <- dir_create_if_not_exists(file.path(x$getOutdir(), "report"))
  if (missing(which)) {
    ## if `which` is unspecified choose `mapFinal` if available,
    ## otherwise try `mapIter`, then try `map1`
    if (is(x$mapFinal, "mapFinal")) {
      report_map_(x, map = "mapFinal", outdir = outdir,
                  block_width = block_width, ...)
    } else if (all(is(object = x$mapIter$`0`$A, class2 = "mapIter"))) {
      report_map_(x, map = "mapIter", outdir = outdir,
                  block_width = block_width, ...)
    } else
      stop("Nothing to report")
  } else {
    which <- match.arg(tolower(which), c("mapFinal", "mapIter"))
    report_map_(x, which, outdir, block_width = block_width, ...)
  }
}
# debug
#map = "mapFinal"
report_map_ <- function(x, map, outdir, block_width, ...) {
  map <- match.arg(tolower(map), c("mapfinal", "mapiter"))
  ref <-  Biostrings::BStringSet(x$getRefSeq())
  names(ref) <- strsplitN(names(ref), "~", 1, fixed = TRUE)
  addins <- list(...)$addins
  if (!is.null(addins)) {
    addins <- Biostrings::readBStringSet(addins)
    names(addins) <- strsplitN(names(addins), "~", 1, fixed = TRUE)
  }
  haps <- x$getLatestRef()

  ## Write html alignment file
  aln_file <-  paste(map, "aln", x$getLrdType(), x$getLrMapper(), "unchecked",
                     sep = ".")
  # get all seqs as stringset
  seqs <- unlist(Biostrings::BStringSetList(haps))
  names(seqs) <- names(haps)
  # if (is.null(addins)){
    # seqs <- c(ref, seqs, alt )
  # } else {
  seqs <- c(ref, seqs, addins)
  # }
  browse_align(seqs, file = file.path(outdir, aln_file), openURL = FALSE)

  ## Write consensus FASTA files

  for (hptype in x$getHapTypes()) {
    hap_file <- paste(map, hptype, x$getLrdType(), x$getLrMapper(),
                      "unchecked.fa", sep = ".")
    seq <- haps[[hptype]]
    names(seq) <- paste(names(seq), " LOCUS=",
                                   x$getLocus(), ";REF=", x$getReference())
    Biostrings::writeXStringSet(
      seq,
      filepath = file.path(outdir, hap_file),
      format = "fasta"
    )
  }

  ## Write Pairwise or Multiple Alignment
  if (length(x$getHapTypes()) == 2){
    pair_file <- paste(map, "aln", x$getLrdType(), x$getLrMapper(), "unchecked",
                       "psa", sep = ".")
    aln <- Biostrings::pairwiseAlignment(pattern = haps[[x$getHapTypes()[[1]]]],
                                         subject = haps[[x$getHapTypes()[[2]]]],
                                         type = "global")
    Biostrings::writePairwiseAlignments(aln, file.path(outdir, pair_file),
                                        block.width = block_width)

  } else if (length(x$getHapTypes()) == 1){
    aln_file <- paste(map, "aln", x$getLrdType(), x$getLrMapper(), "unchecked",
                      "fa", sep = ".")
    Biostrings::writeXStringSet(
      haps[[1]],
      filepath = file.path(outdir, aln_file),
      format = "fasta")
  } else if (length(x$getHapTypes()) > 2){
    aln_file <- paste(map, "aln", x$getLrdType(), x$getLrMapper(), "unchecked",
                      "msa", sep = ".")

    aln <- DECIPHER::AlignSeqs(Biostrings::DNAStringSet(
      foreach(h = haps, .combine=c) %do% h),
                               verbose = FALSE)
    writeMSA(aln , file = file.path(outdir, aln_file))
    # Biostrings::write.phylip(Biostrings::DNAMultipleAlignment(aln),
                             # file.path(outdir, aln_file))
  }

  if (map == "mapfinal") {
    ## Report problematic Variants
    probvar_file <- paste("problems", x$getLrdType(), x$getLrMapper(), "tsv",
                          sep = ".")
    vars <- x$consensus$problematic_variants %>%
      dplyr::arrange(pos, haplotype)
    readr::write_tsv(vars, path = file.path(outdir, probvar_file),
                     append = FALSE, col_names = TRUE)
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
#' \code{report_checked_consensus} will create a \code{checked} directory
#' containing the consensus sequences \code{\{A,B\}.\{readtype\}.\{mapper\}.checked.fa}
#' and a html alignment file \code{aln.\{readtype\}.\{mapper\}.checked.html}.
#' @family DR2S mapper functions
#' @export
#x <- dpb1_2.rep

report_checked_consensus <- function(x, which = "mapFinal") {
  map <- match.arg(tolower(which), c("mapfinal", "mapiter"))
  ending <- ifelse(length(x$getHapTypes()) == 2,
                   "psa",
                   ifelse(length(x$getHapTypes()) == 1,
                          "fa",
                          "msa"))
  pairfile_unchecked <- paste(map, "aln", x$getLrdType(), x$getLrMapper(),
                              "unchecked", ending, sep = ".")
  pairfile_checked   <- paste(map, "aln", x$getLrdType(), x$getLrMapper(),
                              "checked", ending, sep = ".")
  pairfile_checked   <- normalizePath(file.path(x$getOutdir(), "report",
                                                pairfile_checked),
                                      mustWork = FALSE)

  if (!file.exists(pairfile_checked)) {
    msg <- sprintf("The file '%s' must be saved as '%s' before invoking this function",
                   pairfile_unchecked, pairfile_checked)
    stop(msg, call. = FALSE)
  }

  outdir <- dir_create_if_not_exists(file.path(x$getOutdir(), "checked"))
  rs <- readPairFile(pairfile_checked)
  seqs <- Biostrings::DNAStringSet(
    sapply(rs, function(s) Biostrings::DNAString(gsub("-", "", s))))

  ## Alignment
  ref <- x$getRefSeq()
  names(ref) <- strsplitN(names(ref), "~", 1, fixed = TRUE)

  seqsAll <- c(ref, seqs)
  aln_file <-  paste("aln", x$getLrdType(), x$getLrMapper(), "html", sep = ".")
  browse_align(seqsAll, file = file.path(outdir, aln_file), openURL = FALSE)

  ## Export FASTA
  files <- sapply(1:length(seqs), function(sq) {
    file <- paste(names(seqs[sq]), x$getLrdType(), x$getLrMapper(), "fa", sep = ".")
    seq <- seqs[sq]
    names(seq) <- paste0(names(seq), " LOCUS=", x$getLocus(), ";REF=",
                              x$getReference())
    Biostrings::writeXStringSet(
      seq,
      filepath = file.path(outdir, file),
      format = "fasta"
    )
    file
    }
  )
}

#' @export
check_alignment_file <- function(x, which = "mapFinal", where = 0,
                                 editor = "xdg-open") {
  which <- match.arg(tolower(which), c("mapfinal", "mapiter"))
  ending <- ifelse(length(x$getHapTypes()) == 2,
                   "psa",
                   ifelse(length(x$getHapTypes()) == 1,
                          "fa",
                          "msa"))
  pairfile_unchecked <- paste(which, "aln", x$getLrdType(), x$getLrMapper(),
                              "unchecked", ending, sep = ".")
  pairfile_unchecked <- normalizePath(file.path(x$getOutdir(), "report",
                                                pairfile_unchecked),
                                      mustWork = FALSE)
  assertthat::assert_that(
    file.exists(pairfile_unchecked),
    assertthat::is.readable(pairfile_unchecked)
  )
  pairfile_checked <- paste(which, "aln", x$getLrdType(), x$getLrMapper(),
                            "checked", ending, sep = ".")
  pairfile_checked <- normalizePath(file.path(x$getOutdir(), "report",
                                              pairfile_checked),
                                    mustWork = FALSE)
  if (!file.exists(pairfile_checked)) {
    file.copy(pairfile_unchecked, pairfile_checked, overwrite = FALSE)
  }
  where <- 23 + 4*floor((where - 1)/80)
  editor(pairfile_checked, where, use_editor = editor)
}

check_report_status <- function(x){
  if (x$getReportStatus()) {
    currentCall <- strsplit(strsplit(deparse(sys.call()), "\\$")[[1]][2],
                            "\\(")[[1]][1]
    flog.info(paste0("%s: Reporting already done! Nothing to do.",
              "Exit safely for downstream analysis ..."),
              currentCall, name = "info")
  }
  return(x$getReportStatus())
}
#hptype = "B"
refineAlignment <- function(x, hptype){
  ## Overide default arguments
  args <- x$getOpts("refineMapping")
  if (!is.null(args)) {
    env  <- environment()
    list2env(args, envir = env)
  }

  ## stop if no shortreads provided
  if (is.null(x$getConfig("shortreads"))) {
    flog.warn(" Cannot refine mapping. No shortreads provided", name = "info")
    return(invisible(x))
  }


  if (is.null(x$consensus$refine)) {
    x$consensus$refine <- list(
      mapgroup = list(),
      bamfile = list(),
      consensus = list(),
      ref = list()
    )
  }
  reftag    <- "refine"
  outdir    <- dir_create_if_not_exists(file.path(x$getOutdir(), reftag))
  readpathLR  <- x$mapFinal$lreads[hptype]
  refpath   <- .getUpdatedSeqs(x, hptype)
  x$consensus$refine$ref[[hptype]] <- refpath
  names(refpath) <- hptype
  readpathSR <- x$mapFinal$sreads[hptype]


  ## Remap long reads to the same reference sequences as short reads
  flog.info(" Refine mapping for haplotype %s ...", hptype, name = "info" )
  mapgroupLR <- paste0("LR", hptype)
  maptagLR   <- paste("refine", mapgroupLR, x$getLrdType(),
                      x$getLrMapper(), sep = ".")
  ## Mapper
  map_fun <- x$getLrMapFun()
  ## Run mapper
  flog.info("   Mapping long reads against final consensus sequence ...",
            name = "info")
  samfile <- map_fun(
    reffile  = refpath,
    readfile = readpathLR,
    allele   = mapgroupLR,
    readtype = x$getLrdType(),
    refname  = hptype,
#     opts     = list(
#       x = "map-pb -a --cs -c",
# #-c -H -Q -Y",
#
#       k = 10,
#       w = 3,
#       N = 500
#       ),
    force    = TRUE,
    outdir   = outdir
  )

  ## Run bam - sort - index pipeline
  flog.info("  Indexing ...", name = "info")
  bamfile <- bam_sort_index(
    samfile = samfile,
    reffile = refpath,
    force   = TRUE
  )
  x$consensus$refine$bamfile[[mapgroupLR]] = bamfile

  # ## Calculate pileup from graphmap produced SAM file
  # flog.info("  Piling up ...", name = "info")
  # pileup <- Pileup(
  #   bamfile,
  #   x$getThreshold(),
  #   include_deletions = TRUE,
  #   include_insertions = TRUE
  # )
  ## Map short reads
  if (!is.null(unlist(readpathSR))){
    mapgroupSR <- paste0("SR", hptype)
    maptagSR   <- paste("refine", mapgroupSR, x$getLrdType(),
                        x$getSrMapper(), sep = ".")

    readfiles <- unlist(readpathSR)
    ## Mapper
    map_fun <- x$getSrMapFun()

    ## Run mapper
    flog.info("  Mapping short reads against final consensus ...",
              name = "info")
    samfile <- map_fun(
      reffile  = refpath,
      readfile = readfiles,
      allele   = mapgroupSR,
      readtype = x$getSrdType(),
      refname  = hptype,
      # opts     = list(
      #   x      = "sr",
      #   k      = 6,
      #   w = 2,
      #   A = 2,
      #   B = 8,
      #   O = "12,32",
      #   E = "2,1",
      #   r = 50,
      #   p = 0.5,
      #   N = 500
      # ),
      outdir   = outdir,
      force    = TRUE
    )
    ## Run bam - sort - index pipeline
    flog.info("  Indexing ...", name = "info")
    bamfile <- bam_sort_index(
      samfile,
      refpath,
      force = TRUE
    )
    x$consensus$refine$bamfile[[mapgroupSR]] = bamfile

  }
  ## Calculate pileup from graphmap produced SAM file
  flog.info("  Piling up ...", name = "info")
  pileup <- Pileup(
    bamfile,
    x$getThreshold(),
    include_deletions = TRUE,
    include_insertions = TRUE
  )
  pileup <- pileup_include_read_ids(pileup)
  pileup <- pileup_include_insertions(pileup)

  # calc new consensus
  cseq <- conseq(pileup$consmat, paste0("refine", hptype), "ambig",
                 exclude_gaps = FALSE, threshold = x$getThreshold())
  x$consensus$refine$consensus[[hptype]] <- cseq
  invisible(x)
}





# Helpers -----------------------------------------------------------------

.getUpdatedSeqs <- function(x, hptype){
  ending <- ifelse(length(x$getHapTypes()) == 2, "psa", "msa")
  pairfile_checked <- paste("mapfinal.aln", x$getLrdType(), x$getLrMapper(),
                              "checked", ending, sep = ".")
  pairfile_checked <- normalizePath(file.path(x$getOutdir(), "report",
                                                pairfile_checked),
                                      mustWork = FALSE)
  rs <- readPairFile(pairfile_checked)
  seqs <- Biostrings::DNAStringSet(
    sapply(rs, function(s) Biostrings::DNAString(gsub("-", "", s))))

  ## Export FASTA
  seq <- seqs[paste0("hap", hptype)]
  file <- paste(names(seq), x$getLrdType(), x$getLrMapper(), "refined", "fa", sep = ".")
  names(seq) <- paste0(names(seq), " LOCUS=", x$getLocus(), ";REF=",
                          x$getReference())
  seqPath <- file.path(outdir, file)
  Biostrings::writeXStringSet(
    seq,
    filepath = seqPath,
    format = "fasta"
  )
  seqPath
}

readPairFile <- function(pairfile) {
  if (endsWith(pairfile, "psa")){
    rs <- readLines(pairfile)
    ## This assignment relies on the premise that hapA is always used!
    ## Usually this should be true, bcs it is the cluster with the most reads
    rsA <- rs[grepl("^hapA", rs)]
    rsB <- rs[grepl("^hap[B-Z]", rs)]
    hap <- c(
      Biostrings::DNAStringSet(collapse_pair_lines_(rsA)),
      Biostrings::DNAStringSet(collapse_pair_lines_(rsB))
    )
    names(hap) <- c("hapA", "hapB")
  } else if (endsWith(pairfile, "msa")){
    hap <- readMSA(pairfile)
  } else if (endsWith(pairfile, "fa")){
    hap <- Biostrings::readDNAStringSet(pairfile)
    names(hap) <- "hapA"
  }

  ## Check for ambiguous bases
  seqLetters <- Biostrings::uniqueLetters(hap)
  ambigLetters <- seqLetters[which(!seqLetters %in% VALID_DNA(include = "del"))]
  if (!length(ambigLetters) == 0){
    ambigPositions <- sapply(ambigLetters, function(x)
      unlist(Biostrings::vmatchPattern(x, hap)))
    msg <- sapply(1:length(ambigPositions), function(x)
      extractAmbigLetters(ambigPositions, names(ambigPositions)[x]))
    flog.info(msg, name = "info")
    stop("Check reported reference! Ambiguous positions were found")
  }


  hap
}

collapse_pair_lines_ <- function(x) {
  paste0(vapply(strsplit(x, split = "\\s+"), `[[`, 3, FUN.VALUE = ""),
         collapse = "")
}

#' @export
readMSA <- function(file) {
  seqs <- c()
  rows <- scan(file, what = "", sep = "\n", strip.white = TRUE,
               quiet = TRUE, blank.lines.skip = TRUE)

  for (i in seq_len(length(rows))){
    if (!grepl("^[A-Za-z]", rows[i]))
      next
    line <-unlist(strsplit(x = rows[i], "\\s+"))
    name <- line[1]
    if (is.null(seqs) || is.na(seqs[name]))
      seqs[name] <- ""
    seq <- paste0(line[3:(length(line)-1)], collapse = "")
    seqs[name] <- paste0(seqs[name], seq, collapse = "")
  }
  seqs <- Biostrings::DNAStringSet(seqs)
  seqs
}
## HELPER ##
extractAmbigLetters <- function(irange, ambigLetter){
  msg <- sprintf("Found %s; decide for %s at positions: \n", ambigLetter,
                 paste(strsplit(CODE_MAP()[ambigLetter], "")[[1]],
                       collapse = " or "))
  ambigPositionLetter <- unlist(irange[[ambigLetter]])

  msg %<<% paste0(sapply(1:length(ambigPositionLetter), function(x)
    sprintf("%s: %s",
            names(ambigPositionLetter[x]),
            ambigPositionLetter[x])), collapse = "\n") %<<% "\n"
}



## Write a multiple sequence alignment in a phylip like format. Formatted for
## easily accessing positions and differences
#' @export
writeMSA <- function(aln, file="", block.width = 50){
  # ToDo add check for length
  if (!is(aln, "DNAStringSet"))
    stop("'aln' must be a DNAStringSet object and each sequence of same length")
  if (!is.numeric(block.width))
    stop("'block.width' must be a single number")

  get_cons_position <- function(pos) {
   if (length(unique(pos)) == 1) {
     return(".")
   } else if ("-" %in% pos){
     return("+")
   } else {
     return("|")
   }
  }

  aln_mat <-as.matrix(aln)
  aln_stat_line <- Biostrings::BStringSet(paste0(apply(aln_mat, 2, function(x) get_cons_position(x)), collapse = ""))
  alignment <- c(Biostrings::BStringSet(aln), aln_stat_line)

  ## Write Header
  cat("#=======================================\n", file=file)
  cat("#\n", file=file, append = TRUE)
  cat("# Aligned_sequences: ", length(aln)," \n", file=file, append = TRUE)
  sapply(1:length(aln), function(x) cat(sprintf("# %s: %s\n", x, names(aln[x])), file=file, append = TRUE))
  cat("#\n#\n", file=file, append = TRUE)
  cat("#=======================================\n", file=file,append = TRUE)

  # Write sequences
  lstart <- 1L
  lend <- block.width
  alignment_length <- Biostrings::width(alignment[1])
  start_width <- nchar(as.character(1 + alignment_length))
  name_width <- max(20L, nchar(names(alignment)))
  nblock <- alignment_length %/% block.width
  if (alignment_length %% block.width != 0L)
    nblock <- nblock + 1L
  for (i in seq_len(nblock)) {
    to <- i * block.width
    from <- to - block.width + 1L
    if (to > alignment_length)
      to <- alignment_length
    names <- names(alignment)
    strings <- Biostrings::subseq(alignment, from, to)

    ## Split the seq every 10 chars
    sp <- "(.{10})"
    addSp <- "\\1 "
    a <- sapply(1:(length(alignment)-1), function(x) {
      cat(format(names[x], width = name_width), " ",
          format(lstart, justify = "right", width = start_width), " ",
          gsub(sp, addSp, Biostrings::toString(strings[x])), " ",
          format(lend, justify = "right"), "\n",
          sep = "", file = file, append = TRUE)
    })
    cat(format(" ", width = name_width), " ",
          format(" ", justify = "right", width = start_width), " ",
          gsub(sp, addSp, Biostrings::toString(strings[length(alignment)])),
        "\n\n", sep = "", file = file, append = TRUE)

    lstart <- lend + 1
    lend <- lstart + block.width - 1
  }
}
