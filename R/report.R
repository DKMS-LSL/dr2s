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
  aln_file <-  paste(map, "aln", x$getLrdType(), x$getMapper(), "unchecked",
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
    hap_file <- paste(map, hptype, x$getLrdType(), x$getMapper(),
                      "unchecked.fa", sep = ".")
    Biostrings::writeXStringSet(
      haps[[hptype]],
      filepath = file.path(outdir, hap_file),
      format = "fasta"
    )
  }

  ## Write Pairwise Alignment
  if (length(x$getHapTypes()) == 2){
    pair_file <- paste(map, "aln", x$getLrdType(), x$getMapper(), "unchecked",
                       "psa", sep = ".")
    aln <- Biostrings::pairwiseAlignment(pattern = haps[[x$getHapTypes()[[1]]]],
                                         subject = haps[[x$getHapTypes()[[2]]]],
                                         type = "global")
    Biostrings::writePairwiseAlignments(aln, file.path(outdir, pair_file),
                                        block.width = block_width)

  } else if (length(x$getHapTypes()) > 2){
    aln_file <- paste(map, "aln", x$getLrdType(), x$getMapper(), "unchecked",
                      "msa", sep = ".")
    aln <- DECIPHER::AlignSeqs(unlist(Biostrings::DNAStringSetList(haps)),
                               verbose = FALSE)
    Biostrings::write.phylip(Biostrings::DNAMultipleAlignment(aln),
                             file.path(outdir, aln_file))
  }

  if (map == "mapfinal") {
    ## Report problematic Variants
    probvar_file <- paste("problems", x$getLrdType(), x$getMapper(), "tsv",
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
#x <- dpb1_3.r


report_checked_consensus <- function(x, which = "mapFinal") {
  map <- match.arg(tolower(which), c("mapfinal", "mapiter"))
  ending <- ifelse(length(x$getHapTypes()) == 2, "psa", "msa")
  pairfile_unchecked <- paste(map, "aln", x$getLrdType(), x$getMapper(),
                              "unchecked", ending, sep = ".")
  pairfile_checked   <- paste(map, "aln", x$getLrdType(), x$getMapper(),
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
  hapA <- Biostrings::DNAStringSet(gsub("-", "", rs["hapA"]))
  hapB <- Biostrings::DNAStringSet(gsub("-", "", rs["hapB"]))

  ## Alignment
  ref <- x$getRefSeq()
  names(ref) <- strsplitN(names(ref), "~", 1, fixed = TRUE)

  seqs <- c(ref, hapA, hapB)
  aln_file <-  paste("aln", x$getLrdType(), x$getMapper(), "html", sep = ".")
  browse_align(seqs, file = file.path(outdir, aln_file), openURL = FALSE)

  ## Export FASTA
  hapA_file <- paste("A", x$getLrdType(), x$getMapper(), "fa", sep = ".")
  Biostrings::writeXStringSet(
    hapA,
    filepath = file.path(outdir, hapA_file),
    format = "fasta"
  )

  hapB_file <- paste("B", x$getLrdType(), x$getMapper(), "fa", sep = ".")
  Biostrings::writeXStringSet(
    hapB,
    filepath = file.path(outdir, hapB_file),
    format = "fasta"
  )
}

#' @export
check_alignment_file <- function(x, which = "mapFinal", where = 0,
                                 editor = "xdg-open") {
  which <- match.arg(tolower(which), c("mapfinal", "mapiter"))
  ending <- ifelse(length(x$getHapTypes()) == 2, "psa", "msa")
  pairfile_unchecked <- paste(which, "aln", x$getLrdType(), x$getMapper(),
                              "unchecked", ending, sep = ".")
  pairfile_unchecked <- normalizePath(file.path(x$getOutdir(), "report",
                                                pairfile_unchecked),
                                      mustWork = FALSE)
  assertthat::assert_that(
    file.exists(pairfile_unchecked),
    assertthat::is.readable(pairfile_unchecked)
  )
  pairfile_checked <- paste(which, "aln", x$getLrdType(), x$getMapper(),
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
    flog.info("%s: Reporting already done! Nothing to do.
              Exit safely for downstream analysis ...",
              currentCall, name = "info")
  }
  return(x$getReportStatus())
}

# Helpers -----------------------------------------------------------------

readPairFile <- function(pairfile) {

  if (endsWith(pairfile, "psa")){
    rs <- readLines(pairfile)
    rsA <- rs[grepl("^hapA", rs)]
    rsB <- rs[grepl("^hapB", rs)]
    hap <- c(
      Biostrings::DNAStringSet(collapse_pair_lines_(rsA)),
      Biostrings::DNAStringSet(collapse_pair_lines_(rsB))
    )
    names(hap) <- c("hapA", "hapB")
  } else if (endsWith(pairfile, "msa")){
    rs <- Biostrings::readDNAMultipleAlignment(pairfile, format = "phylip")
    hap <- Biostrings::DNAStringSet(rs)
  }
  hap
}

collapse_pair_lines_ <- function(x) {
  paste0(vapply(strsplit(x, split = "\\s+"), `[[`, 3, FUN.VALUE = ""),
         collapse = "")
}
