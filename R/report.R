#' @export
report.DR2S <- function(x, which, block_width = 80, ...) {

  args <- x$getOpts("report")
  if (!is.null(args)) {
    env  <- environment()
    list2env(args, envir = env)
  }

  outdir <- dir_create_if_not_exists(file.path(x$getOutdir(), "report"))
  if (missing(which)) {
    ## if `which` is unspecified choose `map3` if available,
    ## otherwise try `map2`, then try `map1`
    if (is(x$consensus, "ConsList")) {
      report_map_(x, map = "map3", outdir = outdir, block_width = block_width, ...)
    } else if (all(sapply(x$map2, is, "map2"))) {
      report_map_(x, map = "map2", outdir = outdir, block_width = block_width, ...)
    } else if (all(sapply(x$map1, is, "map1"))) {
      report_map_(x, map = "map1", outdir = outdir, block_width = block_width, ...)
    } else
      stop("Nothing to report")
  } else {
    which <- match.arg(tolower(which), c("map3", "map2", "map1"))
    report_map_(x, which, outdir, block_width = block_width, ...)
  }
}

report_map_ <- function(x, map, outdir, block_width, ...) {
  map <- match.arg(tolower(map), c("map3", "map2", "map1"))
  ref <-  Biostrings::BStringSet(x$getRefSeq())
  names(ref) <- strsplitN(names(ref), "~", 1, fixed = TRUE)
  alt <-  Biostrings::BStringSet(x$getAltSeq())
  if (length(alt) > 0) {
    names(alt) <- strsplitN(names(alt), "~", 1, fixed = TRUE)
  }

  addins <- list(...)$addins
  if (!is.null(addins)) {
    addins <- Biostrings::readBStringSet(addins)
    names(addins) <- strsplitN(names(addins), "~", 1, fixed = TRUE)
  }

  if (map == "map1") {
    if (x$hasMap1Alternate()) {
      hapA <- x$map1$A$conseq$merged
      hapB <- x$map1$B$conseq$merged
    } else {
      hapA <- x$map1$A$conseq$reference
      hapB <- x$map1$B$conseq$reference
    }
  } else if (map == "map2") {
    hapA <- x$map2$A$conseq
    hapB <- x$map2$B$conseq
  } else if (map == "map3") {
    hapA <- x$consensus$seq["HapA"]
    hapB <- x$consensus$seq["HapB"]
  }

  ## Write html alignment file
  aln_file <-  paste(map, "aln", x$getLrdType(), x$getMapper(), "unchecked.html", sep = ".")
  browse_align(c(ref, hapA, alt, hapB, addins), file = file.path(outdir, aln_file), openURL = FALSE)

  ## Write consensus FASTA files
  hapA_file <- paste(map, "A", x$getLrdType(), x$getMapper(), "unchecked.fa", sep = ".")
  Biostrings::writeXStringSet(
    hapA,
    filepath = file.path(outdir, hapA_file),
    format = "fasta"
  )

  hapB_file <- paste(map, "B", x$getLrdType(), x$getMapper(), "unchecked.fa", sep = ".")
  Biostrings::writeXStringSet(
    hapB,
    filepath = file.path(outdir, hapB_file),
    format = "fasta"
  )

  ## Write Pairwise Alignment
  pair_file <- paste(map, "aln", x$getLrdType(), x$getMapper(), "unchecked.pair", sep = ".")
  aln <- Biostrings::pairwiseAlignment(pattern = hapA, subject = hapB, type = "global")
  Biostrings::writePairwiseAlignments(aln, file.path(outdir, pair_file), block.width = block_width)

  if (map == "map3") {
    ## Report problematic Variants
    probvar_file <- paste("problems", x$getLrdType(), x$getMapper(), "tsv", sep = ".")
    vars <- x$consensus$problematic_variants %>%
      dplyr::arrange(pos, haplotype)
    readr::write_tsv(vars, path = file.path(outdir, probvar_file), append = FALSE, col_names = TRUE)
  }

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
report_checked_consensus <- function(x, which = "map3") {
  map <- match.arg(tolower(which), c("map3", "map2", "map1"))
  pairfile_unchecked <- paste(map, "aln", x$getLrdType(), x$getMapper(), "unchecked.pair", sep = ".")
  pairfile_checked   <- paste(map, "aln", x$getLrdType(), x$getMapper(), "checked.pair", sep = ".")
  pairfile_checked   <- normalizePath(file.path(x$getOutdir(), "report", pairfile_checked), mustWork = FALSE)

  if (!file.exists(pairfile_checked)) {
    msg <- sprintf("The file '%s' must be saved as '%s' before invoking this function",
                   pairfile_unchecked, pairfile_checked)
    stop(msg, call. = FALSE)
  }

  outdir <- dkms::dir_create_if_not_exists(file.path(x$getOutdir(), "checked"))
  rs <- readPairFile(pairfile_checked)
  hapA <- Biostrings::DNAStringSet(gsub("-", "", rs["HapA"]))
  hapB <- Biostrings::DNAStringSet(gsub("-", "", rs["HapB"]))

  ## Alignment
  ref <- x$getRefSeq()
  names(ref) <- strsplitN(names(ref), "~", 1, fixed = TRUE)
  alt <- x$getAltSeq()
  if (!is.null(alt)) {
    names(alt) <- strsplitN(names(alt), "~", 1, fixed = TRUE)
  }
  seqs <- c(ref, Biostrings::BStringSet(hapA), alt, Biostrings::BStringSet(hapB))
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
check_alignment_file <- function(x, which = "map3", where = 0) {
  which <- match.arg(tolower(which), c("map3", "map2", "map1"))
  pairfile_unchecked <- paste(which, "aln", x$getLrdType(), x$getMapper(), "unchecked", "pair", sep = ".")
  pairfile_unchecked <- normalizePath(file.path(x$getOutdir(), "report", pairfile_unchecked), mustWork = FALSE)
  assertthat::assert_that(
    file.exists(pairfile_unchecked),
    assertthat::is.readable(pairfile_unchecked)
  )
  pairfile_checked <- paste(which, "aln", x$getLrdType(), x$getMapper(), "checked", "pair", sep = ".")
  pairfile_checked <- normalizePath(file.path(x$getOutdir(), "report", pairfile_checked), mustWork = FALSE)
  if (!file.exists(pairfile_checked)) {
    file.copy(pairfile_unchecked, pairfile_checked, overwrite = FALSE)
  }
  where <- 23 + 4*floor((where - 1)/80)
  subl(pairfile_checked, where)
}


# Helpers -----------------------------------------------------------------


readPairFile <- function(pairfile) {
  rs <- readLines(pairfile)
  rsA <- rs[grepl("^HapA", rs)]
  rsB <- rs[grepl("^HapB", rs)]
  hap <- c(
    Biostrings::DNAStringSet(collapse_pair_lines_(rsA)),
    Biostrings::DNAStringSet(collapse_pair_lines_(rsB))
  )
  names(hap) <- c("HapA", "HapB")
  hap
}

collapse_pair_lines_ <- function(x) {
  paste0(vapply(strsplit(x, split = "\\s+"), `[[`, 3, FUN.VALUE = ""), collapse = "")
}


# experimental ------------------------------------------------------------


#' @export
DR2Sseqs <- function(ncores, verbose = FALSE) {
  stopifnot(
    requireNamespace("parallel", quietly = TRUE),
    requireNamespace("doParallel", quietly = TRUE)
  )
  outf <- if (verbose)  "" else "/dev/null"
  if (ncores == 1) {
    foreach::registerDoSEQ()
  } else {
    cl0 <- parallel::makePSOCKcluster(ncores, outfile = outf)
    parallel::clusterEvalQ(cl0, {
      library("DR2S", quietly = TRUE)
      library("Biostrings", quietly = TRUE)
      invisible(NULL)
    })
    doParallel::registerDoParallel(cl0)
  }
  pkgs <- c("DR2S", "Biostrings")
  list(
    fetch = function(path, locus = NULL, verbose = FALSE) {
      locus <- if (is.null(locus)) {
        ".+"
      } else normalise_locus(locus)
      pattern <- paste0("^", locus, "\\.(nanopore||pacbio)\\.(alt|ref|cons)(\\.(alt|multialign))?$")
      all_dirs <- list.dirs(path, recursive = TRUE, full.names = TRUE)
      dirs <- normalizePath(all_dirs[grep(pattern, basename(all_dirs))], mustWork = TRUE)
      # d = dirs[2]
      rs <- foreach(d = dirs, .combine = "c", .errorhandling = "remove",
                    .packages = pkgs, .verbose = verbose) %dopar% {
        x <- DR2S::read_dr2s(d)
        refs <- DR2S:::extract_references(x)
        cons <- DR2S:::extract_consensus(x)
        Biostrings::DNAStringSet(c(refs[unique(names(refs))], cons[order(names(cons))]))
      }
      rs[unique(names(rs))]
    },
    stopCluster = function() {
      if (ncores > 1) {
        parallel::stopCluster(cl0)
      }
    }
  )
}

extract_references <- function(x) {
  locus <- normalise_locus(x$getLocus())
  sample <- x$getSampleId()
  ref <-  Biostrings::BStringSet(x$getRefSeq())
  refname <- strsplitN(names(ref), "~", 1, fixed = TRUE)
  allele  <- strsplitN(refname, "[ *]", 2)
  if (allele == "consensus") {
    names(ref) <- paste0("seqtype=cns|allele=", allele, "|sample=", sample, "|locus=", locus, "|")
  } else {
    names(ref) <- paste0("seqtype=ref|cnstype=reference|allele=", allele, "|sample=", sample, "|locus=", locus, "|")
  }
  alt <-  Biostrings::BStringSet(x$getAltSeq())
  if (length(alt) > 0) {
    altname <- strsplitN(names(alt), "~", 1, fixed = TRUE)
    allele  <- strsplitN(altname, "[ *]", 2)
    names(alt) <- paste0("seqtype=ref|cnstype=alternate|allele=", allele, "|sample=", sample, "|locus=", locus, "|")
  }
  c(ref, alt)
}

extract_consensus <- function(x) {
  lrdtype  <- x$getLrdType()
  reftype  <- x$getReftype()
  haprefA  <- x$getARefType()
  haprefB  <- x$getBRefType()
  locus    <- DR2S:::normalise_locus(x$getLocus())
  sample   <- x$getSampleId()
  nlreadsA <- length(x$partition$hpl$A)
  nlreadsB <- length(x$partition$hpl$B)
  if (lrdtype == "nanopore") {
    tag <- unique(strsplitN(as.character(x$partition$hpl$A), "_", 1, from = "end"))
    if (all(c("template", "complement") %in% tag)) {
      lrdtype <- paste0(lrdtype, "_1D")
    } else if (tag == "2d") {
      lrdtype <- paste0(lrdtype, "_2D")
    }
  }
  map2_A_cons <- x$map2$A$conseq
  names(map2_A_cons) <- name_map2(names(map2_A_cons), reftype, lrdtype, haprefA, sample, locus)
  map2_B_cons <- x$map2$B$conseq
  names(map2_B_cons) <- name_map2(names(map2_B_cons), reftype, lrdtype, haprefB, sample, locus)
  map3_A_cons <- x$consensus$seq["HapA"]
  names(map3_A_cons) <- name_map3(names(map3_A_cons), reftype, lrdtype, haprefA, sample, locus, nlreadsA)
  map3_B_cons <- x$consensus$seq["HapB"]
  names(map3_B_cons) <- name_map3(names(map3_B_cons), reftype, lrdtype, haprefB, sample, locus, nlreadsB)
  names(c(map2_A_cons, map2_B_cons, map3_A_cons, map3_B_cons))
}

name_map2 <- function(name, reftype, lrdtype, hapref, sample, locus) {
  reftype <- switch(reftype,
                    cons = "cons",
                    ref  = "ref1",
                    alt  = "ref2",
                    ref.alt = "ref.ref"
  )
  parts <- strsplit(name, ".", fixed = TRUE)[[1]][6]
  parts <- c("map", "LR", paste(locus, reftype, substr(parts, 1, 1), sep = "."),
             lrdtype, reftype, hapref, substring(parts, 2), sample, locus)
  names(parts) <- c("seqtype", "cnstype", "haplotype", "lrdtype", "reftype",
                    "hapref", "nlreads", "sample", "locus")
  paste0(paste(paste(names(parts), parts, sep = "="), collapse = "|"), "|")
}

name_map3 <- function(name, reftype, lrdtype, hapref, sample, locus, nlreads) {
  reftype <- switch(reftype,
                    cons = "cons",
                    ref  = "ref1",
                    alt  = "ref2",
                    ref.alt = "ref.ref"
  )
  parts <- c("map", "SR",  paste(locus, reftype, substr(name, 4, 4), sep = "."),
             lrdtype, reftype, hapref, nlreads, sample, locus)
  names(parts) <- c("seqtype", "cnstype", "haplotype", "lrdtype", "reftype",
                    "hapref", "nlreads", "sample", "locus")
  paste0(paste(paste(names(parts), parts, sep = "="), collapse = "|"), "|")
}
