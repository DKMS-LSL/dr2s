variant_matrix <- function(aln) {
  stopifnot(is(aln, "XStringSet"))
  seqnames <- names(aln)
  mat <- Biostrings::as.matrix(aln)
  mat_ <- factor(mat)
  dim(mat_) <- dim(mat)
  lvls <- levels(mat_)
  nbins <- length(lvls)
  i <- which(apply(mat_, 2, function(x) {
    sum(tabulate(x, nbins) > 0) > 1
  }))
  vmat <- mat_[, i]
  levels(vmat) <- sub("-", ".", levels(vmat))
  dimnames(vmat) <- list(seqnames, i)
  vmat
}

snp_gap_matrix <- function(mat) {
  gap <- which(levels(mat) == ".")
  gap_idx <- apply(mat, 2, function(col) gap %in% col)
  list(
    snp = mat[, !gap_idx],
    gap = mat[, gap_idx]
  )
}

snp_error <- function(mat) {
  ref <- mat[grep("^seqtype=ref", tolower(rownames(mat))), ]
  snp_pos <- colnames(mat)
  mat_ <- sweep(mat, 2, ref, `!=`)
  list(
    err_num = rowSums(mat_),
    err_pos = apply(mat_, 1, function(i) paste0(snp_pos[i], collapse = ":"))
  )
}

indel <- function(x) {
  i <- diff(as.integer(x)) == 1L
  beginnings <- x[!c(FALSE, i)]
  endings <- x[!c(i, FALSE)]
  runs <- ifelse(beginnings == endings, 1, endings - beginnings + 1)
  structure(x, runs = rep.int(seq_along(runs), runs), class = "indel")
}

`[.indel` <- function(x, i, j, drop = FALSE) {
  stopifnot(is.numeric(i))
  as.vector(x)[attr(x, "runs") == i]
}

length.indel <- function(x) {
  length(unique(attr(x, "runs")))
}

iter.indel <- function(x) {
  i <- 0L
  len <- length(x)

  nextEl <- function() {
    if (i == len) {
      stop("StopIteration")
    }

    i <<- i + 1L
    x[i]
  }

  hasNxt <- function() {
    if (i < len) {
      TRUE
    } else {
      FALSE
    }
  }

  structure(
    list(nextElem = nextEl, hasNext = hasNxt),
    class = c("indel", "abstractiter", "iter")
  )
}

indel_error <- function(mat) {
  ref <- mat[grep("^seqtype=ref", tolower(rownames(mat))), ]
  pos <- as.numeric(colnames(mat))
  ipos <- indel(pos)
  mat_ <- foreach(idl = iter.indel(ipos), .combine = "cbind") %do% {
    i <- which(pos %in% idl)
    if (all(ref[i] == ".")) {
      apply(mat[, i, drop = FALSE], 1, function(x) any(x != 1))
    } else {
      apply(mat[, i, drop = FALSE], 1, function(x) any(x != as.numeric(ref[i])))
    }
  }
  indel_start <- as.character(as.vector(ipos)[match(unique(attr(ipos, "runs")), attr(ipos, "runs"))])
  list(
    indel_num = unname(rowSums(mat_)),
    indel_pos = unname(apply(mat_, 1, function(i) paste0(indel_start[i], collapse = ":")))
  )
}

#' @keywords internal
#' @export
analyseDR2SError <- function(seqs) {
  seqnames <- names(seqs)

  ## add a unique identifier (md5 hash) to this block of sequences
  hash <- digest::digest(paste0(seqnames, collapse = ""))
  seqnames <- paste0(seqnames, hash)
  names(seqs) <- seqnames
  ##
  ##
  cns <- grep("seqtype=cns", seqnames)
  alt <- grep("cnstype=alternate", seqnames)
  seqs0 <- if (length(cns) > 0 || length(alt) > 0) {
    seqs[-union(cns, alt)]
  } else seqs
  aln <- DECIPHER::AdjustAlignment(DECIPHER::AlignSeqs(
    seqs0, verbose = FALSE, gapOpening = -12, gapExtension = -2,
    iterations = 2, refinements = 1, restrict = c(-500, 2, 10)
  ))
  mat <- DR2S:::snp_gap_matrix(DR2S:::variant_matrix(aln))
  prt <- partition(partition_reads(mat[["snp"]]))
  prt0 <- prt %>%
    dplyr::mutate(
      sample  = DR2S:::tagstring_match(read, "sample"),
      locus   = DR2S:::tagstring_match(read, "locus"),
      seqtype = DR2S:::tagstring_match(read, "seqtype"), # sequence type: 'map', 'ref', or 'cns'
      cnstype = DR2S:::tagstring_match(read, "cnstype"), # consensus type: 'LR' (longread-based consensus), 'SR' (shortread-based consensus)
      lrdtype = DR2S:::tagstring_match(read, "lrdtype"), # longread type: 'pacbio' or 'nanopore'
      nlreads = DR2S:::tagstring_match(read, "nlreads"), # number of longreads used
      reftype = ifelse(
        DR2S:::tagstring_match(read, "reftype") == "ref1" | DR2S:::tagstring_match(read, "reftype") == "ref2",
        "ref", DR2S:::tagstring_match(read, "reftype")), # reference type: 'cons', 'ref', or 'ref.ref'
      hapref  = DR2S:::tagstring_match(read, "hapref"),  # 'ref', 'alt'
      hash    = DR2S:::hash_match(read, unique = FALSE),
      cnstype = ifelse(!is.na(DR2S:::tagstring_match(read, "allele")), DR2S:::tagstring_match(read, "allele"), cnstype),
      reftype = factor(reftype, levels = c("cons", "ref", "ref.ref"), ordered = TRUE),
      hapref  = factor(hapref, levels = c("alt", "ref"), ordered = TRUE)
    ) %>%
    dplyr::select(seqtype, sample, locus, haplotype, cnstype, lrdtype, nlreads,
                  reftype, hapref, weight, hash, seqname = read) %>%
    dplyr::arrange(haplotype, desc(seqtype), cnstype, reftype, hapref)
  hapA <- dplyr::filter(prt0, haplotype == "A")
  hapB <- dplyr::filter(prt0, haplotype == "B")
  smat0 <- mat[["snp"]][hapA$seqname, ]
  smat1 <- mat[["snp"]][hapB$seqname, ]
  gmat0 <- mat[["gap"]][hapA$seqname, ]
  gmat1 <- mat[["gap"]][hapB$seqname, ]
  rs <- dplyr::bind_rows(
    dplyr::bind_cols(hapA, snp_error(mat = smat0), indel_error(gmat0)),
    dplyr::bind_cols(hapB, snp_error(mat = smat1), indel_error(gmat1))
  ) %>%
    dplyr::mutate(num_snps = NCOL(smat0)) %>%
    dplyr::arrange(haplotype, desc(seqtype), cnstype, reftype)
  aln <- Biostrings::DNAStringSetList(aln[rs$seqname])
  names(aln) <- unique(rs$hash)
  structure(list(
    errors = rs,
    aln    = aln
  ), class = c("DR2Serror", "list"))
}

tagstring_match <- function(x, tag) {
  pat <- paste0(tag, "=(.+?(?=\\|))")
  m <- regexpr(pat, x, perl = TRUE)
  start <- attr(m, "capture.start")
  end <- attr(m, "capture.start") + attr(m, "capture.length") - 1
  substr(x, start, end) %|ch|% NA_character_
}

hash_match <- function(x, unique = FALSE) {
  hash <- strsplitN(x, split = "|", n = 1, from = "end", fixed = TRUE)
  if (unique) {
    unique(hash)
  } else hash
}

make_display_name <- function(x) {
  locus     <- unique(DR2S:::tagstring_match(x, "locus"))
  sample    <- unique(DR2S:::tagstring_match(x, "sample"))
  seqtype   <- DR2S:::tagstring_match(x, "seqtype")
  allele    <- DR2S:::tagstring_match(x, "allele") %|na|% ""
  allele    <- ifelse(nzchar(allele), paste0(locus, "*", allele, " ", sample), "")
  cnstype   <- DR2S:::tagstring_match(x, "cnstype") %|na|% ""
  lrdtype   <- DR2S:::tagstring_match(x, "lrdtype") %|na|% ""
  reftype   <- DR2S:::tagstring_match(x, "reftype") %|na|% ""
  nlreads   <- DR2S:::tagstring_match(x, "nlreads") %|na|% ""
  gsub(" +", " ", paste(seqtype, allele, cnstype, lrdtype, reftype, nlreads, sep = " "))
}

#' @export
print.DR2Serror <- function(x) {
  cat("An object of class", sQuote(class(x)[1]), "with", length(x), "error tables:\n")
  x$errors %>%
    dplyr::select(
      seqtype, sample, locus, haplotype, cnstype, lrdtype, nlreads,
      reftype, hapref, err_num, indel_num, num_snps, err_pos, indel_pos
    ) %>%
    dplyr::arrange(locus, sample, haplotype, desc(seqtype), cnstype, reftype, hapref) %>%
    print()
}

#' @export
c.DR2Serror <- function(..., recursive = FALSE) {
  x <- list(...)
  structure(list(
    errors = dplyr::bind_rows(lapply(x, `[[`, "errors")) %>%
      dplyr::arrange(locus, sample, haplotype, desc(seqtype), cnstype, reftype, hapref),
    aln = do.call(c, lapply(x, `[[`, "aln"))
  ), class = c("DR2Serror", "list"))
}

#' @export
`[.DR2Serror` <- function(x, i, j, drop = FALSE) {
  hash <- dplyr::distinct(dplyr::select(x$errors, hash))[i, ]
  structure(list(
    errors = dplyr::left_join(hash, x$errors, "hash") ,
    aln    = x$aln[names(x$aln) %in% hash$hash]
  ), class = c("DR2Serror", "list"))
}

#' @export
length.DR2Serror <- function(x) {
  dplyr::n_distinct(x$errors$hash)
}

#' @export
browseDR2Serror <- function(x, i = NULL, show_aln = TRUE) {
  if (!is.null(i)) {
    x <- x[i]
  }
  by_hash <- dplyr::group_by(x$errors, hash)
  invisible(dplyr::do(by_hash, dplyr::select(
    ., seqtype, sample, locus, haplotype, cnstype, lrdtype, nlreads,
    reftype, hapref, err_num, indel_num, num_snps, err_pos, indel_pos) %>%
      dplyr::arrange(locus, sample, haplotype, desc(seqtype), cnstype, reftype, hapref) %>%
      print()
  ))

  if (show_aln) {
    invisible(foreach(seqs = x$aln) %do% {
      names(seqs) <- DR2S:::make_display_name(names(seqs))
      browse_seqs(seqs)
    })
  }

  invisible(x)
}



