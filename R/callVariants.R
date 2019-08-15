.callVariants <- function(x, threshold = NULL) {
  if (is.null(threshold))
    threshold <- x$getThreshold()
  variants <- lapply(x$getHapTypes(), function(hp, x, threshold) {
    envir <- list()
    if (x$hasShortreads()) {
      lr <- consmat(x$remap$LR[[hp]]$pileup, prob = FALSE)
      sr <- consmat(x$remap$SR[[hp]]$pileup, prob = FALSE)
      rs <- .equaliseConsmat(lrm = lr, srm = sr)
      envir$LR <- rs$lrm
      envir$SR <- rs$srm
      referencePath <- x$absPath(x$remap$SR[[hp]]$refpath)
      reference <- Biostrings::readDNAStringSet(referencePath)
      envir$ref <- unname(strsplit1(as.character(reference), ""))
    } else {
      envir$LR <- consmat(x$remap$LR[[hp]]$pileup, prob = FALSE)
      envir$SR <- NULL
      referencePath <- x$absPath(x$remap$LR[[hp]]$refpath)
      reference <- Biostrings::readDNAStringSet(referencePath)
      envir$ref <- unname(strsplit1(as.character(reference), ""))
    }
    apos <- foreach(rt = c("LR", "SR"), .combine = c) %do% {
      ## positions not matching the consensus
      dism <- .noRefMatch(envir[[rt]], envir[["ref"]])
      ## ambiguous positions
      amb <- .ambiguousPositions(envir[[rt]], threshold, FALSE)
      sort(c(amb, dism))
    }
    envir$pos <- unique(sort(apos))
    envir$threshold <- threshold
    compact(.walkVariants(envir))
  }, x = x, threshold = threshold)
  names(variants) <- x$getHapTypes()
  dplyr::arrange(.getVariants(variants), as.numeric(.data$pos), .data$haplotype)
}

.walkVariants <- function(envir) {
  lapply(envir$pos, function(pos) {
    pos <- ifelse(!is.null(envir$SR),
                  pos + offsetBases(envir$SR),
                  pos + offsetBases(envir$LR))
    envir$currentPos <- pos
    # message(pos)
    disambiguateVariant(x = yield(envir, pos), threshold = envir$threshold)
  })
}

# Helpers -----------------------------------------------------------------

.expandLongreadConsmat <- function(lrm, srm) {
  m <- as.matrix(lrm)
  if (NCOL(lrm) == 5) {
    m <- cbind(m, `+` = 0)
  }
  if (length(ins(srm)) > 0) {
    insert <- matrix(c(0, 0, 0, 0, stats::median(rowSums(m)), 0), ncol = 6)
    myIns <- sort(ins(srm))
    myIns <- myIns[myIns < nrow(lrm)]
    INSit <- itertools::ihasNext(iterators::iter(myIns))
    while (itertools::hasNext(INSit)) {
      i <- iterators::nextElem(INSit)
      m <- rbind(m[seq_len(i - 1), ], insert, m[i:NROW(m), ])
    }
    if (!NROW(m) == NROW(srm)) {
      warning("SR and LR of different length! Check problem file")
      if (NROW(m) < NROW(srm)) {
        flog.info(" fill longreads with gaps from %s to %s",
                NROW(m), NROW(srm), name = "info")
        add <- ((NROW(m) + 1):NROW(srm))
        m <- rbind(m, srm[add,])
        m[add,] <- rep.int(0, 6*length(add))
      } else if (NROW(m) > NROW(srm)) {
        flog.info(" fill shortreads with gaps from %s to %s",
                NROW(srm), NROW(m), name = "info")
        add <- ((NROW(srm) + 1):NROW(m))
        srm <- rbind(srm, m[add,])
        srm[add,] <- rep.int(0, 6*length(add))
      }
    }
    stopifnot(NROW(m) == NROW(srm))
    dimnames(m) <- dimnames(srm)
  }
  lrm <- consmat(m, freq = FALSE)
  lrm
}


.equaliseConsmat <- function(lrm, srm) {
  if (NROW(srm) != NROW(lrm) + length(ins(srm))) {
    flog.warn("SR and LR are of different length", name = "info")
  }
  sm <- as.matrix(srm)
  if (NCOL(srm) == 5) {
    sm <- cbind(sm, `+` = 0)
  }
  # if (length(ins(lrm)) > 0) {
  #   insert <- matrix(c(0, 0, 0, 0, stats::median(rowSums(sm)), 0), ncol = 6)
  #   myIns  <- sort(ins(lrm))
  #   myIns  <- myIns[myIns < NROW(srm) + length(ins(lrm))]
  #   INSit  <- itertools::ihasNext(iterators::iter(myIns))
  #   if (length(ins(srm)) > 0) {
  #     while (itertools::hasNext(INSit)) {
  #       i  <- iterators::nextElem(INSit)
  #       if (i %in% ins(srm))
  #         next
  #       ins(srm)[ins(srm) > i] <- ins(srm)[ins(srm) > i] + 1
  #       sm <- rbind(sm[seq_len(i - 1), ], insert, sm[i:NROW(sm), ])
  #     }
  #   } else {
  #     while (itertools::hasNext(INSit)) {
  #       i  <- iterators::nextElem(INSit)
  #       sm <- rbind(sm[seq_len(i - 1), ], insert, sm[i:NROW(sm), ])
  #     }
  #   }
  # }
  dimnames(sm) <- list(pos = as.character(1:NROW(sm)),
                       nucleotide = c("G", "A", "T", "C", "-", "+"))

  lm <- as.matrix(lrm)
  if (NCOL(lrm) == 5) {
    lm <- cbind(lm, `+` = 0)
  }
  if (length(ins(srm)) > 0) {
    insert <- matrix(c(0, 0, 0, 0, stats::median(rowSums(lm)), 0), ncol = 6)
    myIns  <- sort(ins(srm))
    
    ## Don't add the length(of ins(srm) to NROW(lrm) for getting valid 
    ## insertions. We don't know the length of each ins and so its useless
    # myIns  <- myIns[myIns < NROW(lrm) + length(ins(srm))]
    myIns  <- myIns[myIns < NROW(lrm)]
    INSit  <- itertools::ihasNext(iterators::iter(myIns))
    if (length(ins(lrm)) > 0) {
      while (itertools::hasNext(INSit)) {
        i  <- iterators::nextElem(INSit)
        if (i %in% ins(lrm))
          next
        ins(lrm)[ins(lrm) > i] <- ins(lrm)[ins(lrm) > i] + 1
        lm <- rbind(lm[seq_len(i - 1), ], insert, lm[i:NROW(lm), ])
      }
    } else {
      while (itertools::hasNext(INSit)) {
        i  <- iterators::nextElem(INSit)
        lm <- rbind(lm[seq_len(i - 1), ], insert, lm[i:NROW(lm), ])
      }
    }
  }
  dimnames(lm) <- list(pos = as.character(1:NROW(lm)),
                       nucleotide = c("G", "A", "T", "C", "-", "+"))

  if (NROW(sm) != NROW(lm)) {
    flog.warn("SR and LR of different length! Check problem file", name = "info")
    if (NROW(lm) < NROW(sm)) {
      flog.info(" fill longreads with gaps from %s to %s",
                NROW(sm), NROW(sm), name = "info")
      add <- ((NROW(lm) + 1):NROW(sm))
      lm <- rbind(lm, sm[add, ])
      lm[add, ] <- rep.int(0, 6*length(add))
    } else if (NROW(lm) > NROW(sm)) {
      flog.info(" fill shortreads with gaps from %s to %s",
                NROW(sm) + 1, NROW(lm), name = "info")
      add <- ((NROW(sm) + 1):NROW(lm))
      sm  <- rbind(sm, lm[add, ])
      sm[add, ] <- rep.int(0, 6*length(add))
    }
  }

  stopifnot(NROW(sm) == NROW(lm))

  sm.out <- consmat(sm, freq = FALSE)
  ins(sm.out) <- ins(srm)

  lm.out <- consmat(lm, freq = FALSE)
  ins(lm.out) <- ins(lrm)

  list(srm = sm.out, lrm = lm.out)
}


yield <- function(envir, pos) {
  lr <- envir$LR
  sr <- envir$SR
  ref <- envir$ref
  structure(list(
    ## class: variant
    variant = NULL,
    ## Original consensus matrix at the position of the variant
    lr = lr[pos, ],
    sr = sr[pos, ],
    ref = ref[pos],
    ## Disambiguated consensus matrix at the position of the variant
    lr_ = NULL,
    sr_ = NULL,
    ##
    offsetBases = c(
      lr = offsetBases(lr),
      sr = ifelse(is.null(sr), 0, offsetBases(sr)),
      ref = ifelse(is.null(sr), offsetBases(lr), offsetBases(sr))
    ),
    ##
    haplotype = envir$haplotype
    ##
  ), class = c("variantList", "list"))
}
print.variantList <- function(x, threshold = 0.2, ...) {
  if (!is.null(x$variant)) {
    print(x$variant)
  } else {
    print("No variant called at position.")
  }
}

#' Disambiguate a variant
#'
#' @param x A \code{variantList}.
#' @param threshold Minimum frequency for a variant.
#'
#' @return A \code{variantList}.
#' @keywords internal
#' @examples
#' ###
#'
disambiguateVariant <- function(x, threshold) {
  assert_that(is(x, "variantList"))
  warningMsg <- ""

  ## Start with long reads. They have to be there
  lr <- as.matrix(x$lr)
  varl <- .filterVariant(cm = lr, threshold, ignoreInsertions = FALSE)
  ref <- x$ref

  ## Look if its ambiguous
  if (length(varl) > 1) {
    ## Dont look for deletions in longreads
    # ## Look for deletions
    if (all(names(varl) %in% VALID_DNA("none"))) {
      warningMsg <- warningMsg %<<% "|Ambiguous position in long reads"
    } else if ("+" %in% names(varl)) {
      warningMsg <- warningMsg %<<% "|Insertion signal in long reads"
    }
    ## Check for > 2 alleles
    if (length(varl) > 2) {
      warningMsg <- warningMsg %<<% "|More than two long read variants"
      lrBases <- names(varl)[1:2]
    }
    if ("-" %in% names(varl) & !is.null(ref) & names(varl)[names(varl) != "-"] != ref) {
      warningMsg <- warningMsg %<<% "|long read variant not matching the reference"
      lrBases <- c(names(varl), names(varl))
    }
    ## Set the ambiguous bases
    lrBases <- names(varl)
  } else if (sum(lr) > 0) {
    if (!is.null(ref) & names(varl) != ref) {
      if (names(varl) == "-") {
        lrBases <- c(names(varl), names(varl))
      } else {
        warningMsg <- warningMsg %<<% "|long read variant not matching the reference"
        lrBases <- c(names(varl), names(varl))
      }
    } else {
      lrBases <- c(names(varl), names(varl))
      if (is.null(ref))
        warningMsg <- warningMsg %<<% "|no reference "
    }
  } else {
    lrBases <- rep(NA, 2)
    warningMsg <- warningMsg %<<% "|no longreads"
  }

  if (!is.null(x$sr)) {
    sr <- as.matrix(x$sr)
    vars <- .filterVariant(cm = sr, threshold, ignoreInsertions = FALSE)
    if (sum(sr) == 0) {
      warningMsg <- warningMsg %<<% "|no shortreads"
      srBases <- c("-", "-")
    } else {
      ## State which reads are ambiguous
      if (length(vars) > 1) {
        srBases <- names(vars)
        ## Check for deletions
        if ("-" %in% names(vars)) {
          if ((sr[vars["-"]]/sum(sr[vars])) %|na|% 0  > max(threshold, 0.3)) {
            warningMsg <- warningMsg %<<% "|Gap in short reads"
          }
        } else if ("+" %in% names(vars)) {
          warningMsg <- warningMsg %<<% "|Insertion signal in short reads"
        } else {
          warningMsg <- warningMsg %<<% "|Ambiguous position in short reads"
          ## Mismatch in variant bases between long and short reads
          ## not the same order or base
          if (length(varl) == length(vars)) {
            if (any(!varl == vars)) {
              ## No base matches
              if (length(intersect(names(varl), names(vars))) == 0) {
                warningMsg <- warningMsg %<<%
                  "|No intersect between long and short read variants"
              } else if (any(!sort(varl) == sort(vars))) { # different order
                  warningMsg <- warningMsg %<<%
                    "|Major/minor variant different in long and short reads"
              } else if (length(varl) > length(vars)) {
                  warningMsg <- warningMsg %<<%
                    "|Variant in long but not in short reads"
              } else if (length(vars) > length(varl)) {
                  warningMsg <- warningMsg %<<%
                    "|Variant in short but not in long reads"
              } else {
                warningMsg <- warningMsg %<<%
                  "|Mismatch between long and short read variants"
              }
            }
          }
  
          ## Check for more than two alleles
          if (length(vars) > 2) {
            warningMsg <- warningMsg %<<% "|More than two short read variants"
            srBases <- names(vars)[1:2]
          }
        }
      } else if (names(vars) != ref) {
        ## Check if it is a gap
        if (names(vars) == "-") {
          warningMsg <- warningMsg %<<% "|Only gap in short read variant, not matching the reference"
        } else {
          warningMsg <- warningMsg %<<% "|short read variant not matching the reference"
        }
        srBases <- c(names(vars), names(vars))
      } else {
        ## Use only non-gap variants from only longreads.
        warningMsg <- ifelse("-" %in% names(varl),
                             "",
                             ifelse(length(varl) > 1, "|Variant only in long reads", "")) %<<% warningMsg
        srBases <- c(names(vars), names(vars))
      }
    }
  } else {
    srBases <- rep(NA, 2)
  }

  # names(bases) <- c("ref", "alt")
  warningMsg <- sub("|", "", trimws(warningMsg), fixed = TRUE)
  #  IGNORE IF ONLY GAPS????
  if (nzchar(warningMsg))
    return(.variant(lrBases = lrBases, srBases = srBases, 
                          lrSupport = attr(varl, "support"), 
                          srSupport = attr(vars, "support"),
                          warning = warningMsg, vlist = x))
  NULL
}

.filterVariant <- function(cm, threshold, ignoreInsertions = TRUE) {
  ## ignore insertions
  if (ignoreInsertions) {
    bases <- VALID_DNA("del")
  } else {
    bases <- VALID_DNA("indel")
  }
  cm <- cm[, bases, drop = FALSE]
  cmf <- sweep(cm, 1, .rowSums(cm, NROW(cm), NCOL(cm)), `/`)
  vars <- which(apply(cmf, 2, function(col) all(col > threshold)))
  vars <- vars[order(cmf[vars], decreasing = TRUE)]
  attr(vars, "support") <-  max(cmf[vars])
  vars
}

.variant <- function(lrBases, srBases, lrSupport, srSupport, warning, vlist) {
  refbase <- ifelse(any(is.na(srBases)), lrBases[[1]], srBases[[1]])
  altbase <- ifelse(any(is.na(srBases)), lrBases[[2]], srBases[[2]])
  if (!is.null(vlist$sr)) {
    pos <- as.integer(row.names(vlist$sr))
    cm <- rbind(as.matrix(vlist$lr), as.matrix(vlist$sr))
    rownames(cm) <- c("LR", "SR")
  } else {
    pos <- as.integer(row.names(vlist$lr))
    cm <- as.matrix(vlist$lr)
    rownames(cm) <- "LR"
  }

  structure(
    c(refbase, altbase),
    class       = "variant",
    lrBases     = lrBases,
    srBases     = srBases,
    lrSupport = lrSupport,
    srSupport = srSupport,
    position    = pos,
    warning     = warning,
    haplotype   = vlist$haplotype,
    cm          = cm,
    offsetBases = ifelse(!is.null(vlist$sr),
                    vlist$offsetBases[["sr"]],
                    vlist$offsetBases[["lr"]]),
    readsRefSr  = vlist$sr[,refbase],
    readsRefLr  = vlist$lr[,refbase],
    readsAltSr  = vlist$sr[,altbase],
    readsAltLr  = vlist$lr[,altbase]
    ## here we have to remove the offsets we incurred at previous
    ## ambiguous positions to match the original polymorhic positions
    ## in the pileup
  )
}

#' @export
print.variant <- function(x, ...) {
  cat(sprintf("Variant at position [%s]\n",
              paste0(attr(x, "position") - attr(x, "offsetBases"),
                     collapse = "~")), sep = "")

  m <- attr(x, "cm")
  cat("\nConsensus matrix:\n")
  print(m, right = TRUE, quote = FALSE)

  cat("\nReference: ", sQuote(x[[1]]),
      "; #Sreads: ", attr(x, "readsRefSr"),
      "; #Lreads: ", attr(x, "readsRefLr"),
      "\n", sep = "")
  cat("\nAlternate: ", sQuote(x[[2]]),
      "; #Sreads: ", attr(x, "readsAltSr"),
      "; #Lreads: ", attr(x, "readsAltLr"),
      "\n", sep = "")

  if (nzchar(attr(x, "warning"))) {
    cat("\nWARNING: ", attr(x, "warning"), "\n", sep = "")
  }
}