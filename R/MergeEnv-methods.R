yield <- function(envir, ...) UseMethod("yield")
yield.HapEnv <- function(envir, pos = NULL) {
  lr <- envir$LR
  sr <- envir$SR

  pos <- if (is.null(pos)) envir$pos else pos
  structure(list(
    ## class: variant
    variant = NULL,
    ## Original consensus matrix at the position of the variant
    lr = lr[pos, ],
    sr = sr[pos, ],
    ## Disambiguated consensus matrix at the position of the variant
    lr_ = NULL,
    sr_ = NULL,
    ##
    offsetBases = c(
      lr = offsetBases(lr),
      sr = ifelse(is.null(sr), 0, offsetBases(sr))
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
    m <- rbind(x$lr, x$sr)
    grp <- "Pos [" %<<% colon(unique(rownames(m))) %<<% "]"
    rownames(m) <- if(NROW(m) == 2) (c("LR", "SR")) else  ("LR")

    cat(sprintf("%+10s", grp), sep = "\n")
    print(m, right = TRUE, quote = FALSE)

    if (!is.null(x$lr_)) {
      m_ <- rbind(x$lr_, x$sr_)
      rownames(m_) <- c("LR", "SR")

      cat("\nDisambiguated consensus matrix:\n")
      print(m_, right = TRUE, quote = FALSE)
    }
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
  varl <- .filterVariant(cm = lr, threshold)

  ## Look if its ambiguous
  if (length(varl) > 1) {
    ## Look for deletions
    if ("-" %in% names(varl)) {
      if (lr[varl["-"]]/sum(lr[varl]) > max(threshold/(2/3), 0.3)) {
        warningMsg <- warningMsg %<<% "|Gap overrepresented in long reads"
      }
    } else {
      warningMsg <- warningMsg %<<% "|Ambiguous position in long reads"
    }
    ## Check for > 2 alleles
    if (length(varl) > 2) {
      warningMsg <- warningMsg %<<% "|More than two long read variants"
      srBasel <- names(varl)[1:2]
    }
    ## Set the ambiguous bases
    lrBases <- names(varl)
  } else {
    warningMsg <- warningMsg %<<% "|Variant only in short reads"
    lrBases <- c(names(varl), NA)
  }

  if (!is.null(x$sr)) {
    sr <- as.matrix(x$sr)
    vars <- .filterVariant(cm = sr, threshold)

    ## State which reads are ambiguous
    if (length(vars) > 1 ) {
      srBases <- names(vars)
      ## Check for deletions
      if ("-" %in% names(vars)) {
        if ((sr[vars["-"]]/sum(sr[vars])) %|na|% 0  > max(threshold/2*3, 0.3)) {
          warningMsg <- warningMsg %<<% "|Gap in short reads"
        }
      } else {
        warningMsg <- warningMsg %<<% "|Ambiguous position in short reads"

        ## Spurious insertion in short reads that should have been set to zero.
        if ("+" %in% names(vars)) {
          warningMsg <- warningMsg %<<% "|Insertion signal in short reads"
        }

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
    } else {
      ## Use only non-gap variants from only longreads.
      warningMsg <- ifelse("-" %in% names(varl),
                           "",
                           "|Variant only in long reads") %<<% warningMsg
      srBases <- c(names(vars), NA)
    }
  } else {
    srBases <- c(NA, NA)
  }

  # names(bases) <- c("ref", "alt")
  warningMsg <- sub("|", "", trimws(warningMsg), fixed = TRUE)
  #  IGNORE IF ONLY GAPS????
  if (nzchar(warningMsg))
    x$variant <- .variant(lrBases = lrBases, srBases = srBases,
                          warning = warningMsg, vlist = x)
  x
}

.filterVariant <- function(cm, threshold) {
  ## ignore insertions
  cm <- cm[, VALID_DNA("del"), drop = FALSE]
  cmf <- sweep(cm, 1, .rowSums(cm, NROW(cm), NCOL(cm)), `/`)
  which(apply(cmf, 2, function(col) all(col > threshold)))
}

.variant <- function(lrBases, srBases, warning, vlist) {
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
    class = "variant",
    lrBases = lrBases,
    srBases = srBases,
    position = pos,
    warning = warning,
    haplotype = vlist$haplotype,
    cm = cm,
    offsetBases = ifelse(!is.null(vlist$sr),
                    vlist$offsetBases[["sr"]],
                    vlist$offsetBases[["lr"]]),
    readsRefSr = vlist$sr[,refbase],
    readsRefLr = vlist$lr[,refbase],
    readsAltSr = vlist$sr[,altbase],
    readsAltLr = vlist$lr[,altbase]
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

`.update<-` <- function(envir, pos = NULL, value) UseMethod(".update<-")
`.update<-.HapEnv` <- function(envir, pos = NULL, value) {
  pos <- if (is.null(pos)) envir$pos else pos
  envir$LR <- .doUpdate_(cm = envir$LR, pos)
  if (!is.null(envir$SR)) envir$SR <- .doUpdate_(cm = envir$SR, pos)
  envir$variants <- c(envir$variants, list(value$variant))
  envir$currentVariant <- value$variant
  invisible(envir)
}

.doUpdate_ <- function(cm, pos) {
  pos_ <- pos
  if (anyNA(cm)) {
    d <- dim(cm)
    omit <- seq_along(cm)[is.na(cm)]
    omit <- unique(((omit - 1) %% d[1L]) + 1L)
    cm <- cm[-omit, , drop = FALSE]
    rownames(cm) <- as.character(seq_len(NROW(cm)))
    i <- which((ins_ <- ins(cm)) == pos - offsetBases(cm))
    if (length(i) == 1) {
      run_ <- attr(ins_, "run")
      ins_ <- ins_[!run_ == run_[i]]
      run_ <- run_[!run_ == run_[i]]
      attr(ins_, "run") <- run_
    }
    ins(cm) <- ins_
    offsetBases(cm) <- offsetBases(cm)
  }
  ## update rowsums
  n(cm) <- .rowSums(cm, NROW(cm), NCOL(cm))
  cm
}

## TODO: rm lr_ and sr_ ????
