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
    offset = c(
      lr = offset(lr),
      sr = ifelse(is.null(sr), 0, offset(sr))
    ),
    ##
    haplotype = envir$haplotype
    ##
  ), class = c("variant_list", "list"))
}
print.variant_list <- function(x, threshold = 0.2, ...) {
  if (!is.null(x$variant)) {
    print(x$variant)
  } else {
    m <- rbind(x$lr, x$sr)
    grp <- paste0("Pos [", colon(unique(rownames(m))), "]")
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
#' @param x A \code{variant_list}.
#' @param threshold Minimum frequency for a variant.
#' @param min_overrep_lr Issue a warning if the majority variant in the
#' long reads is only overrepresented by less than \code{min_overrep_lr}.
#' @param max_overrep_sr Issue a warning if the majority variant in the
#' short reads is overrepresented by more than \code{max_overrep_sr}.
#'
#' @return A \code{variant_list}.
#' @keywords internal
#' @examples
#' ###
#'
                                 # min_overrep_lr = 1.5
                                 # max_overrep_sr = 2
#x <- a
disambiguate_variant <- function(x,
                                 threshold,
                                 min_overrep_lr = 2,
                                 max_overrep_sr = 2) {
  stopifnot(is(x, "variant_list"))

  warning_msg <- ""
  
  ## Start with long reads. They have to be there
  lr <- as.matrix(x$lr)
  varl <- .filterVariant(cm = lr, threshold)
  
  ## Look if its ambiguous
  if (length(varl) > 1 ){
    warning_msg <- warning_msg %<<% "|Ambiguous position in long reads"
  }

  ## Look for deletions
  if ("-" %in% names(varl)) { 
    if (lr[varl["-"]]/sum(lr[varl]) > threshold/(2/3)) {
      warning_msg <- warning_msg %<<% "|Gap overrepresented in long reads"
    }
  }
  
  ## Check for > 2 alleles
  if (length(varl) > 2) {
    warning_msg <- warning_msg %<<% "|More than two long read variants"
  }

  ## Set the ambiguous bases
  bases <- names(varl) 
  
  if (!is.null(x$sr)){
    sr <- as.matrix(x$sr)
    vars <- .filterVariant(cm = sr, threshold)
    
    ## State which reads are ambiguous
    if (length(vars) > 1 ){
      warning_msg <- warning_msg %<<% "|Ambiguous position in short reads"
    }

    ## Spurious insertion in short reads that should have been set to zero.
    if ("+" %in% names(vars)) {
      warning_msg <- warning_msg %<<% "|Insertion signal in short reads"
    }
    
    ## Mismatch in variant bases between long and short reads
    ## not the same order or base
    if (any(!varl == vars)) {
      ## No base matches
      if (length(intersect(names(varl), names(vars))) == 0) {
        warning_msg <- warning_msg %<<% 
          "|No intersect between long and short read variants"
      } else if (any(!sort(varl) == sort(vars))) { # different order
          warning_msg <- warning_msg %<<% 
            "|Major/minor variant different in long and short reads"
      } else if (length(varl) > length(vars)) { 
          warning_msg <- warning_msg %<<% 
            "|Variant in long but not in short reads"
      } else if (length(vars) > length(varl)) { 
          warning_msg <- warning_msg %<<% 
            "|Variant in short but not in long reads"
      } else {
        warning_msg <- warning_msg %<<% 
          "|Mismatch between long and short read variants"
      }
    }
    
    ## Check for deletions
    if ("-" %in% names(vars)) { 
      if ((sr[vars["-"]]/sum(sr[vars])) %|na|% 0  > threshold/(2/3)) {
        warning_msg <- warning_msg %<<% "|Gap overrepresented in short reads"
      }
    }
    
    ## Check for more than two alleles
    if (length(vars) > 2) {
      warning_msg <- warning_msg %<<% "|More than two short read variants"
    }
    
    ## Set the ambiguous bases to short read variants if not fewer
    if (!length(vars) < length(varl))  
      bases <- names(vars)
  }
  
  # names(bases) <- c("ref", "alt")
  warning_msg <- sub("|", "", trimws(warning_msg), fixed = TRUE)
  
  x$variant <- variant(bases = bases, warning = warning_msg, vlist = x)
  # x <- x$variant
  x
  # } 
  # x$variant <- variant(bases = bases, warning = warning_msg, vlist = x)
    # x$variant <- variant(bases = bases, proportion = unname(1 - vaf_lr[2]), margin_of_error = margin_of_error, warning = warning_msg, vlist = x)
    # ## Spurious deletion signal in long reads
    # if (length(varl) > 2 && "-" %in% names(varl)) {
    #   lrm <- lr[, varl]
    #   if (names(which.min(lrm/sum(lrm))) == "-") {
    #     ## it should be safe to silently remove the deletion if the signal
    #     ## for it is the weakest among the variants.
    #     varl <- varl[!names(varl) %in% "-"]
    #   }
    # }

    # vlr <- lr[, varl] + 1L ## ensure that no div/zero can occur
    # ## total number of long reads
    # N_lr <- sum(vlr)
    # ## Variant allele frequency
    # vaf_lr <- sort(vlr/N_lr, decreasing = TRUE)
    # ## Error margins
    # standard_error <- sqrt(prod(vaf_lr)/N_lr)
    # margin_of_error <- 1.96*standard_error + (0.5/N_lr)
    # overrep_lr <- vaf_lr[1]/vaf_lr[2]
    # if (overrep_lr < min_overrep_lr) {
    #   fmt <- "|Variant [%s] is only %s-fold overrepresented in long reads"
    #   warning_msg <- warning_msg %<<% sprintf(fmt, names(overrep_lr), round(overrep_lr, 1))
    # }
    ## build variant
    # bases <- names(vaf_lr)
    # names(bases) <- c("ref", "alt")
    # warning_msg <- sub("|", "", trimws(warning_msg), fixed = TRUE)
    # x$lr_ <- lr
    # x$sr_ <- NULL
    # x$variant <- variant(bases = bases, proportion = unname(1 - vaf_lr[2]), margin_of_error = margin_of_error, warning = warning_msg, vlist = x)
  # }
  x
}

.filterVariant <- function(cm, threshold) {
  cmf <- sweep(cm, 1, .rowSums(cm, NROW(cm), NCOL(cm)), `/`)
  which(apply(cmf, 2, function(col) all(col > threshold)))
}

variant <- function(bases, warning, vlist) {
  refbase <- bases[[1]]
  altbase <- bases[[2]]
  base_ref_ <- if (altbase == "-") c(refbase, "+") else refbase
  base_alt_ <- if (refbase == "-") c(altbase, "+") else altbase
  if (!is.null(vlist$sr)){
    pos <- as.integer(row.names(vlist$sr))
    cm <- rbind(as.matrix(vlist$lr), as.matrix(vlist$sr))
    rownames(cm) <- c("LR", "SR")
  } else {
    pos <- as.integer(row.names(vlist$lr))
    cm <- as.matrix(vlist$lr)
    rownames(cm) <- "LR"
  }
  
  structure(
    bases,
    class = "variant",
    position = pos,
    warning = warning,
    haplotype = vlist$haplotype,
    cm = cm,
    offset = ifelse(!is.null(vlist$sr),
                    vlist$offset[["sr"]],
                    vlist$offset[["lr"]]),
    reads_ref_sr = vlist$sr[,refbase],
    reads_ref_lr = vlist$lr[,altbase],
    reads_alt_sr = vlist$sr[,refbase],
    reads_alt_lr = vlist$lr[,altbase]
    ## here we have to remove the offsets we incurred at previous
    ## ambiguous positions to match the original polymorhic positions
    ## in the pileup
  )
}

#' @export
print.variant <- function(x, ...) {
  cat(sprintf("Variant at position [%s]\n", 
              paste0(attr(x, "position") - attr(x, "offset"), 
                     collapse = "~")), sep = "")

  m <- attr(x, "cm")
  cat("\nConsensus matrix:\n")
  print(m, right = TRUE, quote = FALSE)

  cat("\nReference: ", sQuote(x[[1]]), 
      "; #Sreads: ", attr(x, "reads_ref_sr"), 
      "; #Lreads: ", attr(x, "reads_ref_lr"), 
      "\n", sep = "")
  cat("\nAlternate: ", sQuote(x[[2]]), 
      "; #Sreads: ", attr(x, "reads_alt_sr"), 
      "; #Lreads: ", attr(x, "reads_alt_lr"), 
      "\n", sep = "")

  if (nzchar(attr(x, "warning"))) {
    cat("\nWARNING: ", attr(x, "warning"), "\n", sep = "")
  }
}

`update<-` <- function(envir, pos = NULL, value) UseMethod("update<-")
`update<-.HapEnv` <- function(envir, pos = NULL, value) {
  pos <- if (is.null(pos)) envir$pos else pos
  envir$LR <- do_update_(cm = envir$LR, pos)
  if (!is.null(envir$SR)) envir$SR <- do_update_(cm = envir$SR, pos)
  envir$variants <- c(envir$variants, list(value$variant))
  envir$current_variant <- value$variant
  invisible(envir)
}

do_update_ <- function(cm, pos) {
  pos_ <- pos
  if (anyNA(cm)) {
    d <- dim(cm)
    omit <- seq_along(cm)[is.na(cm)]
    omit <- unique(((omit - 1) %% d[1L]) + 1L)
    cm <- cm[-omit, , drop = FALSE]
    rownames(cm) <- as.character(seq_len(NROW(cm)))
    i <- which((ins_ <- ins(cm)) == pos - offset(cm))
    if (length(i) == 1) {
      run_ <- attr(ins_, "run")
      ins_ <- ins_[!run_ == run_[i]]
      run_ <- run_[!run_ == run_[i]]
      attr(ins_, "run") <- run_
    }
    ins(cm) <- ins_
    offset(cm) <- offset(cm) - n_
  }
  ## update rowsums
  n(cm) <- .rowSums(cm, NROW(cm), NCOL(cm))
  cm
}

