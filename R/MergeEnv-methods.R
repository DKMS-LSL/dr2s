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
    haplotype = envir$haplotype,
    ## 95% confidence limit of allelic balance
    balance_upper_confint = envir$balance_upper_confint
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
                                 min_overrep_lr =5,
                                 max_overrep_sr = 2) {
  stopifnot(is(x, "variant_list"))

  warning_msg <- ""
  lr <- as.matrix(x$lr)
  varl <- filter_variant(cm = lr, threshold)

  if (!is.null(x$sr)){
    sr <- as.matrix(x$sr)
    vars <- filter_variant(cm = sr, threshold)

    ## Spurious insertion in short reads that should have been set to zero.
    if ("+" %in% names(vars)) {
      warning_msg <- warning_msg %<<% "|Insertion signal in short reads"
      # vars <- vars[!names(vars) %in% "+"]
      # sr[, "+"] <- 0
    }

    ## Mismatch in variant bases between long and short reads
    if (length(varl) > 1 && any(!names(varl) %in% names(vars))) {
      if (length(intersect(names(varl), names(vars))) == 0) {
        warning_msg <- warning_msg %<<% "|No intersect between between long and short read variants"
      } else if ("-" %in% names(varl)) {
        if (lr[varl["-"]]/sum(lr[varl]) > threshold/(2/3)) {
          warning_msg <- warning_msg %<<% "|Gap overrepresented in long reads"
        }
      } else {
        warning_msg <- warning_msg %<<% "|Mismatch between long and short read variants"
      }
    }

    vlr <- lr[, vars] + 1L ## ensure that no div/zero can occur
    vsr <- sr[, vars]

    ## total number of long reads
    N_lr <- sum(vlr)
    N_sr <- sum(vsr)
    ## Variant allele frequency
    vaf_lr <- sort(vlr/N_lr, decreasing = TRUE)
    vaf_sr <- sort(vsr/N_sr, decreasing = TRUE)
    if (length(vaf_sr) > 2) {
      warning_msg <- warning_msg %<<% "|More than two short read variants"
      vaf_lr <- vaf_lr[1:2]
      vaf_sr <- vaf_lr[1:2]
    }
    ## Error margins
    standard_error <- sqrt(prod(vaf_sr)/N_sr)
    margin_of_error <- 1.96*standard_error + (0.5/N_lr)

    overrep_sr <- vaf_sr[1]/vaf_sr[2]
    # if (overrep_sr > max_overrep_sr) {
    #   fmt <- "|Variant [%s] is %s-fold overrepresented in short reads"
    #   warning_msg <- warning_msg %<<% sprintf(fmt, names(overrep_sr), round(overrep_sr, 1))
    # }
    if (!is.na(overrep_sr) && overrep_sr < min_overrep_lr){
      fmt <- "|Variant [%s] is only %s-fold overrepresented in short reads; Keep it ambiguous"
      warning_msg <- warning_msg %<<% sprintf(fmt, names(overrep_sr), round(overrep_sr, 1))
    }
    overrep_lr <- vaf_lr[1]/vaf_lr[2]
    if (!is.na(overrep_lr) && overrep_lr < min_overrep_lr) {
      fmt <- "|Variant [%s] is only %s-fold overrepresented in long reads"
      warning_msg <- warning_msg %<<% sprintf(fmt, names(overrep_lr), round(overrep_lr, 1))
    }

    ## build variant
    bases <- names(vaf_sr)
    names(bases) <- c("ref", "alt")
    warning_msg <- sub("|", "", trimws(warning_msg), fixed = TRUE)
    x$lr_ <- lr
    x$sr_ <- sr

    x$variant <- variant(bases = bases, proportion = unname(1 - vaf_lr[2]), margin_of_error = margin_of_error, warning = warning_msg, vlist = x)
  } else {
    ## Spurious deletion signal in long reads
    if (length(varl) > 2 && "-" %in% names(varl)) {
      lrm <- lr[, varl]
      if (names(which.min(lrm/sum(lrm))) == "-") {
        ## it should be safe to silently remove the deletion if the signal
        ## for it is the weakest among the variants.
        varl <- varl[!names(varl) %in% "-"]
      }
    }

    vlr <- lr[, varl] + 1L ## ensure that no div/zero can occur
    ## total number of long reads
    N_lr <- sum(vlr)
    ## Variant allele frequency
    vaf_lr <- sort(vlr/N_lr, decreasing = TRUE)
    ## Error margins
    standard_error <- sqrt(prod(vaf_lr)/N_lr)
    margin_of_error <- 1.96*standard_error + (0.5/N_lr)
    overrep_lr <- vaf_lr[1]/vaf_lr[2]
    if (overrep_lr < min_overrep_lr) {
      fmt <- "|Variant [%s] is only %s-fold overrepresented in long reads"
      warning_msg <- warning_msg %<<% sprintf(fmt, names(overrep_lr), round(overrep_lr, 1))
    }
    ## build variant
    bases <- names(vaf_lr)
    names(bases) <- c("ref", "alt")
    warning_msg <- sub("|", "", trimws(warning_msg), fixed = TRUE)
    x$lr_ <- lr
    x$sr_ <- NULL
    x$variant <- variant(bases = bases, proportion = unname(1 - vaf_lr[2]), margin_of_error = margin_of_error, warning = warning_msg, vlist = x)
  }
  x
}

filter_variant <- function(cm, threshold) {
  cmf <- sweep(cm, 1, .rowSums(cm, NROW(cm), NCOL(cm)), `/`)
  which(apply(cmf, 2, function(col) all(col > threshold)))
}

variant <- function(bases, proportion, margin_of_error, warning, vlist) {
  refbase <- bases[["ref"]]
  altbase <- bases[["alt"]]
  base_ref_ <- if (altbase == "-") c(refbase, "+") else refbase
  base_alt_ <- if (refbase == "-") c(altbase, "+") else altbase
  if (!is.null(vlist$sr)){
    pos <- as.integer(row.names(vlist$sr))
    cm <- rbind(as.matrix(vlist$lr), vlist$lr_, as.matrix(vlist$sr), vlist$sr_)
    rownames(cm) <- c("LR1", "LR2", "SR1", "SR2")
  } else {
    pos <- as.integer(row.names(vlist$lr))
    cm <- rbind(as.matrix(vlist$lr), vlist$lr_)
    rownames(cm) <- c("LR1", "LR2")
  }
  if(!is.null(vlist$sr)){
   reads_ref = unlist(compact(ids(vlist$sr)[[as.character(pos - vlist$offset[["sr"]])]][base_ref_]))
   reads_alt = unlist(compact(ids(vlist$sr)[[as.character(pos - vlist$offset[["sr"]])]][base_alt_]))
  } else {
    reads_ref = list()
    reads_alt = list()
  }
  structure(
    bases,
    class = "variant",
    position = pos,
    proportion = proportion,
    margin_of_error = margin_of_error,
    warning = warning,
    haplotype = vlist$haplotype,
    cm = cm,
    offset = ifelse(!is.null(vlist$sr),
                    vlist$offset[["sr"]],
                    vlist$offset[["lr"]]),
    ## here we have to remove the offsets we incurred at previous
    ## ambiguous positions to match the original polymorhic positions
    ## in the pileup
    reads_ref =reads_ref,
    reads_alt = reads_alt
  )
}

#' @export
print.variant <- function(x, ...) {
  cat(sprintf("Variant at position [%s]\n", paste0(attr(x, "position") - attr(x, "offset"), collapse = "~")), sep = "")

  if (nrow(attr(x, "cm")) == 4){
    m <- attr(x, "cm")[c(1, 3), ]
    rownames(m) <- c("LR", "SR")
    m_ <- attr(x, "cm")[c(2, 4), ]
    rownames(m_) <- c("LR", "SR")
  } else {
    m <- t(as.matrix(attr(x, "cm")[1, ]))
    rownames(m) <- "LR"
    m_ <- t(as.matrix(attr(x, "cm")[2, ]))
    rownames(m_) <- "LR"
  }
  cat("\nOriginal consensus matrix:\n")
  print(m, right = TRUE, quote = FALSE)
  cat("\nDisambiguated consensus matrix:\n")
  print(m_, right = TRUE, quote = FALSE)

  cat("\nReference: ", sQuote(x[["ref"]]), "; #Sreads: ", length(attr(x, "reads_ref")), "\n", sep = "")
  cat("Alternate: ", sQuote(x[["alt"]]), "; #Sreads: ", length(attr(x, "reads_alt")), "\n", sep = "")
  cat("Shortreads reference frequency: ", round(attr(x, "proportion"), 2), "\n", sep = "")
  cat("Lower 95% limit: ", round(attr(x, "proportion") - attr(x, "margin_of_error"), 2), "\n", sep = "")

  if (nzchar(attr(x, "warning"))) {
    cat("\nWARNING: ", attr(x, "warning"), "\n", sep = "")
  }
}

`update<-` <- function(envir, pos = NULL, value) UseMethod("update<-")
`update<-.HapEnv` <- function(envir, pos = NULL, value) {
  pos <- if (is.null(pos)) envir$pos else pos
  envir$LR <- do_update_(cm = envir$LR, val = value$lr_, pos)
  if (!is.null(envir$SR)) envir$SR <- do_update_(cm = envir$SR, val = value$sr_, pos)
  envir$variants <- c(envir$variants, list(value$variant))
  envir$current_variant <- value$variant
  invisible(envir)
}

do_update_ <- function(cm, val, pos) {
  if ((n_ <- NROW(val)) > 1) {
    pos_ <- pos:(pos + n_ - 1L)
  } else {
    pos_ <- pos
  }
  cm[pos_, ] <- val
  if (anyNA(cm)) {
    ids_ <- ids(cm)
    d <- dim(cm)
    omit <- seq_along(cm)[is.na(cm)]
    omit <- unique(((omit - 1) %% d[1L]) + 1L)
    cm <- cm[-omit, , drop = FALSE]
    rownames(cm) <- as.character(seq_len(NROW(cm)))
    ids(cm) <- ids_
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

