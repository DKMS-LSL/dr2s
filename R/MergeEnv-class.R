#' @include MergeEnv-methods.R
NULL


# Class: MergeEnv ---------------------------------------------------------


#' Constructor for \code{\link[=MergeEnv_]{MergeEnv}} objects.
#'
#' @param x A \code{\link[=DR2S_]{DR2S}} object.
#' @param threshold When do we call a variant a variant.
#'
#' @return A \code{\link[=MergeEnv_]{MergeEnv}} object.
#' @export
#'
#' @examples
#' ###
MergeEnv <- function(x, threshold) {
  assertthat::assert_that(is(x, "DR2S"))
  if (missing(threshold)) {
    threshold <- x$getThreshold()
  }
  MergeEnv_$new(x, threshold)
}

#' Class \code{"MergeEnv"}
#'
#' @docType class
#' @usage MergeEnv(x)
#' @field threshold \code{[numeric]}; when do we call a variant a variant.
#' @field a \code{[HapEnv, environment]}; container for all things haplotype A.
#' @field b \code{[HapEnv, environment]}; container for all things haplotype B.
#' @field x \code{[\link[=DR2S_]{DR2S}]}; the original \code{DR2S} object.
#' @keywords data internal
#' @return Object of \code{\link{R6Class}} representing a MergeEnv.
#' @section Public Methods:
#' \describe{
#' \item{\code{x$init(hapEnv)}}{Intialise a \code{HapEnv}}
#' \item{\code{x$walk(hapEnv, until = FALSE, verbose = FALSE)}}{Perform a SNP
#' walk along the consensus matrices until there are no SNPs to be walked
#' \code{until = FALSE}, until a certain position (\code{until = 1234}), or
#' one step at a time (\code{until = NULL}).}
#' }
MergeEnv_ <- R6::R6Class(
  classname = "MergeEnv",
  public = list(
    threshold = NA_real_,
    hptypes   = NA,
    x = NA,
    initialize = function(x, threshold = x$getThreshold()) {
      assertthat::assert_that(is(x, "DR2S"))
      self$hptypes <- foreach(hptype = x$getHapTypes(), .final = function(h) setNames(h, x$getHapTypes())) %do% {
       structure(as.environment(list(haplotype = hptype)), class = "HapEnv")
      }
      self$threshold = threshold
      self$x = x
    },

    isInitialised = function(hapEnv) {
      hapEnv <- match.arg(hapEnv, self$x$getHapTypes())
      !is.null(self$hptypes[[hapEnv]]$init)
    },
    getHapEnv = function(hapEnv) {
      hapEnv <- match.arg(hapEnv, self$x$getHapTypes())
      self$hptypes[[hapEnv]]
    },
    currentVariant = function(envir) {
      if (!is.null(envir$current_variant)) {
        print(envir$current_variant)
      } else NULL
    },
    print = function(left = 3, right = left, ...) {
      for (hptype in self$hptypes){
        self$showMatrix(hptype, left = left, right = right)
        cat("\n")
      }
    }
  )
)

MergeEnv_$set("public", "init", function(hapEnv) {
  hapEnv <- match.arg(hapEnv, self$x$getHapTypes())
  envir <- self$hptypes[[hapEnv]]
  envir$SR <- self$x$mapFinal$pileup[[paste0("SR", hapEnv)]]$consmat ## consmat short reads
  lr <- self$x$mapFinal$pileup[[paste0("LR", hapEnv)]]$consmat ## consmat long reads

  envir$LR = expand_longread_consmat(lrm = lr, srm = envir$SR)
  envir$variants = list()
  envir$POSit = hlatools::ihasNext(iter(apos <- ambiguous_positions(envir$SR, self$threshold)))
  envir$pos = 1L
  envir$current_variant = NULL
  envir$init = TRUE
  envir$balance = apply(as.matrix(envir$SR[apos, ]), 1, function(x) {
    tmp <- sort(x, decreasing = TRUE)[1:2]
    tmp[1] / tmp[2]
  })
  envir$balance_upper_confint <- sum(mean(envir$balance, na.rm = TRUE), 1.96*sd(envir$balance, na.rm = TRUE) , na.rm = TRUE)
})

## self$walk_one() ####
MergeEnv_$set("public", "walk_one", function(hapEnv, verbose = FALSE) {
  hapEnv <- match.arg(hapEnv, self$x$getHapTypes())
  if (!self$isInitialised(hapEnv)) {
    self$init(hapEnv)
  }

  envir <- self$getHapEnv(hapEnv)
  private$step_through(envir)
  if (verbose) {
    self$currentVariant(envir)
    cat("\nConsensus sequence:\n")
    self$showConsensus(envir)
  }
  invisible(self)
})

## self$walk() ####
# self <- menv
# hapEnv <- "A"
MergeEnv_$set("public", "walk", function(hapEnv, verbose = FALSE) {
  hapEnv <- match.arg(hapEnv, self$x$getHapTypes())
  if (!self$isInitialised(hapEnv)) {
    self$init(hapEnv)
  }
  envir <- self$getHapEnv(hapEnv)
  while (private$step_through(envir)) {
    if (verbose) {
      self$currentVariant(envir)
      cat("\nConsensus sequence:\n")
      self$showConsensus(envir)
    }
  }

  invisible(self)
})


## private$step_through() ####
MergeEnv_$set("private", "step_through", function(envir) {
  if (!hlatools::hasNext(envir$POSit)) {
    return(FALSE)
  }
  envir$pos <- iterators::nextElem(envir$POSit) + offset(envir$SR)
  # debug
  envir$balance_upper_confint
  rs <- disambiguate_variant(yield(envir), threshold = self$threshold)
  update(envir) <- rs
  TRUE
})


## self$showConsensus() ####
MergeEnv_$set("public", "showConsensus", function(envir, pos, left = 6, right = left, offset = 0) {
  # debug
  # envir <- self$a
  #self$a$init()
  if (is.null(envir$init)) {
    cat("Haplotype", envir$haplotype, "not initialised.")
    return(invisible(NULL))
  }
  if (missing(pos)) {
    pos <- envir$pos
  }
  min <- minimum(pos + offset - left, 1)
  lr <- envir$LR[min:(pos + offset + right), , drop = FALSE]
  sr <- envir$SR[min:(pos + offset + right), , drop = FALSE]
  ## Conseq
  lcs <- tolower(make_ambig_consensus_(lr, threshold = self$threshold, exclude_gaps = FALSE, as_string = TRUE))
  substr(lcs, left + 1, left + 1) <- toupper(substr(lcs, left + 1, left + 1))
  scs <- tolower(make_ambig_consensus_(sr, threshold = self$threshold, exclude_gaps = FALSE, as_string = TRUE))
  substr(scs, left + 1, left + 1) <- toupper(substr(scs, left + 1, left + 1))
  show <- sprintf(" Haplotype %s [%s] \nlr: %s\nsr: %s\n\n", envir$haplotype, pos, lcs, scs)
  cat(show)
})

## self$showMatrix() ####
MergeEnv_$set("public", "showMatrix", function(envir, pos, left = 6, right = left, offset = 0) {
  if (is.null(envir$init)) {
    cat("Haplotype", envir$haplotype, "not initialised.")
    return(invisible(NULL))
  }
  if (missing(pos)) {
    pos <- envir$pos
  }
  min <- minimum(pos + offset - left, 1)
  lr <- envir$LR[min:(pos + offset + right), , drop = FALSE]
  sr <- envir$SR[min:(pos + offset + right), , drop = FALSE]
  lcs <- make_ambig_consensus_(lr, threshold = self$threshold, exclude_gaps = FALSE, as_string = TRUE)
  scs <- make_ambig_consensus_(sr, threshold = self$threshold, exclude_gaps = FALSE, as_string = TRUE)
  cat("Haplotype ", envir$haplotype, "\nLong read map position [", pos + offset, "] Consensus [", lcs, "]\n")
  print(lr, n = NROW(lr), noHead = TRUE, transpose = TRUE)
  cat("Short read map position [", pos + offset, "] Consensus [", scs, "]\n")
  print(sr, n = NROW(sr), noHead = TRUE, transpose = TRUE)
})

## self$export() ####
MergeEnv_$set("public", "export", function() {

  cons <- structure(
    c(
      foreach (hptype = self$x$getHapTypes(),
                 .final = function(x) setNames(x, self$x$getHapTypes())) %do% {
        envir <- self$hptypes[[hptype]]
        list(
        matrix       = envir$SR,
        variants     = compact(envir$variants),
        phasemat     = envir$phasemat,
        phasebreaks  = envir$phasebreaks
      )},
      seq = list(foreach(hptype = self$x$getHapTypes(),
                            .final = function(x) setNames(x, self$x$getHapTypes())) %do% {
        sr <- self$hptypes[[hptype]]$SR
        cseq <- conseq(sr, paste0("HAP", hptype), "ambig", exclude_gaps = TRUE, threshold = self$threshold)
        metadata(cseq) <- list()
        cseq
      })
    ),
    class = c("ConsList", "list")
    )

  self$x$setConsensus(cons)
  return(invisible(self$x))
})


# Helpers -----------------------------------------------------------------

expand_longread_consmat <- function(lrm, srm) {
  m <- cbind(as.matrix(lrm), `+` = 0)
  if (length(ins(srm)) > 0) {
    insert <- matrix(c(0, 0, 0, 0, median(rowSums(m)), 0), ncol = 6)
    my_ins <- sort(ins(srm))
    my_ins <- my_ins[my_ins < nrow(lrm)]
    INSit <- hlatools::ihasNext(iter(my_ins))
    while (hlatools::hasNext(INSit)) {
      i <- iterators::nextElem(INSit)
      m <- rbind(m[1:(i - 1), ], insert, m[i:NROW(m), ])
    }
    stopifnot(NROW(m) == NROW(srm))
    dimnames(m) <- dimnames(srm)
  }
  lrm <- consmat(m, freq = FALSE)
  lrm
}

