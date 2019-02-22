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
  assert_that(is(x, "DR2S"))
  if (missing(threshold)) {
    threshold <- x$getThreshold()
  }
  MergeEnv_$new(x, threshold)
}

#' Class \code{"MergeEnv"}
#'
#' @docType class
#' @usage MergeEnv(x, threshold)
#' @field threshold \code{[numeric]}; when do we call a variant a variant.
#' @field x \code{[\link[=DR2S_]{DR2S}]}; the original \code{DR2S} object.
#' @field threshold When do we call a variant a variant.
#' @keywords data internal
#' @return Object of \code{\link[R6]{R6Class}} representing a MergeEnv.
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
      assert_that(is(x, "DR2S"))
      self$hptypes <- foreach(hptype = x$getHapTypes(),
                              .final = function(h)
                                setNames(h, x$getHapTypes())) %do% {
       structure(as.environment(list(haplotype = hptype)), class = "HapEnv")
      }
      #threshold <- max(c(threshold, 0.2))
      # threshold <- max(c(threshold, 0.3))
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
      if (!is.null(envir$currentVariant)) {
        print(envir$currentVariant)
      } else NULL
    },
    print = function(left = 3, right = left, ...) {
      for (hapEnv in self$hptypes) {
        self$showMatrix(envir = hapEnv, left = left, right = right)
        cat("\n")
      }
    }
  )
)

## self$init() ####
MergeEnv_$set("public", "init", function(hapEnv) {
  hapEnv <- match.arg(hapEnv, self$x$getHapTypes())
  envir <- self$hptypes[[hapEnv]]

  if (self$x$hasShortreads()) {
    lr <- consmat(self$x$remap$LR[[hapEnv]]$pileup, prob = FALSE)
    sr <- consmat(self$x$remap$SR[[hapEnv]]$pileup, prob = FALSE)
    rs <- .equaliseConsmat(lrm = lr, srm = sr)
    envir$LR <- rs$lrm
    envir$SR <- rs$srm
    referencePath <- self$x$absPath(self$x$remap$SR[[hapEnv]]$refpath)
    reference <- Biostrings::readDNAStringSet(referencePath)
    envir$ref <- unname(strsplit1(as.character(reference), ""))
  } else {
    envir$LR <- consmat(self$x$remap$LR[[hapEnv]]$pileup, prob = FALSE)
    envir$SR <- NULL
    referencePath <- self$x$absPath(self$x$remap$LR[[hapEnv]]$refpath)
    reference <- Biostrings::readDNAStringSet(referencePath)
    envir$ref <- unname(strsplit1(as.character(reference), ""))
  }
  apos <- foreach(rt = c("LR", "SR"), .combine = c) %do% {
    ## positions not matching the consensus
    dism <- .noRefMatch(envir[[rt]], envir[["ref"]])
    ## ambiguous positions
    amb <- .ambiguousPositions(envir[[rt]], self$threshold, TRUE)
    sort(c(amb, dism))
  }
  apos <- unique(sort(apos))
  envir$POSit = itertools::ihasNext(iterators::iter(apos))
  envir$variants = list()
  envir$pos = 1L
  envir$currentVariant = NULL
  envir$init = TRUE
})

## self$walkOne() ####
MergeEnv_$set("public", "walkOne", function(hapEnv, verbose = FALSE) {
  hapEnv <- match.arg(hapEnv, self$x$getHapTypes())
  if (!self$isInitialised(hapEnv)) {
    self$init(hapEnv)
  }

  envir <- self$getHapEnv(hapEnv)
  private$stepThrough(envir)
  if (verbose) {
    self$currentVariant(envir)
    cat("\nConsensus sequence:\n")
    self$showConsensus(envir)
  }
  invisible(self)
})

## self$walk() ####
MergeEnv_$set("public", "walk", function(hp, verbose = FALSE) {
  hp <- match.arg(hp, self$x$getHapTypes())
  if (!self$isInitialised(hp)) {
    self$init(hp)
  }
  envir <- self$getHapEnv(hp)
  # while (stepThrough(envir)) {
  while (private$stepThrough(envir)) {
    if (verbose) {
      self$currentVariant(envir)
      cat("\nConsensus sequence:\n")
      self$showConsensus(envir)
    }
  }

  invisible(self)
})

## private$stepThrough() ####
MergeEnv_$set("private", "stepThrough", function(envir) {
  # stepThrough <- function(envir) {
  if (!itertools::hasNext(envir$POSit)) {
    return(FALSE)
  }
  envir$pos <- ifelse(!is.null(envir$SR),
                      iterators::nextElem(envir$POSit) + offsetBases(envir$SR),
                      iterators::nextElem(envir$POSit) + offsetBases(envir$LR))
  # message(envir$pos)
  # envir$pos <- 962
  # envir$pos <- 4136
  # p  <- envir$pos
  # x  <- yield(envir)
  rs <- disambiguateVariant(x = yield(envir), threshold = self$threshold)
  .update(envir) <- rs
  TRUE
})

## self$showConsensus() ####
MergeEnv_$set("public", "showConsensus", function(envir,
                                                  pos, left = 6, right = left,
                                                  offsetBases = 0) {
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
  min <- min(pos + offsetBases - left, 1)
  lr <- envir$LR[min:(pos + offsetBases + right), , drop = FALSE]
  sr <- envir$SR[min:(pos + offsetBases + right), , drop = FALSE]
  ## Conseq
  lcs <- tolower(.makeAmbigConsensus_(lr, threshold = self$threshold,
                                      suppressAllGaps = FALSE, asString = TRUE))
  ## here
  substr(lcs, left + 1, left + 1) <- toupper(substr(lcs, left + 1, left + 1))
  scs <- tolower(.makeAmbigConsensus_(sr, threshold = self$threshold,
                                      suppressAllGaps = FALSE, asString = TRUE))
  substr(scs, left + 1, left + 1) <- toupper(substr(scs, left + 1, left + 1))
  show <- sprintf(" Haplotype %s [%s] \nlr: %s\nsr: %s\n\n",
                  envir$haplotype, pos, lcs, scs)
  cat(show)
})

## self$showMatrix() ####
MergeEnv_$set("public", "showMatrix", function(envir, pos, left = 6,
                                               right = left, offsetBases = 0) {
  if (is.null(envir$init)) {
    cat("Haplotype", envir$haplotype, "not initialised.")
    return(invisible(NULL))
  }
  if (missing(pos)) {
    pos <- envir$pos
  }
  min <- min(max(pos + offsetBases - left, 1), 1)
  lr <- envir$LR[min:(pos + offsetBases + right), , drop = FALSE]
  sr <- envir$SR[min:(pos + offsetBases + right), , drop = FALSE]
  lcs <- .makeAmbigConsensus_(lr, threshold = self$threshold,
                              suppressAllGaps = FALSE, asString = TRUE)
  if (!is.null(sr)) {
    scs <- .makeAmbigConsensus_(sr, threshold = self$threshold,
                                suppressAllGaps = FALSE, asString = TRUE)
  }
  cat("Haplotype ", envir$haplotype,
      "\nLong read map position [", pos + offsetBases,
      "] Consensus [", lcs, "]\n")
  print(lr, n = NROW(lr), noHead = TRUE, transpose = TRUE)
  if (!is.null(sr)) {
    cat("Short read map position [", pos + offsetBases, "] Consensus [",
        scs, "]\n")
    print(sr, n = NROW(sr), noHead = TRUE, transpose = TRUE)
  }
})

## self$export() ####
MergeEnv_$set("public", "export", function() {
  cons <- structure(
    c(
      foreach(hp = self$x$getHapTypes(),
              .final = function(x) stats::setNames(x, self$x$getHapTypes())) %do% {
                envir <- self$hptypes[[hp]]
                list(
                  matrix = if (!is.null(envir$SR)) envir$SR else envir$LR,
                  variants = compact(envir$variants)
                )
              },
      ## consensus for checking with ambigs
      seq = list(foreach(hp = self$x$getHapTypes(),
                         .final = function(x)
                           stats::setNames(x, self$x$getHapTypes())) %do% {
                             cmat <- if (!is.null(self$hptypes[[hp]]$SR)) {
                               self$hptypes[[hp]][["SR"]]
                             } else {
                               consmat(self$x$mapFinal$LR[[hp]]$pileup)
                             }
                             cseq <- conseq(x = cmat, name = "hap" %<<% hp, type = "ambig", suppressAllGaps = TRUE, threshold = self$threshold)
                             S4Vectors::metadata(cseq) <- list()
                             cseq
                           }),
      ## consensus for remapping without ambiguities
      noAmbig = list(foreach(hp = self$x$getHapTypes(),
                             .final = function(x)
                               setNames(x, self$x$getHapTypes())) %do% {
                                 cmat <- if (!is.null(self$hptypes[[hp]]$SR)) {
                                   self$hptypes[[hp]]$SR
                                 } else {
                                   consmat(self$x$mapFinal$LR[[hp]]$pileup)
                                 }
                                 cseq <- conseq(cmat, "hap" %<<% hp, "prob", suppressAllGaps = TRUE)
                                 S4Vectors::metadata(cseq) <- list()
                                 cseq
                               })
    ),
    class = c("ConsList", "list")
  )
  self$x$setConsensus(cons)
  return(invisible(self$x))
})


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

