#' @export
#debug
#x <- readDR2S("~/bioinf/DR2S/KIR/20171130/3DL3/out/7203477/")
# x <- dr2s
# threshold <- x$getThreshold()
# cache = TRUE
# library(foreach)
# library(futile.logger)
#x <- dr2s
polish.DR2S <- function(x,
                        threshold = NULL,
                        checkHpCount = TRUE,
                        cache = TRUE) {
  flog.info("Step 5: polish ...", name = "info")
  ## Set this threshold to double the usual.
  ## Probably better to do this only for lr
  if (is.null(threshold)) {
    threshold <- max(x$getThreshold(), 0.2)
  }

  ## Check if reporting is finished and exit safely for downstream analysis
  if (x$getReportStatus()) {
    currentCall <- strsplit1(deparse(sys.call()), "\\.")[1]
    flog.info("%s: Reporting already done! Nothing to do." %<<%
                " Exit safely for downstream analysis.",
              currentCall, name = "info")
    return(invisible(x))
  }
  assert_that(
    is(x$mapFinal, "mapFinal"),
    all(unlist(foreach(i = x$mapFinal$pileup) %do% {
      is(i$consmat, "consmat")
    }))
  )

  args <- x$getOpts("polish")
  if (!is.null(args)) {
    env  <- environment()
    list2env(args, envir = env)
  }

  # get read types and haplotypes
  rtypes  <- names(x$mapFinal$pileup)
  hptypes <- x$getHapTypes()
  menv <- MergeEnv(x, threshold)
  for (hptype in hptypes) {
    menv$init(hptype)
    menv$walk(hptype)
  }
  rs <- menv$export()
  ## Get variants
  vars <- .getVariants(x = rs)

  ## Check homopolymer count; Only check if the count is found in both
  if (checkHpCount) {
    checkHomoPolymerCount(rs)
    vars <- rbind(vars, foreach(hptype = hptypes, .combine = rbind) %do% {
      if (hptype %in% names(rs$mapFinal$homopolymers)) {
        seq <- rs$consensus$seq[[hptype]]
        seqrle <- .seq2rle(seq)
        n <- seqrle$length[seqrle$length > 10]
        modeN <- sort(rs$mapFinal$homopolymers[[hptype]])
        names(n) <- names(modeN)
        if (!all(modeN == n)) {
          missingN <- modeN[which(!n == modeN)]
          varsHP <- tibble::tibble(haplotype = hptype,
                                   pos       = names(missingN),
                                   ref       = "",
                                   alt       = "",
                                   warning   = sprintf(paste(
                                     "Homopolymer at position %s should be of", 
                                     "length %s, but is %s"),
                                     names(missingN), missingN, 
                                     n[names(missingN)]),
                                   refSR     = "",
                                   altSR     = "",
                                   refLR     = "",
                                   altLR     = "")
          varsHP
        }
      }
    })
  }

  rs$consensus$variants = dplyr::arrange(vars, .data$pos, .data$haplotype)

  if (cache)
    rs$cache()

  invisible(rs)
}


# Helpers -----------------------------------------------------------------



.getVariants <- function(x) {
  hptypes <- x$getHapTypes()
  hvars <- lapply(hptypes, function(t) x$consensus[[t]]$variants)
  names(hvars) <- hptypes

  if (all(sapply(hvars, function(x) length(x) == 0))) {
    return(
      dplyr::data_frame(
        haplotype = character(0), pos = integer(0), ref = character(0),
        alt = character(0), warning = character(0), refSR = character(0),
        altSR = character(0), refLR = character(0), altLR = character(0)
      )
    )
  }

  df <- do.call("rbind",
                lapply(hptypes, function(h)
                  do.call("rbind", lapply(hvars[[h]],
                                          .extractVariant_, h = h))))
  df
}
.extractVariant_ <- function(v, h) {
  data.frame(
    haplotype = h,
    pos = attr(v, "position") %||% NA,
    ref = v[[1]] %||% NA,
    alt = v[[2]] %||% NA,
    warning = attr(v, "warning") %||% NA,
    refSR = attr(v, "srBases")[[1]] %||% NA,
    altSR = attr(v, "srBases")[[2]] %||% NA,
    refLR = attr(v, "lrBases")[[1]] %||% NA,
    altLR = attr(v, "lrBases")[[2]] %||% NA,
    stringsAsFactors = FALSE
  )
}
