#' @export
#debug
#x <- readDR2S("~/bioinf/DR2S/KIR/20171130/3DL3/out/7203477/")
# x <- dr2s
# threshold <- x$getThreshold()
# lower_limit <- 0.6
# cache = TRUE
# library(foreach)
# library(futile.logger)
#x <- self
polish.DR2S <- function(x,
                        threshold = x$getThreshold(),
                        lower_limit = 0.80,
                        check_hp_count = TRUE,
                        cache = TRUE) {
  flog.info("Step 5: polish ...", name = "info")

  ## Check if reporting is finished and exit safely for downstream analysis
  if (x$getReportStatus()) {
    currentCall <- strsplit1(deparse(sys.call()), "\\.")[1]
    flog.info(paste0("%s: Reporting already done! Nothing to do.",
                     " Exit safely for downstream analysis."),
              currentCall, name = "info")
    return(invisible(x))
  }

  assertthat::assert_that(
    is(x$mapFinal, "mapFinal"),
    all(unlist(foreach(rtype = x$mapFinal$pileup) %do% {
      is(rtype$consmat, "consmat")
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

  ## Problematic Variants
  vars <- get_problematic_variants(x = rs, lower_limit = 0.6)
  vars <- dplyr::ungroup(vars)

  ## Check homopolymer count; Only check if the count is found in both
  if (check_hp_count) {
    checkHomoPolymerCount(rs)
    vars <- rbind(vars, foreach(hptype = hptypes, .combine = rbind) %do% {
      if (hptype %in% names(rs$mapFinal$homopolymers)) {
        seq <- rs$consensus$seq[[hptype]]
        seqrle <- .seq2rle(seq)
        n <- seqrle$length[seqrle$length > 10]
        modeN <- sort(rs$mapFinal$homopolymers[[hptype]])
        if (!all(modeN == n)) {
          missingN <- modeN[which(!n == modeN)]
          varsHP <- tibble::tibble(haplotype = hptype,
                                   pos       = names(missingN),
                                   ref       = "",
                                   alt       = "",
                                   freq      = "",
                                   lower     = "",
                                   warning   = sprintf(
                                     "Homopolymer at position %s should be %s",
                                                     names(missingN), missingN))
          varsHP
        }
      }
    })
  }

  rs$consensus$problematic_variants = dplyr::arrange(vars, pos, haplotype)

  if (cache)
    rs$cache()

  invisible(rs)
}


# Helpers -----------------------------------------------------------------

# x <- rs

get_problematic_variants <- function(x, lower_limit) {
  stopifnot(is(x, "DR2S"))
  vars <- collect_variants(x)
  vars <- dplyr::group_by(vars, haplotype)
  vars <- dplyr::filter(vars, lower < lower_limit | nzchar(warning))
  dplyr::mutate(vars, freq = round(freq, 3), lower = round(lower, 3))
}


collect_variants <- function(x) {
  hptypes <- x$getHapTypes()
  hvars <- lapply(hptypes, function(t) x$consensus[[t]]$variants)
  names(hvars) <- hptypes

  if (all(sapply(hvars, function(x) length(x) == 0))) {
    return(
      dplyr::data_frame(
        haplotype = character(0), pos = integer(0), ref = character(0),
        alt = character(0), freq = double(0), lower = double(0),
        warning = character(0)
      )
    )
  }

  df <- do.call("rbind",
                lapply(hptypes, function(h)
                  do.call("rbind", lapply(hvars[[h]],
                                          extract_variant_, h = h))))
  df
}

extract_variant_ <- function(v, h) {
  data.frame(
    haplotype = h,
    pos = attr(v, "position") %||% NA,
    ref = v[["ref"]] %||% NA,
    alt = v[["alt"]] %||% NA,
    freq = attr(v, "proportion") %||% NA,
    lower = (attr(v, "proportion") - attr(v, "margin_of_error")) %||% NA,
    warning = attr(v, "warning") %||% NA,
    stringsAsFactors = FALSE
  )
}
