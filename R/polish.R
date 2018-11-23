#' @export
polish.DR2S <- function(x, ...) {
  ## If reporting is already done exit safely
  if (.checkReportStatus(x)) return(invisible(x))
  assert_that(x$hasMapFinal())

  ## Collect start time for polish runstats
  start.time <- Sys.time()

  flog.info("# polish", name = "info")

  ## Get global arguments
  ## Set this threshold to double the usual.
  ## Probably better to do this only for lr
  threshold <- max(x$getThreshold(), 0.2)

  ## Export polish config to function environment
  args <- x$getOpts("polish")
  list2env(args, envir = environment())
  assert_that(
    is.double(threshold),
    exists("checkHpCount") && is.logical(checkHpCount),
    exists("hpCount") && is.numeric(hpCount)
  )

  # get read types and haplotypes
  # rdtypes <- names(x$mapFinal) ## LR, SR
  hptypes <- x$getHapTypes() ## A, B, C, ...
  menv <- MergeEnv(x, threshold)
  ## hp = "D"
  for (hp in hptypes) {
    menv$init(hp)
    menv$walk(hp)
  }
  rs <- menv$export()
  ## Get variants
  vars <- .getVariants(x = rs)

  ## Check homopolymer count; Only check if the count is found in both
  if (checkHpCount) {
    checkHomopolymerCount(rs, hpCount = hpCount)
    vars <- rbind(vars, foreach(hp = hptypes, .combine = rbind) %do% {
      if (hp %in% names(rs$consensus$homopolymers)) {
        seq <- rs$consensus$seq[[hp]]
        seqrle <- .seq2rle(seq)
        # n <- seqrle$length[seqrle$length >= 10] %||% 0
        n <- which(seqrle$length > 8)
        nCount <- seqrle$length[n]
        origPosition <- vapply(n, function(ni) sum(seqrle$lengths[1:((ni) - 1)]), FUN.VALUE = integer(1))
        modeN <- sort(rs$consensus$homopolymers[[hp]])
        #names(n) <- names(modeN)
        if (!all(names(modeN) %in% vapply(origPosition, function(ni) (ni - 5):(ni + 5), FUN.VALUE = integer(11)))) {
          missingN <- modeN[which(!n %in% modeN)]
          varsHP <- tibble::tibble(haplotype = hp,
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
          return(varsHP)
        }
      }
      tibble::tibble(haplotype = "",
                     pos       = "",
                     ref       = "",
                     alt       = "",
                     warning   = "",
                     refSR     = "",
                     altSR     = "",
                     refLR     = "",
                     altLR     = "")
    })
  }

  rs$consensus$variants = dplyr::arrange(vars, .data$pos, .data$haplotype)

  ## set polish runstats
  .setRunstats(x, "polish",
               list(Runtime = format(Sys.time() - start.time)))

  return(invisible(rs))
}


# Helpers -----------------------------------------------------------------


.getVariants <- function(x) {
  hptypes <- x$getHapTypes()
  hvars <- set_names(lapply(hptypes, function(t) x$consensus[[t]]$variants),
                     hptypes)

  if (all(vapply(hvars, function(x) length(x) == 0, logical(1)))) {
    return(
      tibble::data_frame(
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
