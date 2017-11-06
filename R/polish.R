#' @export

#debug
#x <- dpb1_3.r
# threshold <- x$getThreshold()
# lower_limit <- 0.6
# cache = TRUE
# library(foreach)
#x <- a
polish.DR2S <- function(x, threshold = x$getThreshold(), lower_limit = 0.80, cache = TRUE) {

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
  rtypes <- names(x$mapFinal$pileup)
  hptypes <- x$getHapTypes()
  menv <- MergeEnv(x, threshold)
  for (hptype in hptypes){
    menv$init(hptype)
    menv$walk(hptype)
  }

  ## Short read phasing
  ## phasing.R
  for (hptype in hptypes){
    menv$hptypes[[hptype]]$phasemat    <- phasematrix(compact(menv$hptypes[[hptype]]$variants))
    menv$hptypes[[hptype]]$phasebreaks <- phasebreaks(menv$hptypes[[hptype]]$phasemat)
  }

  ## phasing plots
  for (hptype in hptypes){
    if (!is.null(menv$hptypes[[hptype]]$phasemat)){
      file <- file.path(x$getOutdir(), paste("plot4.phasing", hptype, x$getLrdType(), x$getMapper(), "pdf", sep = "."))
      pdf(file, width = 16, height = 10, onefile = TRUE)
      p <- phaseplot(menv$hptypes[[hptype]]$phasemat) +
        ggplot2::ggtitle(paste0("Short read phasing A ", hptype))
      print(p)
      dev.off()
    }
  }

  rs <- menv$export()

  ## Problematic Variants
  vars <- get_problematic_variants(x = rs, lower_limit = .6)
  vars <- dplyr::ungroup(vars)
  rs$consensus$problematic_variants = dplyr::arrange(vars, pos, haplotype)
  
  if (cache)
    rs$cache()
  rs
}


# Helpers -----------------------------------------------------------------

# x <- rs

get_problematic_variants <- function(x, lower_limit) {
  stopifnot(is(x, "DR2S"))
  vars <- collect_variants(x = x)
  vars <- dplyr::group_by(vars, haplotype)
  vars <- dplyr::filter(vars, lower < lower_limit | nzchar(warning))
  dplyr::mutate(vars, freq = round(freq, 3), lower = round(lower, 3))
}


collect_variants <- function(x) {
  hptypes <- x$getHapTypes()
  hvars <- lapply(hptypes, function(t) x$consensus[[t]]$variants)
  names(hvars) <- hptypes


  if (all(sapply(hvars, function(x) length(x) == 0))){
    return(
      dplyr::data_frame(
        haplotype = character(0), pos = integer(0), ref = character(0),
        alt = character(0), freq = double(0), lower = double(0), warning = character(0)
      )
    )
  }

  phasebreaks <- lapply(hptypes, function(t) dplyr::filter(x$consensus[[t]]$phasebreaks, breakpos)$posx)
  names(phasebreaks) <- hptypes

  df <- do.call("rbind",
                  lapply(hptypes, function(h) do.call("rbind", lapply(hvars[[h]],  extract_variant_, h = h))))

  df <- df %>%
      dplyr::group_by(haplotype) %>%
      dplyr::mutate(warning = ifelse(pos  %in% unlist(phasebreaks[haplotype]),
                                     sub("^\\|", "", warning %<<% "|Phasebreak in short reads"),
                                     warning))
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
