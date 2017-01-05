#' @export
polish.DR2S <- function(x, threshold = x$getThreshold(), lower_limit = 0.60, cache = TRUE) {
  assertthat::assert_that(
    is(x$map3, "map3"),
    is(x$map3$pileup$LRA$consmat, "consmat"),
    is(x$map3$pileup$SRA$consmat, "consmat"),
    is(x$map3$pileup$LRB$consmat, "consmat"),
    is(x$map3$pileup$SRB$consmat, "consmat")
  )

  args <- x$getOpts("polish")
  if (!is.null(args)) {
    env  <- environment()
    list2env(args, envir = env)
  }

  menv <- MergeEnv(x, threshold)
  menv$init("A")
  menv$walk("A")
  menv$init("B")
  menv$walk("B")

  ## Short read phasing
  ## phasing.R
  menv$a$phasemat    <- phasematrix(compact(menv$a$variants))
  menv$a$phasebreaks <- phasebreaks(menv$a$phasemat)
  menv$b$phasemat    <- phasematrix(compact(menv$b$variants))
  menv$b$phasebreaks <- phasebreaks(menv$b$phasemat)

  ## phasing plots
  if (!is.null(menv$a$phasemat)) {
    fileA <-file.path(x$getOutdir(), paste("plot4.phasing.A", x$getLrdType(), x$getMapper(), "pdf", sep = "."))
    pdf(fileA, width = 16, height = 10, onefile = TRUE)
    pa <- phaseplot(menv$a$phasemat) +
      ggplot2::ggtitle("Short read phasing A")
    print(pa)
    dev.off()
  }

  if (!is.null(menv$b$phasemat)) {
    fileB <- file.path(x$getOutdir(), paste("plot4.phasing.B", x$getLrdType(), x$getMapper(), "pdf", sep = "."))
    pdf(fileB, width = 16, height = 10, onefile = TRUE)
    pb <- phaseplot(menv$b$phasemat) +
      ggplot2::ggtitle("Short read phasing B")
    print(pb)
    dev.off()
  }

  rs <- menv$export()

  ## Problematic Variants
  vars <- get_problematic_variants(x = rs, lower_limit = lower_limit)
  vars <- dplyr::ungroup(vars)
  rs$consensus$problematic_variants = dplyr::arrange(vars, pos, haplotype)

  if (cache)
    rs$cache()
  rs
}


# Helpers -----------------------------------------------------------------


get_problematic_variants <- function(x, lower_limit) {
  stopifnot(is(x, "DR2S"))
  vars <- collect_variants(x = x$consensus)
  vars <- dplyr::group_by(vars, haplotype)
  vars <- dplyr::filter(vars, lower < lower_limit | nzchar(warning))
  dplyr::mutate(vars, freq = round(freq, 3), lower = round(lower, 3))
}

collect_variants <- function(x) {
  avars <- x$A$variants
  bvars <- x$B$variants

  if (length(avars) == 0 && length(bvars) == 0) {
    return(
      dplyr::data_frame(
        haplotype = character(0), pos = integer(0), ref = character(0),
        alt = character(0), freq = double(0), lower = double(0), warning = character(0)
      )
    )
  }

  dfa <- do.call("rbind", lapply(avars, extract_variant_, h = "A"))
  ## check for phasing problems
  aphasebreak <- dplyr::filter(x$A$phasebreaks, breakpos)$posx
  dfa <- dplyr::mutate(dfa, warning = ifelse(
    pos %in% aphasebreak,
    sub("^\\|", "", warning %<<% "|Phasebreak in short reads"),
    warning
  ))

  dfb <- do.call("rbind", lapply(bvars, extract_variant_, h = "B"))
  ## check for phasing problems
  bphasebreak <- dplyr::filter(x$B$phasebreaks, breakpos)$posx
  dfb <- dplyr::mutate(dfb, warning = ifelse(
    pos %in% bphasebreak,
    sub("^\\|", "", warning %<<% "|Phasebreak in short reads"),
    warning
  ))

  dplyr::as_data_frame(rbind(dfa, dfb))
}

extract_variant_ <- function(v, h) {
  data.frame(
    haplotype = h,
    pos = attr(v, "position"),
    ref = v[["ref"]],
    alt = v[["alt"]],
    freq = attr(v, "proportion"),
    lower = attr(v, "proportion") - attr(v, "margin_of_error"),
    warning = attr(v, "warning"),
    stringsAsFactors = FALSE
  )
}
