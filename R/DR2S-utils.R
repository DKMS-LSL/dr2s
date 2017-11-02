#' @export
read_dr2s <- function(path) {
  outdir <- normalizePath(path, mustWork = TRUE)
  rds_path <- dir(outdir, pattern = "^DR2S(YAML)?.*.rds", full.names = TRUE)
  if (length(rds_path) == 0) {
    warning("No DR2S analysis object found", call. = TRUE, immediate. = TRUE)
    return(NULL)
  }
  rs <- lapply(rds_path, readRDS)
  ## update outdir
  rs <- lapply(rs, function(x) {
    x$setConfig("outdir", outdir)
    x
  })
  if (length(rs) == 1) {
    rs[[1]]
  } else rs
}

#' @export
# x <- dedk.report
#position = 5000
run_igv <- function(x, position, ...) {
  assertthat::assert_that(
    is(x, "DR2S"),
    is(x$consensus, "ConsList")
  )

  if (!nzchar(exe_path <- normalizePath(Sys.which("igv"), mustWork = FALSE))) {
    stop("No IGV executable found on the PATH", call. = FALSE)
  }

  basedir <- normalizePath(x$getOutdir(), mustWork = TRUE)
  if (.Platform$OS.type == "windows") {
    fsep <- "\\"
  } else {
    fsep <- "/"
  }
  if (missing(position)) {
    position <- as.integer(x$consensus$problematic_variants[1, "pos"])
  }
  for (hp in x$getHapTypes()){
    igv <- file.path(basedir, paste0("igv", hp, ".txt"), fsep = )
    ref   <- file.path(basedir, hp, basename(x$mapIter[[as.character(x$getIterations())]][[hp]]$seqpath))
    bamLR <- file.path(basedir, "final", basename(x$mapFinal$bamfile[[paste0("LR", hp)]]))
    bamSR <- file.path(basedir, "final", basename(x$mapFinal$bamfile[[paste0("SR", hp)]]))
    chr <- sub(">", "", readLines(ref, 1))
    batch <- list(
      new = "",
      genome = ref,
      load   = bamLR,
      load   = bamSR,
      goto   = paste0(chr, ":", position - 50, "-", position + 50),
      sort   = "position"
    )
    if (.Platform$OS.type != "windows") {
      cat(sprintf("%s %s", names(batch), batch), sep = "\n", file = igv)
      system(paste0(shQuote(exe_path), " -g ", shQuote(ref), " ", shQuote(bamSR), " &"))
#      system(paste0(shQuote(exe_path), " -b ", shQuote(igv), " &"))#, ...)
    } else if (.Platform$OS.type == "windows") {
      cat(sprintf("%s %s", names(batch), batch), sep = "\n", file = igv)
      msg <- sprintf(
        "\nOpen an instances of the command prompt and paste in the command:\n%s\n",
        paste0(shQuote(exe_path), " -b ", shQuote(igv))
      )
      message(msg)
    }
  }


  # batchA <- list(
  #   new = "",
  #   genome = refA,
  #   load   = bamALR,
  #   load   = bamASR,
  #   goto   = paste0(chrA, ":", position - 50, "-", position + 50),
  #   sort   = "position"
  # )
  #
  # batchB <- list(
  #   new = "",
  #   genome = refB,
  #   load   = bamBLR,
  #   load   = bamBSR,
  #   goto   = paste0(chrB, ":", position - 50, "-", position + 50),
  #   sort   = "position"
  # )

}
