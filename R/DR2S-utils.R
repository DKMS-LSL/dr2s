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
  igv_batch_files <- file.path(basedir, paste0("igv", c("A", "B"), ".txt"), fsep = )
  igv1 <- igv_batch_files[1]
  igv2 <- igv_batch_files[2]

  refA   <- file.path(basedir, "A", basename(x$mapFinal$aref))
  bamALR <- file.path(basedir, "merged", basename(x$mapFinal$bamfile$LRA))
  bamASR <- file.path(basedir, "merged", basename(x$mapFinal$bamfile$SRA))
  refB   <- file.path(basedir, "B", basename(x$mapFinal$bref))
  bamBLR <- file.path(basedir, "merged", basename(x$mapFinal$bamfile$LRB))
  bamBSR <- file.path(basedir, "merged", basename(x$mapFinal$bamfile$SRB))
  chrA <- sub(">", "", readLines(refA, 1))
  chrB <- sub(">", "", readLines(refB, 1))

  if (missing(position)) {
    position <- as.integer(x$consensus$problematic_variants[1, "pos"])
  }

  batchA <- list(
    new = "",
    genome = refA,
    load   = bamALR,
    load   = bamASR,
    goto   = paste0(chrA, ":", position - 50, "-", position + 50),
    sort   = "position"
  )

  batchB <- list(
    new = "",
    genome = refB,
    load   = bamBLR,
    load   = bamBSR,
    goto   = paste0(chrB, ":", position - 50, "-", position + 50),
    sort   = "position"
  )

  if (.Platform$OS.type != "windows") {
    cat(sprintf("%s %s", names(batchA), batchA), sep = "\n", file = igv1)
    system(paste0(shQuote(exe_path), " -b ", shQuote(igv1), " &"), ...)
    cat(sprintf("%s %s", names(batchB), batchB), sep = "\n", file = igv2)
    system(paste0(shQuote(exe_path), " -b ", shQuote(igv2), " &"), ...)
  } else if (.Platform$OS.type == "windows") {
    cat(sprintf("%s %s", names(batchA), batchA), sep = "\n", file = igv1)
    cat(sprintf("%s %s", names(batchB), batchB), sep = "\n", file = igv2)
    msg <- sprintf(
      "\nOpen two instances of the command prompt and paste in the commands:\n%s\n%s\n",
      paste0(shQuote(exe_path), " -b ", shQuote(igv1)),
      paste0(shQuote(exe_path), " -b ", shQuote(igv2))
    )
    message(msg)
  }
}
