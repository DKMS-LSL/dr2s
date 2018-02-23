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
  conf_log(outdir = outdir, logName = "info")
  if (length(rs) == 1) {
    rs[[1]]
  } else rs
}


finish_cn1 <- function(x){
  flog.info("Set only allele to A", name = "info")
  x$setHapTypes(c("A"))

  flog.info("Write mapInit data to mapFinal", name = "info")
  x$mapFinal = structure(
    list(
      dir     = x$getOutdir(),
      sreads  = list(A = x$getSrdDir()),
      lreads  = list(A = x$getLrdDir()),
      ref     = list(A = x$mapInit$SR2$ref),
      bamfile = list(),
      pileup  = list(),
      tag     = list(),
      seqpath = list()
    ), class = c("mapFinal", "list")
  )
  if (!is.null(x$mapInit$SR1)) {
    # Write shortread data to mapFinal
    mapgroupSR <- "SRA"
    x$mapFinal$pileup[[mapgroupSR]] = x$mapInit$SR2$pileup
    x$mapFinal$tag[[mapgroupSR]] = x$mapInit$SR2$tag
  }

  # Write longread data to mapFinal
  mapgroupLR <- "LRA"
  x$mapFinal$pileup[[mapgroupLR]] = x$mapInit$pileup
  x$mapFinal$tag[[mapgroupLR]] = x$mapInit$tag

  flog.info("Get latest consensus from last mapping ...", name = "info")

  if (!is.null(x$mapInit$SR1)) {
    cseq <- conseq(x$mapInit$SR2$pileup$consmat, "mapFinalA", "ambig", exclude_gaps = TRUE, threshold = x$getThreshold())
  } else {
    cseq <- conseq(x$mapInit$pileup$consmat, "mapFinalA", "ambig", exclude_gaps = TRUE, threshold = x$getThreshold())
  }
  x$mapFinal$seq$A <- cseq

  flog.info("Polish and look for inconsistencies between shortreads and longreads ...", name = "info")
  polish(x)
  flog.info("Report consensus sequence and potential problematic variants", name = "info")
  report(x)

  return(x)
}

#' @export
# x <- dedk.report
#position = 5000
run_igv <- function(x, position, map = "mapFinal", open_now = TRUE, ...) {
  assertthat::assert_that(
    is(x, "DR2S")
  )
  if (!nzchar(exe_path <- normalizePath(Sys.which("igv"), mustWork = FALSE))) {
    stop("No IGV executable found on the PATH", call. = FALSE)
  }
  map <- match.arg(map, c("mapInit",
                          "mapFinal",
                          "mapIter",
                          "refine"))


  basedir <- normalizePath(x$getOutdir(), mustWork = TRUE)
  igvdir <- file.path(basedir, ".igv")
  dir_create_if_not_exists(igvdir)
  if (.Platform$OS.type == "windows") {
    fsep <- "\\"
  } else {
    fsep <- "/"
  }
  if (missing(position)) {
    position <- ifelse(NROW(x$consensus$problematic_variants) == 0,
                       50,
                       as.integer(x$consensus$problematic_variants[1, "pos"]))
  }
  igvConfigs <- list()
  # hp = "A"
  haptypes <- ifelse(map == "mapInit", "Init", x$getHapTypes())
  for (hp in haptypes){
    igv <- file.path(igvdir, paste0("igv", hp, map, ".xml"))
    if (map == "mapFinal") {
      ref   <- file.path(hp, basename(x$mapIter[[as.character(x$getIterations())]][[hp]]$seqpath))
      bamLR <- file.path("final", basename(x$mapFinal$bamfile[[paste0("LR", hp)]]))
      if (!is.null(x$mapFinal$sreads[[hp]]))
        bamSR <- file.path("final", basename(x$mapFinal$bamfile[[paste0("SR", hp)]]))
    } else if (map == "refine") {
      if (!is.null(x$consensus$refine$ref[[hp]])) {
        ref   <- x$consensus$refine$ref[[hp]]
        bamLR <- x$consensus$refine$bamfile[[paste0("LR", hp)]]
        if (!is.null(x$mapFinal$sreads[[hp]]))
          bamSR <- x$consensus$refine$bamfile[[paste0("SR", hp)]]
      } else {
        ref   <- file.path(hp, basename(x$mapIter[[as.character(x$getIterations())]][[hp]]$seqpath))
        bamLR <- file.path("final", basename(x$mapFinal$bamfile[[paste0("LR", hp)]]))
        if (!is.null(x$mapFinal$sreads[[hp]]))
          bamSR <- file.path("final", basename(x$mapFinal$bamfile[[paste0("SR", hp)]]))
      }
    } else if (map == "mapInit") {
      if (x$getPartSR()) {
        ref <- unname(x$mapInit$SR2$seqpath)
        bamSR <- x$mapInit$SR2$bamfile
      } else {
        ref <- x$getRefSeq()
      }
      bamLR <- x$mapInit$bamfile

    } else if (map == "mapIter") {
      ref <- x$mapIter[[x$getIterations()]][[hp]]$seqpath
      bamLR <- x$mapIter[[x$getIterations()]][[hp]]$bamfile
    }
    chr <- strsplit(sub(">", "", readLines(file.path(basedir, ref), 1)), "\\s+")[[1]][1]
    locus <- paste0(chr, ":", min(c((abs(position - 50)),0)), "-", position + 50)

    xml <- XML::xmlTree()
    suppressWarnings(xml$addTag("Global",
                                attrs = c(genome = file.path("..", ref),
                                          locus = locus),
                                close = FALSE))
    xml$addTag("Resources", close = FALSE)
    xml$addTag("Resource", attrs = c(path = file.path("..", bamLR)))
    if (exists("bamSR"))
      xml$addTag("Resource", attrs = c(path = file.path("..", bamSR)))
    xml$closeTag()
    xml$closeTag()
    XML::saveXML(xml, file = igv)

    igvConfigs[[hp]] <- x$relPath(igv)
  }
  igvCommand <- file.path(basedir, paste0("IGV_", map, ".sh"))
  cmds <- paste0("igv ", igvConfigs, " &")
  write(cmds, igvCommand)
  Sys.chmod(igvCommand, mode = "775")

  if (open_now) {
    if (.Platform$OS.type != "windows") {
      cwd <- getwd()
      setwd(basedir)
      system(igvCommand)
      setwd(cwd)
    } else if (.Platform$OS.type == "windows") {
      msg <- sprintf(
        "\nOpen an instances of the command prompt and paste in the command:\n%s\n",
        paste0(shQuote(exe_path), " ", shQuote(file.path(basedir, igv)))
      )
      message(msg)
    }
  }
}
