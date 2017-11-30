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
   # Write shortread data to mapFinal
   mapgroupSR <- "SRA"
   x$mapFinal$pileup[[mapgroupSR]] = x$mapInit$SR2$pileup
   x$mapFinal$tag[[mapgroupSR]] = x$mapInit$SR2$tag

   # Write longread data to mapFinal
   mapgroupLR <- "LRA"
   x$mapFinal$pileup[[mapgroupLR]] = x$mapInit$pileup
   x$mapFinal$tag[[mapgroupLR]] = x$mapInit$tag

   flog.info("Get latest consensus from last shortreadm apping ...", name = "info")
   cseq <- conseq(x$mapInit$SR2$pileup$consmat, "mapFinalA", "ambig", exclude_gaps = TRUE, threshold = x$getThreshold())
   x$mapFinal$seq$A <- x$mapInit$SR2$conseq

   flog.info("Polish and look for inconsistencies between shortreads and longreads ...", name = "info")
   polish(x)
   flog.info("Report consensus sequence and potential problematic variants", name = "info")
   report(x)

   return(x)
}

#' @export
# x <- dedk.report
#position = 5000
run_igv <- function(x, position, ...) {
  assertthat::assert_that(
    is(x, "DR2S"),
    is(x$consensus, "ConsList"),
    requireNamespace("futile.logger", quietly = TRUE)

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
    position <- ifelse(NROW(x$consensus$problematic_variants) == 0, 50, as.integer(x$consensus$problematic_variants[1, "pos"]))
  }

  for (hp in x$getHapTypes()){
    igv <- file.path(basedir, paste0("igv", hp, ".xml"), fsep = )
    ref   <- file.path(basedir, hp, basename(x$mapIter[[as.character(x$getIterations())]][[hp]]$seqpath))
    bamLR <- file.path(basedir, "final", basename(x$mapFinal$bamfile[[paste0("LR", hp)]]))
    bamSR <- file.path(basedir, "final", basename(x$mapFinal$bamfile[[paste0("SR", hp)]]))
    chr <- sub(">", "", readLines(ref, 1))
    locus <- paste0(chr, ":", position - 50, "-", position + 50)

    xml <- XML::xmlTree()
    xml$addTag("Global", attrs = c(genome = ref, locus = locus), close = FALSE)
    xml$addTag("Resources", close = FALSE)
    xml$addTag("Resource", attrs = c(path = bamLR))
    xml$addTag("Resource", attrs = c(path = bamSR))
    xml$closeTag()
    xml$closeTag()
    XML::saveXML(xml, file = igv)

    if (.Platform$OS.type != "windows") {
      system(paste0(shQuote(exe_path), " ", shQuote(igv), " &"))#, ...)
    } else if (.Platform$OS.type == "windows") {
      msg <- sprintf(
        "\nOpen an instances of the command prompt and paste in the command:\n%s\n",
        paste0(shQuote(exe_path), " ", shQuote(igv))
      )
      message(msg)
    }
  }
}
