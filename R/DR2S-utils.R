#' Read a cached DR2S object
#' @param path the path to the folder storing the DR2S object. Usually the
#' \code{outdir} of a previous run.
#'
#' @examples
#' \dontrun{
#' library(DR2S)
#' path <- "path"
#' dr2s <- readDR2S(path)
#' }
#' @export
readDR2S <- function(path) {
  outdir <- normalizePath(path, mustWork = TRUE)
  rdsPath <- dir(outdir, pattern = "^DR2S(YAML)?.*.rds", full.names = TRUE)
  if (length(rdsPath) == 0) {
    warning("No DR2S analysis object found", call. = TRUE, immediate. = TRUE)
    return(NULL)
  }
  rs <- lapply(rdsPath, readRDS)
  ## update outdir
  rs <- lapply(rs, function(x) {
    x$setConfig("outdir", outdir)
    x
  })
  .confLog(outdir = outdir, logName = "info")
  if (length(rs) == 1) {
    rs[[1]]
  } else rs
}


.finishCn1 <- function(x, plot = TRUE) {
  ## Initiate indenter
  indent <- indentation(1)
  indent2 <- incr(indent)

  flog.info("%sSet single allele to <A>", indent(), name = "info")

  x$setHapTypes("A")
  hp <- x$getHapTypes()
  maplabel <- "mapFinal"
  outdir   <- .dirCreateIfNotExists(self$absPath(maplabel))

  flog.info("%sWrite mapInit LR data into mapFinal LR", indent(), name = "info")
  x$mapFinal$LR[[hp]] = x$mapInit
  x$mapFinal$LR[[hp]]$meta$maplabel <- maplabel
  ##
  ## Construct consensus sequence
  ##
  refname  <- "LR" %<<% hp
  consname <- maplabel %<<% ".consensus." %<<% refname
  conspath <- file.path(outdir, consname %<<% ".fa")
  names(conspath) <- self$relPath(conspath)
  pileupLR <- x$mapFinal$LR[[hp]]$pileup
  flog.info("%sConstruct consensus <%s>", indent2(), names(conspath), name = "info")
  conseq <- .writeConseq(x = pileupLR, name = consname, type = "ambig",
                         threshold = 1/4, suppressAllGaps = TRUE,
                         replaceIndel = "", conspath = conspath)
  x$mapFinal$LR[[hp]]$conspath <- x$relPath(conspath)

  if (x$hasShortreads()) {
    flog.info("%sWrite mapInit SR data into mapFinal SR", indent(), name = "info")
    x$mapFinal$SR[[hp]] = meta(x$mapInit, "SR2")
    x$mapFinal$SR[[hp]]$meta$maplabel <- maplabel
    ##
    ## Construct consensus sequence
    ##
    refname  <- "SR" %<<% hp
    consname <- maplabel %<<% ".consensus." %<<% refname
    conspath <- file.path(outdir, consname %<<% ".fa")
    names(conspath) <- self$relPath(conspath)
    pileupSR <- x$mapFinal$SR[[hp]]$pileup
    flog.info("%sConstruct consensus <%s>", indent2(), names(conspath), name = "info")
    conseq <- .writeConseq(x = pileupSR, name = consname, type = "ambig",
                           threshold = 1/4, suppressAllGaps = TRUE,
                           replaceIndel = "", conspath = conspath)
    x$mapFinal$SR[[hp]]$conspath <- x$relPath(conspath)
  }

  if (plot) {
    flog.info("%sPlot MapFinal summary", indent(), name = "info")
    ## Coverage and base frequency
    readtypes <- if (self$hasShortreads()) c("LR", "SR") else "LR"
    ## readtype = "LR"
    plotlist <- foreach(readtype = readtypes) %do% {
      suppressWarnings(self$plotMapFinalSummary(iteration = "final",
                                                readtype = readtype, thin = 0.25, width = 20))
    }
    plotRows <- 1
    hptypes  <- hp
    p <- cowplot::plot_grid(plotlist = plotlist, nrow = plotRows, labels = readtypes)
    cowplot::save_plot(p, filename = self$absPath("plot.MapFinal.pdf"),
                       base_width = 12*length(hptypes),
                       title = paste(self$getLocus(), self$getSampleId(),
                                     sep = "." ),
                       base_height = 3*length(readtypes))
    cowplot::save_plot(p, filename = self$absPath(".plots/plot.MapFinal.svg"),
                       base_width = 12*length(hptypes),
                       base_height = 3*length(readtypes))
  }

  flog.info(
    "Polish and check for inconsistencies between shortreads and longreads", name = "info")
  polish(x)

  flog.info(
    "Report consensus sequence and potential problematic variants", name = "info")
  report(x)

  return(invisible(x))
}

#' Create IGV xml config files for directly open the longread and shortread
#' mapping in an IGV session. It creates also ".sh" and ".bat" scripts for
#' easy opening the IGV instance.
#' @param x A \code{\link{DR2S_}} object.
#' @param position The position where IGV focuses on startup.
#' @param map Which mapping should be opened. One of "mapInit", "mapIter",
#' "mapFinal" or "refine".
#' @param open whether to open an IGV instance now.
#' @export
createIgvConfigs <- function(x, position, map = "mapFinal", open = TRUE) {
  assert_that(is(x, "DR2S"))
  map <- match.arg(map, c("mapInit", "mapIter", "mapFinal", "refine"))
  basedir <- normalizePath(x$getOutdir(), mustWork = TRUE)
  igvdir <- file.path(basedir, ".pplib")
  .dirCreateIfNotExists(igvdir)
  .dirCreateIfNotExists(file.path(basedir, "win"))
  if (.Platform$OS.type == "windows") {
    fsep <- "\\"
  } else {
    fsep <- "/"
  }
  if (missing(position)) {
    position <- ifelse(NROW(x$consensus$variants) == 0,
                       50,
                       as.integer(x$consensus$variants[1, "pos"]))
  }
  igvConfigs <- list()
  # hp = "Init"
  # hp = "A"
  hptypes <- if (map == "mapInit") "Init" else x$getHapTypes()
  for (hp in hptypes) {
    igv <- file.path(igvdir, "igv" %<<% hp %<<% map %<<% ".xml")
    ##
    if (map == "mapInit") {
      if (x$hasShortreads()) {
        ref <- refpath(meta(x$mapInit, "SR2"))
        bamSR <- bampath(meta(x$mapInit, "SR2"))
      } else {
        ref <- refpath(x$mapInit)
      }
      bamLR <- bampath(x$mapInit)
    }
    else if (map == "mapIter") {
      ref <- refpath(x$mapIter[[max(names(x$mapIter))]][[hp]])
      bamLR <- bampath(x$mapIter[[max(names(x$mapIter))]][[hp]])
    }
    else if (map == "mapFinal") {
      ref   <- refpath(x$mapIter[[max(names(x$mapIter))]][[hp]])
      bamLR <- bampath(x$mapFinal$LR[[hp]])
      if (x$hasShortreads()) {
        bamSR <- bampath(x$mapFinal$SR[[hp]])
      }
    }
    else if (map == "refine") {
      if (!is.null(x$consensus$refine$ref[[hp]])) {
        ref   <- x$consensus$refine$ref[[hp]]
        bamLR <- x$consensus$refine$bamfile[["LR" %<<% hp]]
        if (x$hasShortreads())
          bamSR <- x$consensus$refine$bamfile[["SR" %<<% hp]]
      } else {
        ref   <- refpath(x$mapIter[[max(names(x$mapIter))]][[hp]])
        bamLR <- bampath(x$mapFinal$LR[[hp]])
        if (x$hasShortreads()) {
          bamSR <- bampath(x$mapFinal$SR[[hp]])
        }
      }
    }
    ##
    chr <- strsplit1(sub(">", "", readLines(file.path(basedir, ref), 1)), "\\s+")[1]
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
    xml <- XML::xmlTree()
    suppressWarnings(xml$addTag("Global",
                                c(genome = file.path("..",
                                                     gsub("/", "\\\\", ref),
                                                     fsep = "\\"),
                                  locus = locus),
                                close = FALSE))
    xml$addTag("Resources", close = FALSE)
    xml$addTag("Resource", c(path = file.path("..",
                                              gsub("/", "\\\\", bamLR),
                                              fsep = "\\")))
    if (exists("bamSR"))
      xml$addTag("Resource", c(path = file.path("..",
                                                gsub("/", "\\\\", bamSR),
                                                fsep = "\\")))
    xml$closeTag()
    xml$closeTag()
    XML::saveXML(xml, file = gsub(".xml", ".win.xml", igv))
    igvConfigs[[hp]] <- x$relPath(igv)
  }
  igvCommand <- file.path(basedir, "runIGV_" %<<% map %<<% ".sh")
  cmds <- "igv " %<<% igvConfigs %<<% " &"
  write(cmds, igvCommand)
  Sys.chmod(igvCommand, mode = "775")

  ## for windows
  igvCommandWin <- file.path(basedir, "win", "runIGV_" %<<% map %<<% ".bat")
  cmdsWin <- "START P:\\IGV_2.4.10\\jre1.8.0_131\\bin\\javaw -Xmx2G " %<<%
    "-jar P:\\IGV_2.4.10\\lib\\igv.jar  " %<<% "%~dp0..\\" %<<%
    gsub("/", "\\\\", gsub(".xml", ".win.xml", igvConfigs))
  write(cmdsWin, igvCommandWin)

  if (open) {
    if (.Platform$OS.type != "windows") {
      cwd <- getwd()
      setwd(basedir)
      system(igvCommand)
      setwd(cwd)
    } else {
      message("\nOpen IGV manually")
    }
  }
}


#' Create subsampled bam files and fasta index files for igv.js instances.
#'
#' @param reference The reference fasta file used for mapping
#' @param bamfile The bamfile which should be viewed in IGV.js
#' @param ... Additional parameters passed to \code{\link{subSampleBam}}.
#' @export
createIgvJsFiles <- function(reference, bamfile, outdir, ...) {
  assert_that(
    file.exists(bamfile),
    endsWith(bamfile, ".bam"),
    file.exists(reference),
    endsWith(reference, ".fa"))
  ## Subsample the bam file
  resultList <- subSampleBam(bamfile = bamfile, ...)#fragmentReads = TRUE, sampleSize = 100)#, ...)
  ## Index the reference
  indexFa(reference)
  resultList$referenceFile <- .cropPath(outdir, reference)
  resultList$original <- .cropPath(outdir, resultList$original)
  resultList$sampled <- .cropPath(outdir, resultList$sampled)
  structure(
    resultList,
   class = c("igvjs", "list")
  )
}
