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


.finishCn1 <- function(x) {
  flog.info("Set only allele to A", name = "info")
  x$setHapTypes("A")

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

  flog.info("Get consensus from final mapping ...", name = "info")

  if (!is.null(x$mapInit$SR1)) {
    cseq <- conseq(consmat(x$mapInit$SR2$pileup), "mapFinalA", "ambig",
                   excludeGaps = TRUE, threshold = x$getThreshold())
  } else {
    cseq <- conseq(consmat(x$mapInit$pileup), "mapFinalA", "ambig",
                   excludeGaps = TRUE, threshold = x$getThreshold())
  }
  x$mapFinal$seq$A <- cseq

  flog.info(
    "Polish and look for inconsistencies between shortreads and longreads ...",
    name = "info")
  polish(x)

  flog.info("Report consensus sequence and potential problematic variants",
            name = "info")
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
  assert_that(
    is(x, "DR2S")
  )
  map <- match.arg(map, c("mapInit",
                          "mapFinal",
                          "mapIter",
                          "refine"))

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
  # hp = "A"
  haptypes <- if(map == "mapInit") {
    "Init"
  } else {
    x$getHapTypes()
  }
  for (hp in haptypes){
    igv <- file.path(igvdir, "igv" %<<% hp %<<% map %<<% ".xml")
    if (map == "mapFinal") {
      ref   <- x$mapIter[[as.character(x$getIterations())]][[hp]]$seqpath
      bamLR <- file.path("mapFinal",
                         basename(x$mapFinal$bamfile[["LR" %<<% hp]]))
      if (!is.null(x$mapFinal$sreads[[hp]]))
        bamSR <- file.path("mapFinal",
                           basename(x$mapFinal$bamfile[["SR" %<<% hp]]))
    } else if (map == "refine") {
      if (!is.null(x$consensus$refine$ref[[hp]])) {
        ref   <- x$consensus$refine$ref[[hp]]
        bamLR <- x$consensus$refine$bamfile[["LR" %<<% hp]]
        if (!is.null(x$mapFinal$sreads[[hp]]))
          bamSR <- x$consensus$refine$bamfile[["SR" %<<% hp]]
      } else {
        ref   <- x$mapIter[[as.character(x$getIterations())]][[hp]]$seqpath
        bamLR <- file.path("mapFinal",
                           basename(x$mapFinal$bamfile[["LR" %<<% hp]]))
        if (!is.null(x$mapFinal$sreads[[hp]]))
          bamSR <- file.path("mapFinal",
                             basename(x$mapFinal$bamfile[["SR" %<<% hp]]))
      }
    } else if (map == "mapInit") {
      if (x$hasShortreads()) {
        ref <- unname(x$mapInit$SR2$seqpath)
        bamSR <- x$mapInit$SR2$bamfile
      } else {
        ref <- x$getRefSeq()
      }
      bamLR <- x$mapInit$bamfile

    } else if (map == "mapIter") {
      ref <- x$mapIter[[as.character(x$getIterations() - 1)]][[hp]]$seqpath
      bamLR <- x$mapIter[[as.character(x$getIterations())]][[hp]]$bamfile
    }
    if (!is(ref, "DNAStringSet")) {
      chr <- strsplit1(sub(">", "",
                           readLines(file.path(basedir, ref), 1)), "\\s+")[1]
    } else {
      chr <- names(ref)
    }
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
