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
    stop("No DR2S analysis object found", call. = TRUE, immediate. = TRUE)
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

.finishCn1 <- function(x){
  flog.info("Set only allele to A", name = "info")
  hptype <- "A"
  x$setHapTypes(c(hptype))
  
  ## Write cons for mapIter
  mapIterDir <- .dirCreateIfNotExists(
    normalizePath(file.path(x$getOutdir(), "mapIter", hptype), mustWork = FALSE))
  maplabel <- "mapIter0"
  consname <- maplabel %<<% ".consensus." %<<% hptype
  conspath <- file.path(mapIterDir, consname %<<% ".fa")
  names(conspath) <- x$relPath(conspath)
  flog.info("Construct consensus <%s>", names(conspath), name = "info")
  conseq <- .writeConseq(x = x$mapInit$pileup$consmat, name = consname, type = "prob",
                         threshold = NULL, suppressAllGaps = TRUE,
                         replaceIndel = "", conspath = conspath)
  readfile <- dot(c(
    "hap" %<<% hptype, x$getLrdType(), x$getLrdMapper(),
    "n" %<<% "fastq", "gz"))
  readpath <- file.path(mapIterDir, readfile)
  file.copy(x$getLongreads(), readpath)
  ##
  x$mapIter$`0`[[hptype]] <- MapList_(
    ## mapdata
    readpath  = x$relPath(readpath),
    refpath   = refpath(x$mapInit),
    bampath   = bampath(x$mapInit),
    conspath  = x$relPath(conspath),
    pileup    = x$mapInit$pileup,
    stats     = list(),
    ## required metadata
    maplabel  = maplabel,
    refname   = refname(x$mapInit),
    mapper    = meta(x$mapInit, "mapper"),
    opts      = meta(x$mapInit, "opts")
  )
  x$setHomozygous(TRUE)
  # flog.info("Write mapInit data to mapIter", name = "info")
  # x$mapFinal$LR$A = x$mapInit
  # if (!is.null(x$mapInit$meta$SR2)) {
  #   # Write shortread data to mapFinal
  #   x$mapFinal$SR$A = x$mapInit$meta$SR2
  # }

  # flog.info("Get latest consensus from last mapping ...", name = "info")

  ## Construct consensus sequence
  # map <- if (is.null(x$mapFinal$SR)) {
  #   x$mapFinal$LR$A
  # } else {
  #   x$mapFinal$SR$A
  # }
  # maplabel <- "mapInit"
  # consname <- maplabel %<<% ".consensus." %<<% refname(map)
  # outdir   <- .dirCreateIfNotExists(x$absPath(maplabel))
  # conspath <- file.path(outdir, consname %<<% ".fa")
  # names(conspath) <- x$relPath(conspath)
  # flog.info("Construct consensus <%s>", names(conspath), name = "info")
  # conseq <- .writeConseq(x = map$pileup, name = consname, type = "ambig",
  #                        threshold = x$getThreshold(), suppressAllGaps = TRUE,
  #                        replaceIndel = "", conspath = conspath)
  # if (is.null(x$mapFinal$SR)) {
  #   x$mapFinal$LR$A$conspath  = x$relPath(conspath)
  # } else {
  #   x$mapFinal$SR$A$conspath  = x$relPath(conspath)
  # }
  # flog.info("Report consensus sequence and potential problematic variants",
  #           name = "info")
  # x <- report(x)
  return(invisible(x))
}

#' Create IGV xml config files for directly open the longread and shortread
#' mapping in an IGV session. It creates also ".sh" and ".bat" scripts for
#' easy opening the IGV instance.
#' @param x A \code{\link{DR2S_}} object.
#' @param position The position where IGV focuses on startup.
#' @param map Which mapping should be opened. One of "mapInit", "mapIter",
#' "mapFinal" or "remap".
#' @param open whether to open an IGV instance now.
#' @export
createIgvConfigs <- function(x, position, map = "mapFinal", open = TRUE) {
  assert_that(
    is(x, "DR2S")
  )
  map <- match.arg(map, c("mapInit",
                          "mapFinal",
                          "mapIter",
                          "remap"))

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
  haptypes <- if (map == "mapInit") {
    "Init"
  } else {
    x$getHapTypes()
  }
  for (hp in haptypes) {
    igv <- file.path(igvdir, "igv" %<<% hp %<<% map %<<% ".xml")
    if (map == "mapFinal") {
      ref   <- x$mapIter[[as.character(x$getIterations())]][[hp]]$seqpath
      bamLR <- file.path("mapFinal",
                         basename(x$mapFinal$bamfile[["LR" %<<% hp]]))
      if (!is.null(x$mapFinal$sreads[[hp]]))
        bamSR <- file.path("mapFinal",
                           basename(x$mapFinal$bamfile[["SR" %<<% hp]]))
    } else if (map == "remap") {
        ref   <- x$remap$LR[[hp]]$refpath
        bamLR <- x$remap$LR[[hp]]$bampath
        if (x$hasShortreads()) 
          bamSR <- x$remap$SR[[hp]]$bampath
    } else if (map == "mapInit") {
      if (x$hasShortreads()) {
        ref <- unname(x$mapInit$refpath)
        bamSR <- x$mapInit$meta$SR2$bampath
      } else {
        ref <- x$mapInit$refpath
      }
      bamLR <- x$mapInit$bampath

    } else if (map == "mapIter") {
      ref <- x$mapIter[[as.character(x$getIterations())]][[hp]]$conspath
      bamLR <- x$mapIter[[as.character(x$getIterations())]][[hp]]$bampath
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
  #mapper <- readDR2S("~/bioinf/DR2S/KIR/181128_3DL3CmitONT/outDevel/KIR3DL3/1397771")
  #bamfile <- mapper$absPath(mapper$mapFinal$SR$A$bampath)
  assert_that(
    file.exists(bamfile),
    endsWith(bamfile, ".bam"),
    file.exists(reference),
    endsWith(reference, ".fa"))
  ## Subsample the bam file
  resultList <- subSampleBam(bamfile = bamfile, ...)#fragmentReads = TRUE, sampleSize = 100)#, ...)
  ## Index the reference
  Rsamtools::indexFa(reference)
  resultList$referenceFile <- .cropPath(outdir, reference)
  resultList$original <- .cropPath(outdir, resultList$original)
  resultList$sampled <- .cropPath(outdir, resultList$sampled)
  structure(
    resultList,
   class = c("igvjs", "list")
  )
}
