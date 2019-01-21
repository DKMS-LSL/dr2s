mapInitSR <- function(self, opts = list(), includeDeletions = TRUE,
                      includeInsertions = TRUE, callInsertions = TRUE,
                      callInsertionThreshold = 0.15, clip = FALSE,
                      distributeGaps = FALSE, removeError = TRUE, topx = 0,
                      outdir, clean, ...) {

  ## get indenter
  indent <- list(...)$indent %||% indentation()

  ## Get flags
  forceMapping <- self$getForceMapping()
  microsat     <- self$getMicrosatellite()

  ## For all mappings
  maplabel <- "mapInit1"
  mapfun   <- self$getSrdMapFun()
  readtype <- self$getSrdType()
  readfile <- self$getShortreads()

  ## Primary mapping
  reffile  <- self$getRefPath()
  refname  <- self$getReference()
  pileup   <- mapReads(
    mapfun = mapfun, maplabel = maplabel, reffile = reffile, refname = refname,
    readfile = readfile, readtype = readtype, opts = opts, outdir = outdir,
    includeDeletions = includeDeletions, includeInsertions = includeInsertions,
    callInsertions = TRUE, callInsertionThreshold = callInsertionThreshold,
    clip = FALSE, distributeGaps = FALSE, removeError = TRUE, topx = 0,
    clean = clean, indent = indent, ...)#minMapq = 50)

  ## Check if the coverage is somewhat equally distributed
  if (max(rowSums(consmat(pileup, freq = FALSE))) /
      quantile(rowSums(consmat(pileup, freq = FALSE)), 0.75) > 5) {
    plotfile <- self$absPath("plot.MapInit.SR.problem.pdf")
    .checkCoverage(pileup, forceMapping, plotfile, maplabel, indent = incr(indent))
  }

  ## Construct primary initial consensus sequence
  consname <- refname %<<% ".consensus"
  conspath <- file.path(outdir, strip(maplabel %<<% "." %<<% consname %<<% ".fa", "_"))
  names(conspath) <- self$relPath(conspath)
  flog.info("%sConstruct consensus <%s>", indent(), names(conspath), name = "info")
  conseq <- .writeConseq(x = pileup, name = consname, type = "prob",
                         threshold = NULL, suppressAllGaps = TRUE,
                         replaceIndel = "N", conspath = conspath)

  if (microsat) {
    flog.info("%sRemap shortreads to extended reference", indent(), name = "info")
    pileup <- mapReads(
      mapfun = mapfun, maplabel = maplabel, reffile = conspath, refname = consname,
      readfile = readfile, readtype = readtype, opts = opts, outdir = outdir,
      includeDeletions = includeDeletions, includeInsertions = includeInsertions,
      callInsertions = TRUE, callInsertionThreshold = callInsertionThreshold,
      clip = FALSE, distributeGaps = FALSE, removeError = TRUE, topx = 0,
      clean = clean, indent = indent, ...)#minMapq = 50)

    # Construct secondary initial consensus sequence
    consname <- consname %<<% ".2"
    conspath <- file.path(outdir, strip(maplabel %<<% "." %<<% consname %<<% ".fa", "_"))
    names(conspath) <- self$relPath(conspath)
    flog.info("%sConstruct consensus <%s>", indent(), names(conspath), name = "info")
    conseq <- .writeConseq(x = pileup, name = consname, type = "prob",
                           threshold = NULL, suppressAllGaps = TRUE,
                           replaceIndel = "N", conspath = conspath)
  }

  SR1 = MapList_(
    ## mapdata
    readpath = self$relPath(readfile),
    refpath  = self$relPath(reffile),
    bampath  = self$relPath(bampath(pileup)),
    conspath = self$relPath(conspath),
    pileup   = pileup,
    stats    = list(coverage = .coverage(pileup)),
    ## metadata
    maplabel = maplabel,
    refname  = refname,
    mapper   = self$getSrdMapper(),
    opts     = opts)

  ## Second mapping to infer polymorphic positions
  ## from same reference as longreads
  reffile <- self$absPath(conspath(SR1))
  refname <- consname(SR1)
  maplabel <- "mapInit2"
  pileup  <- mapReads(
    mapfun = mapfun, maplabel = maplabel, reffile = reffile, refname = refname,
    readfile = readfile, readtype = readtype, opts = opts, outdir = outdir,
    includeDeletions = TRUE, includeInsertions = FALSE, callInsertions = FALSE,
    callInsertionThreshold = callInsertionThreshold, clip = FALSE,
    distributeGaps = FALSE, removeError = TRUE, topx = 0, clean = clean,
    indent = indent, ...)#minMapq = 50)

  SR2 = MapList_(
    ## mapdata
    refpath  = self$relPath(reffile),
    readpath = self$relPath(readfile),
    bampath  = self$relPath(bampath(pileup)),
    conspath = NULL,
    pileup   = pileup,
    stats    = list(coverage = .coverage(pileup)),
    ## metadata
    maplabel = maplabel,
    refname  = refname,
    mapper   = self$getSrdMapper(),
    opts     = opts)

  list(SR1 = SR1, SR2 = SR2)
}

