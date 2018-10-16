# Method: MapInit ####
#' @export
mapInit.DR2S <- function(x,
                         opts = list(),
                         optsname = "",
                         partSR = TRUE,
                         threshold = NULL,
                         minBaseQuality = 3,
                         minMapq = 50,
                         maxDepth = 1e4,
                         minNucleotideDepth = 3,
                         includeDeletions = TRUE,
                         includeInsertions = TRUE,
                         microsatellite = FALSE,
                         force = FALSE,
                         fullname = FALSE,
                         filterScores = TRUE,
                         forceMapping = FALSE,
                         topx = 0,
                         createIgv = TRUE,
                         plot = TRUE) {
  x$runMapInit(opts = opts,
               optsname = optsname,
               partSR = partSR,
               threshold = threshold,
               minBaseQuality = minBaseQuality,
               minMapq = minMapq,
               maxDepth = maxDepth,
               minNucleotideDepth = minNucleotideDepth,
               includeDeletions = includeDeletions,
               includeInsertions = includeInsertions,
               microsatellite = microsatellite,
               force = force,
               fullname = fullname,
               filterScores = filterScores,
               forceMapping = forceMapping,
               topx = topx,
               createIgv = createIgv,
               plot = plot)
  invisible(x)
}

DR2S_$set("public", "runMapInit", function(opts = list(),
                                           optsname = "",
                                           partSR = TRUE,
                                           threshold = NULL,
                                           minBaseQuality = 3,
                                           minMapq = 50,
                                           maxDepth = 1e4,
                                           minNucleotideDepth = 3,
                                           includeDeletions = TRUE,
                                           includeInsertions = TRUE,
                                           microsatellite = FALSE,
                                           force = FALSE,
                                           fullname = FALSE,
                                           filterScores = TRUE,
                                           forceMapping = FALSE,
                                           topx = 0,
                                           createIgv = TRUE,
                                           plot = TRUE) {

  flog.info("Step 0: mapInit ...", name = "info")

  # # debug
  # opts = list()
  # optsname = ""
  # partSR = FALSE
  # threshold = 0.20
  # minBaseQuality = 3
  # minMapq = 50
  # maxDepth = 1e4
  # minNucleotideDepth = 3
  # includeDeletions = TRUE
  # includeInsertions = TRUE
  # microsatellite = TRUE
  # force = FALSE
  # fullname = FALSE
  # filterScores = FALSE
  # forceMapping = TRUE
  # plot = TRUE
  # topx <- 0
  # createIgv = TRUE
  # library(ggplot2)
  # library(S4Vectors)
  # library(Rsamtools)
  # library(foreach)
  # library(futile.logger)
  # library(cowplot)
  # self <-dr2s
  # self <- mapper

  ## Overide default arguments
  args <- self$getOpts("mapInit")
  if (!is.null(args)) {
    env  <- environment()
    list2env(args, envir = env)
  }
  if (is.null(threshold)) {
    threshold <- self$getThreshold()
  }

  ## Get all options and prepare mapping
  microsatellite  <- self$getMicrosatellite()
  partSR          <- self$getPartSR()
  forceMapping    <- self$getForceMapping()
  filterScores    <- self$getFilterScores()
  outdir          <- .dirCreateIfNotExists(self$absPath("mapInit"))
  .dirCreateIfNotExists(path = file.path(self$absPath(".plots")))
  igv <- list()
  clean <- TRUE

  if (partSR) {
    mapLabel <- "mapInit1"
    reffile  <- self$getRefPath()
    allele   <- self$getReference()
    readtype <- self$getSrdType()

    maptag  <- paste(mapLabel, paste(litArrows(c(allele, readtype,
                                                 self$getSrMapper(),
                                                 optstring(opts, optsname))),
                                     collapse = " "))

    flog.info(" Map shortreads to provided reference", name = "info")

    pileup <- mapReads(maptag = maptag, reffile = reffile,
                       readfile = self$getShortreads(), threshold = threshold,
                       allele = allele, readtype = readtype, opts = opts,
                       outdir = outdir, minMapq = minMapq, optsname = optsname,
                       minBaseQuality = minBaseQuality, maxDepth = maxDepth,
                       minNucleotideDepth = minNucleotideDepth, force = force,
                       includeDeletions = includeDeletions, clean = clean,
                       includeInsertions = includeInsertions,
                       callInsertions = TRUE, mapFun = self$getSrMapFun())
    # ## TODO: maybe bum this?
    # if (filterScores) {
    #   flog.info(" Filter reads with low alignment scores", name = "info")
    #   ## Run bam - sort - index pipeline
    #   bamfile <- .bamSortIndex(samfile, self$getRefPath(),
    #                             minMapq, force = force, clean = TRUE)
    #   ## Filter Reads
    #   bam <- scanBam(bamfile,
    #                  param = ScanBamParam(tag = "AS",
    #                  what = c("qname", "pos", "cigar")))[[1]]
    #   readfilter <- .filterReads(bam = bam, preserveRefEnds = TRUE)
    #   .fileDeleteIfExists(bamfile)
    #
    #   flog.info(" Write new shortread fastqs to file", name = "info")
    #   fqs <- self$getShortreads()
    #   fqdir <- .dirCreateIfNotExists(file.path(outdir,self$getSrdType()))
    #   # write fastq's
    #   readfile <- c()
    #   readfile <- foreach(fq = fqs, .combine = c) %do% {
    #     srFastqHap = file.path(fqdir, basename(fq))
    #     .writePartFq(fq = fq, srFastqHap = srFastqHap,
    #                   dontUseReads = readfilter)
    #     srFastqHap
    #   }
    #   # set new shortread directory
    #   self$setConfig("filteredShortreads", self$relPath(fqdir))
    #
    #   flog.info(" Map filtered shortreads to provided reference", name = "info")
    #   ## Rerun mapper
    #   flog.info("  Mapping ...", name = "info")
    #   samfile <- mapFun(
    #     reffile  = self$getRefPath(),
    #     readfile = readfile,
    #     allele   = self$getReference(),
    #     readtype = self$getSrdType(),
    #     opts     = opts,
    #     refname  = "",
    #     optsname = optsname,
    #     force    = force,
    #     outdir   = outdir
    #   )
    # }

    ## Check if the coverage is somewhat equally distributed
    if (max(rowSums(pileup$consmat)) /
        quantile(rowSums(pileup$consmat), 0.75) > 5) {
      plotFile <- self$absPath("plot.MapInit.SR.problem.pdf")
      .checkCoverage(pileup, forceMapping, plotFile, maptag)
    }

    # calc initial consensus
    flog.info(" Construct initial consensus from shortreads", name = "info")
    ## Get conseq
    conseqName <- "Init.consensus." %<<%
      sub(".bam", "", basename(pileup$bamfile))
    conseqPath  <- file.path(outdir, conseqName %<<% ".fa")
    conseq <- .getWriteConseq(pileup = pileup, name = "mapInit",
                              type = "prob",  threshold = threshold,
                              forceExcludeGaps = TRUE, conseqPath = conseqPath)

    if (microsatellite) {
      mapLabel  <- "mapInit1.2"
      reffile <- conseqPath
      allele <- conseqName
      readtype <- self$getSrdType()
      maptag  <- paste(mapLabel, paste0(litArrows(c(conseqName, readtype,
                                                    self$getSrMapper(),
                                                    optstring(opts, optsname))),
                                        collapse = " "))

      flog.info(" Refine microsatellites or repeats by extending the reference",
                name = "info")
      flog.info(" Remap shortreads to initial consensus from shortreads",
                name = "info")

      pileup <- mapReads(maptag = maptag, reffile = reffile,
                         readfile = self$getShortreads(), threshold = threshold,
                         allele = allele, readtype = readtype, opts = opts,
                         outdir = outdir, minMapq = minMapq, force = force,
                         minBaseQuality = minBaseQuality, optsname = optsname,
                         minNucleotideDepth = minNucleotideDepth, clean = clean,
                         mapFun = self$getSrMapFun(), maxDepth = maxDepth,
                         includeDeletions = includeDeletions,
                         includeInsertions = includeInsertions,
                         callInsertions = TRUE)


      # Infer initial consensus
      flog.info(" Construct second consensus from shortreads " %<<%
                  "with refined repeats", name = "info")
      conseqName <- "Init.consensus.2" %<<% sub(".bam", "",
                                                basename(pileup$bamfile))
      conseqPath  <- file.path(outdir, conseqName %<<% ".fa")
      conseq <- .getWriteConseq(pileup, name = "mapInit1.2",
                                type = "prob",  threshold = threshold,
                                forceExcludeGaps = TRUE,
                                conseqPath = conseqPath)
    }
    mapInitSR1 = structure(
      list(
        reads   = self$relPath(self$getShortreads()),
        bamfile = self$relPath(pileup$bamfile),
        pileup  = pileup,
        tag     = maptag,
        conseq  = conseq,
        seqpath = self$relPath(conseqPath),
        ref     = conseqName
      ),
      class  = c("mapInit", "list")
    )

    ## Second mapping to infer polymorphic positions
    ## from same reference as longreads
    reffile  <- self$absPath(mapInitSR1$seqpath)
    allele   <- mapInitSR1$ref
    readtype <- self$getSrdType()
    mapLabel <- "mapInit2"
    maptag   <- paste(mapLabel, paste0(litArrows(c(allele, readtype,
                                                   self$getSrMapper(),
                                                   optstring(opts, optsname))),
                                       collapse = " "))

    flog.info(" Remap shortreads to consensus for SNP calling", name = "info")
    pileup <- mapReads(maptag = maptag, reffile = reffile,
                       readfile = self$getShortreads(), threshold = threshold,
                       allele = allele, readtype = readtype, opts = opts,
                       outdir = outdir, minMapq = minMapq, optsname = optsname,
                       minBaseQuality = minBaseQuality, maxDepth = maxDepth,
                       minNucleotideDepth = minNucleotideDepth, force = force,
                       includeDeletions = TRUE, clean = clean,
                       includeInsertions = FALSE,  callInsertions = FALSE,
                       mapFun = self$getSrMapFun())

    ### TODO wrap this command up
    if (createIgv)
      igv[["SR"]] <- createIgvJsFiles(reffile, pileup$bamfile,
                                      self$getOutdir(),
                                      sampleSize = 100)

    mapInitSR2 = structure(
      list(
        reads   = self$relPath(self$getShortreads()),
        bamfile = self$relPath(pileup$bamfile),
        pileup  = pileup,
        tag     = maptag,
        conseq  = conseq,
        seqpath = self$relPath(conseqPath),
        ref     = conseqName
      ),
      class  = c("mapInit", "list")
    )
  }

  ## Map longreads
  if (exists("mapInitSR1")) {
    flog.info(" Map longreads to consensus for clustering",
              name = "info")
    reffile <- self$absPath(mapInitSR1$seqpath)
    refseq  <- mapInitSR1$conseq
    allele  <- mapInitSR1$ref
    conseqName <- allele
    maprepeat  <- FALSE
  } else {
    flog.info(" Map longreads to provided reference for clustering",
              name = "info")
    reffile <- self$getRefPath()
    refseq  <- self$getRefSeq()
    allele  <- self$getReference()
    conseqName <- allele
    maprepeat  <- TRUE
  }
  readfile <- self$getLongreads()
  readtype <- self$getLrdType()
  mapFun   <- self$getLrMapFun()
  mapLabel <- "mapInit"
  maptag   <- paste(mapLabel, paste0(litArrows(c(conseqName, readtype,
                                                 self$getLrMapper(),
                                                 optstring(opts, optsname))),
                                     collapse = " "))

  if (maprepeat) {
    pileup <- mapReads(maptag = maptag, reffile = reffile, readfile = readfile,
                       allele = allele, readtype = readtype, opts = opts,
                       optsname = optsname, force = force, maxDepth = maxDepth,
                       outdir = outdir, minMapq = minMapq, threshold = threshold,
                       minBaseQuality = minBaseQuality, clean = clean,
                       minNucleotideDepth = minNucleotideDepth,
                       includeDeletions = TRUE, includeInsertions = TRUE,
                       callInsertions = TRUE, mapFun = mapFun,
                       distributeGaps = TRUE, refseq = refseq, topx = topx)
    if (!is.null(pileup$reads)) {
      file <- paste(
        self$getSampleId(), self$getLrdType(), paste0("n", topx),
        "fastq", "gz", sep = ".")
      fqout <- .fileDeleteIfExists(file.path(outdir, file))
      fq  <- .extractFastq(pileup$bamfile, pileup$reads)
      ShortRead::writeFastq(fq, fqout, compress = TRUE)
      readfile <- fqout
    }

    conseqName <- "Init.LRconsensus." %<<%
      sub(".bam", "", basename(pileup$bamfile))
    maptag   <- paste(mapLabel, paste0(litArrows(c(conseqName, readtype,
                                                   self$getLrMapper(),
                                                   optstring(opts, optsname))),
                                       collapse = " "))
    self$setConfig("refPath", file.path(basename(outdir), conseqName %<<% ".fa"))
    reffile <- self$getRefPath()
    refseq  <- .getWriteConseq(pileup = pileup, name = "mapInitLR",
                               type = "prob",  threshold = threshold,
                               forceExcludeGaps = TRUE, conseqPath = reffile)
    allele  <- conseqName
  }

  pileup <- mapReads(maptag = maptag, reffile = reffile, readfile = readfile,
                     allele = allele, readtype = readtype, opts = opts,
                     optsname = optsname, force = force, maxDepth = maxDepth,
                     outdir = outdir, minMapq = minMapq, threshold = threshold,
                     minBaseQuality = minBaseQuality, clean = clean,
                     minNucleotideDepth = minNucleotideDepth,
                     includeDeletions = TRUE, includeInsertions = FALSE,
                     callInsertions = FALSE, mapFun = mapFun,
                     distributeGaps = TRUE, refseq = refseq)

  if (createIgv)
    igv[["LR"]] <- createIgvJsFiles(reffile, pileup$bamfile, outdir,
                                    sampleSize = 100, fragmentReads = TRUE)

  self$mapInit = structure(
    list(
      reads   = self$relPath(readfile),
      bamfile = self$relPath(pileup$bamfile),
      pileup  = pileup,
      tag     = maptag,
      SR1     = NULL,
      SR2     = NULL,
      igv     = igv
    ),
    class  = c("mapInit", "list")
  )

  if (partSR) {
    self$mapInit$SR1 <- mapInitSR1
    self$mapInit$SR2 <- mapInitSR2
  }
  createIgvConfigs(x = self, map = "mapInit", open = "FALSE")

  if (plot) {
    flog.info(" Plot MapInit summary ", name = "info")
    ## Coverage and frequency of minor alleles
    p <- self$plotmapInitSummary(
      thin = 0.25,
      width = 2
    )
    plotRows <- ifelse(partSR, 2, 1)
    cowplot::save_plot(self$absPath("plot.MapInit.pdf"),
                       plot = p, ncol = 1, nrow = plotRows,
                       base_aspect_ratio = as.numeric(paste(5, plotRows,sep = ".")),
                       title = paste(self$getLocus(), self$getSampleId(), sep = "." ))
    cowplot::save_plot(self$absPath(".plots/plot.MapInit.svg"),
                       plot = p, ncol = 1, nrow = plotRows,
                       base_aspect_ratio = as.numeric(paste(5, plotRows, sep = ".")))
  }
  return(invisible(self))
})

#' @export
print.mapInit <- function(x, ...) {
  msg <- sprintf("An object of class '%s'\n", class(x)[1])
  msg <- sprintf(
    "%s [Tag]     %s\n [Reads]   %s\n [Bamfile] %s\n [Pileup]\n",
    msg,
    x$tag,
    basename(x$reads),
    basename(x$bamfile)
  )
  cat(msg)
}

## Method: partitionLongReads ####
#' @export
partitionLongReads.DR2S <- function(x,
                                    threshold         = NULL,
                                    skipGapFreq       = 2/3,
                                    distAlleles       = NULL,
                                    noGapPartitioning = FALSE,
                                    selectAllelesBy   = "count",
                                    minClusterSize    = 15,
                                    plot              = TRUE,
                                    ...) {
  x$runPartitionLongReads(threshold         = threshold,
                          skipGapFreq       = skipGapFreq,
                          noGapPartitioning = noGapPartitioning,
                          selectAllelesBy   = selectAllelesBy,
                          minClusterSize    = minClusterSize,
                          distAlleles       = distAlleles,
                          plot              = plot,
                          ...)
  x$runSplitLongReadsByHaplotype(plot = plot)
  x$runExtractLongReads()
  x$runGetPartitionedConsensus()
  invisible(x)
}

DR2S_$set("public",
          "runPartitionLongReads",
          function(threshold = NULL,
                   skipGapFreq = 2/3,
                   distAlleles = NULL,
                   noGapPartitioning = FALSE,
                   selectAllelesBy = "count",
                   minClusterSize = 15,
                   plot = TRUE,
                   ...) {
            # debug
            # threshold = NULL
            # skipGapFreq = 2/3
            # distAlleles = NULL
            # noGapPartitioning = TRUE
            # selectAllelesBy = "distance"
            # minClusterSize = 15
            # plot = TRUE
            # self <- dr2s
            # library(futile.logger)
            # library(assertthat)

            flog.info("Step 1: PartitionLongReads ...", name = "info")
            flog.info(" Partition longreads into haplotypes", name = "info")

            ## Overide default arguments
            args <- self$getOpts("partitionLongReads")
            if (!is.null(args)) {
              env  <- environment()
              list2env(args, envir = env)
            }
            if (is.null(threshold)) {
              threshold <- self$getThreshold()
            }
            if (is.null(distAlleles)) {
              distAlleles <- self$getDistAlleles()
            }
            assert_that(
              self$hasPileup(),
              is.double(skipGapFreq),
              is.double(threshold),
              is.count(distAlleles),
              is.logical(plot)
            )


            ## Get the reference sequence
            if (!is.null(self$mapInit$SR1)) {
              useSR <- TRUE
              refseq <- self$mapInit$SR1$conseq
              flog.info(" Construct SNP matrix from shortreads", name = "info")
            } else {
              useSR <- FALSE
              refseq <- self$getRefSeq()
              flog.info(" Construct SNP matrix from longreads", name = "info")
            }
            ppos <- self$polymorphicPositions(useSR = useSR)

            ## Spurious gaps, especially in longreads can hinder a correct clustering
            ## Remove gap positions for clustering
            if (noGapPartitioning) {
              flog.info(" Use only non-gap positions for clustering",
                        name = "info")
              ppos <- ppos %>%
                dplyr::filter(a1 != "-" & a2 != "-")
            }

            ## Check if already finished because it is a homozygous sample
            if (NROW(ppos) == 0) {
              flog.warn(" No polymorphic positions for clustering! Only single allele?",
                        name = "info")
              flog.info(" Entering polish and report pipeline", name = "info")
              return(invisible(.finishCn1(self)))
            }

            mat <- if (tryCatch(
              !is(self$partition, "PartList"),
              error = function(e)
                TRUE
            ) ||
            !(all(ppos$position %in% colnames(self$partition$mat)) &&
              all(colnames(self$partition$mat) %in% ppos$position))) {
              SNPmatrix(bamfile = self$absPath(self$mapInit$bamfile), refseq = refseq,
                        polymorphicPositions = ppos)
            } else {
              self$partition$mat
            }

            flog.info(" Partition %s longreads over %s SNPs", NROW(mat), NCOL(mat),
                      name = "info")
            prt <- partitionReads(x = mat,
                                  skipGapFreq = skipGapFreq,
                                  deepSplit = 1,
                                  threshold = threshold,
                                  distAlleles = distAlleles,
                                  sortBy = selectAllelesBy,
                                  minClusterSize = minClusterSize)
            ## Set sample haplotypes
            self$setHapTypes(levels(as.factor(PRT(prt))))

            # Check if we have only one cluster and finish the pipeline if so
            if (length(self$getHapTypes()) == 1) {
              flog.warn(" Only one allele left!")
              flog.info(" Entering polish and report pipeline", name = "info")
              return(invisible(.finishCn1(self)))
            }

            self$partition = structure(list(
              mat = mat,
              prt = prt,
              hpl = NULL,
              lmt = NULL
            ),
            class = c("PartList", "list"))
            return(invisible(self))
          })

#' @export
print.PartList <- function(x, ...) {
  msg <- sprintf("An object of class '%s'.\n", class(x)[1])
  msg <- sprintf("%s [Matrix]    %s reads; %s polymorphic positions\n",
                 msg, NROW(x$mat), NCOL(x$mat))
  cat(msg)
  cat(" [Partition] ")
  print(x$prt, ncols = 4, nrows = 6)
  cat(" [Haplotypes] ")
  print(x$hpl)
}

DR2S_$set("public", "runSplitLongReadsByHaplotype", function(plot = TRUE) {

  flog.info(" Split partitioned longreads by score", name = "info")
  ## Check if reporting is already finished and exit safely
  if (.checkReportStatus(self)) return(invisible(self))
  assert_that(self$hasPartition())

  ## Overide default arguments
  args <- self$getOpts("partitionLongReads")
  if (!is.null(args)) {
    env  <- environment()
    list2env(args, envir = env)
  }

  prt  <- partition(self$getPartition())
  haplotypes <- levels(prt$haplotype)
  tag  <- self$getMapTag("init", "LR")

  # Set all limits to NULL
  self$setLimits(rep(NA, length(haplotypes)))
  prts <-  lapply(haplotypes, function(x) prt[prt$haplotype == x,])
  names(prts) <- haplotypes
  scores <- lapply(prts, function(x) x$mcoef)
  lmts <- .optimalPartitionLimits(scores)

  self$setLimits(setNames(as.list(lmts$limits$c), haplotypes))
  self$partition$lmt <- lmts$plt

  # Get only reads within the limit
  reads <- setNames(lapply(names(self$getLimits()), function(x) {
    dplyr::filter(prt, haplotype == x, mcoef >= self$getLimits()[x])
  }), names(self$getLimits()))
  for (hp in haplotypes) {
    flog.info("  %s: Using %s longreads with score > %.2f",
              hp, nrow(reads[[hp]]), self$getLimits()[hp], name = "info")
  }

  # results structure
  self$partition$hpl <- structure(
    setNames(lapply(reads, function(x) {
      structure(
        x$read,
        q = x$mcoef,
        freq = NROW(x)/NROW(dplyr::bind_rows(reads)),
        limit = self$getLimits()[[x$haplotype[1]]]
      )}), haplotypes),
    class = c("HapList", "list"))

  if (plot) {
    ## Plot the consensus sequences from clustering
    .browseSeqs(SQS(self$getPartition()),
                file = self$absPath("partition.fa.html"),
                openURL = FALSE)

    p <- self$plotPartitionSummary(label = tag,
                                   limits = unlist(self$getLimits()))

    cowplot::save_plot(self$absPath("plot.Partition.pdf"), plot = p,
                       title = paste(self$getLocus(), self$getSampleId(), sep = "."),
                       base_height = 12, base_width = 10)
    cowplot::save_plot(self$absPath(".plots/plot.Partition.svg"), plot = p,
                       base_aspect_ratio = 1.2)

    outf  <- self$absPath("plot.Sequence")
    ppos <- SNP(self$getPartition())
    names(ppos) <- seq_along(ppos)
    pwm <- lapply(PWM(self$getPartition()), function(pwm) {
      pwm[pwm < 0.1] <- 0
      pwm
    })
    p <- self$plotSeqLogo(ppos, pwm)
    cowplot::save_plot(filename = self$absPath("plot.Sequence.pdf"),
                       plot        = p,
                       base_width  = 0.4*length(ppos) + 1.4,
                       base_height = 2.5*length(pwm),
                       title       = paste(self$getLocus(), self$getSampleId(),
                                           sep = "." ),
                       units       = "cm",
                       limitsize   = FALSE)
    cowplot::save_plot(filename = self$absPath(".plots/plot.Sequence.svg"),
                       plot        = p,
                       base_width  = 0.4*length(ppos) + 1.4,
                       base_height = 2.5*length(pwm),
                       units       = "cm",
                       limitsize   = FALSE)

  }
  return(invisible(self))
})

#' @export
print.HapList <- function(x, ...) {
  msg <- sprintf("An object of class '%s'.\n", class(x)[1])
  for (haplotype in names(x)) {
    msg <- msg %<<%
      sprintf("%s: n %s; frequency %s; limit %s\n",
              haplotype,
              length(x[haplotype]),
              round(attr(x[[haplotype]], "freq"), 3),
              attr(x[[haplotype]], "limit"))
  }
  cat(msg)
}

DR2S_$set("public", "runExtractLongReads", function() {

  flog.info(" Extract haplotyped longreads", name = "info")

  ## Check if reporting is already finished and exit safely
  if (.checkReportStatus(self)) return(invisible(self))

  assert_that(self$hasHapList())

  ## Overide default arguments
  args <- self$getOpts("partitionLongReads")
  if (!is.null(args)) {
    env  <- environment()
    list2env(args, envir = env)
  }

  ## do this for each haptype
  hptypes <- self$getHapTypes()
  for (hptype in hptypes) {
    dir <- .dirCreateIfNotExists(
      normalizePath(file.path(self$getOutdir(), "mapIter", (hptype)),
                    mustWork = FALSE)
    )
    qnames <- self$getHapList(hptype)

    fq  <- .extractFastq(
      x = self$absPath(self$mapInit$bamfile),
      qnames = qnames
    )
    file <- paste(
      "hap", hptype, self$getLrdType(), self$getLrMapper(),
      "lim" %<<% as.character(100*abs(attr(self$getHapList(hptype), "limit"))),
      "n" %<<% length(fq),
      "fastq", "gz", sep = ".")
    out <- .fileDeleteIfExists(file.path(dir, file))
    ShortRead::writeFastq(fq, out, compress = TRUE)
    self$mapIter$`0`[[hptype]] = structure(
      list(
        dir     = self$relPath(dir),
        reads   = self$relPath(out),
        ref     = NULL,
        bamfile = NULL,
        pileup  = NULL,
        conseq  = NULL,
        seqpath = NULL,
        params  = NULL,
        tag     = NULL
      ),
      class = c("mapIter", "list")
    )
  }

  return(invisible(self))
})


DR2S_$set(
  "public", "runGetPartitionedConsensus",
  function(opts = list()) {
    ## Check if reporting is already finished and exit safely
    if (.checkReportStatus(self)) return(invisible(self))

    # Construct consensus from initial mapping with the clustered reads
    flog.info(" Construct consensus sequences using the mapInit reference",
              name = "info")
    bamfile <- self$absPath(self$mapInit$bamfile)
    ref <- if (self$getPartSR()) self$mapInit$SR1$conseq else self$getRefSeq()
    mat <- .msaFromBam(bamfile, ref, paddingLetter = ".")
    foreach(hptype = self$getHapTypes()) %do% {
      flog.info("  Constructing a consensus for haplotype %s ...",
                hptype, name = "info")
      readIds <- self$getHapList(hptype)
      cmat <- .extractIdsFromMat(mat, readIds)
      conseqName <- "consensus.mapIter.0." %<<% hptype
      conseq <- conseq(cmat, name = conseqName, type = "prob",
                       excludeGaps = FALSE)
      seqpath <- self$absPath(
        file.path(self$mapIter$`0`[[hptype]]$dir, conseqName %<<% ".fa"))
      Biostrings::writeXStringSet(
        Biostrings::DNAStringSet(gsub("[-+]", "", conseq)),
        seqpath)
      self$mapIter$`0`[[hptype]] = structure(
        list(
          dir     = self$mapIter$`0`[[hptype]]$dir,
          reads   = self$mapIter$`0`[[hptype]]$reads,
          ref     = "mapIter0",
          bamfile = NULL,
          pileup  = NULL,
          conseq  = conseq,
          seqpath = self$relPath(seqpath),
          params  = NULL,
          tag     = "mapIter0"
        ),
        class = c("mapIter", "list")
      )
      NULL
    }
  })


## Method: mapIter ####
#' @export
mapIter.DR2S <- function(x,
                         opts = list(),
                         iterations = 1,
                         minBaseQuality = 3,
                         minMapq = 0,
                         maxDepth = 1e4,
                         minNucleotideDepth = 3,
                         includeDeletions = TRUE,
                         includeInsertions = TRUE,
                         gapSuppressionRatio = 2/5,
                         force = FALSE,
                         fullname = FALSE,
                         plot = TRUE) {
  x$runMapIter(opts = opts,
               iterations = iterations,
               minBaseQuality = minBaseQuality,
               minMapq = minMapq,
               maxDepth = maxDepth,
               minNucleotideDepth = minNucleotideDepth,
               gapSuppressionRatio = gapSuppressionRatio,
               force = force,
               fullname = fullname,
               plot = plot)
  invisible(x)
}

DR2S_$set(
  "public", "runMapIter",
  function(opts = list(),
           iterations = 1,
           threshold = NULL,
           minBaseQuality = 3,
           minMapq = 0,
           maxDepth = 1e4,
           minNucleotideDepth = 3,
           gapSuppressionRatio = 1/4,
           force = FALSE,
           fullname = FALSE,
           plot = TRUE) {

    # # debug
    # self <- dr2s
    # opts = list()
    # minBaseQuality = 3
    # minMapq = 0
    # maxDepth = 1e4
    # minNucleotideDepth = 3
    # includeDeletions = TRUE
    # includeInsertions = TRUE
    # gapSuppressionRatio = 2/5
    # force = FALSE
    # fullname = TRUE
    # plot = TRUE
    # iterations = 1
    # ##

    flog.info("Step 2: MapIter ...", name = "info")
    flog.info(" Iterative mapping of partitioned longreads", name = "info")

    ## Check if reporting is already finished and exit safely
    if (.checkReportStatus(self)) return(invisible(self))

    ## Overide default arguments
    args <- self$getOpts("mapIter")
    if (!is.null(args)) {
      env  <- environment()
      list2env(args, envir = env)
    }

    if (is.null(threshold)) {
      threshold <- self$getThreshold()
    }
    hptypes <- self$getHapTypes()
    iterations <- self$getIterations()
    baseoutdir   <- self$absPath("mapInit")

    includeInsertions = ifelse(self$getPartSR(), FALSE, TRUE)
    callInsertions = ifelse(self$getPartSR(), FALSE, TRUE)

    ## Mapper
    mapFun <- self$getLrMapFun()

    for (iteration in seq_len(iterations)) {
      flog.info(" Iteration %s of %s", iteration, iterations, name = "info")

      iterationC <- toString(iteration)
      prevIteration <- self$mapIter[[toString(iteration - 1)]]

      foreach(hptype = hptypes) %do% {
        reffile  <- self$absPath(prevIteration[[hptype]]$seqpath)
        refseq   <- prevIteration[[hptype]]$conseq
        readfile <- self$absPath(prevIteration[[hptype]]$reads)
        readtype <- self$getLrdType()
        allele   <- "mapIter" %<<% iterationC
        optsname <- sprintf("%s", hptype)
        refname  <- prevIteration[[hptype]]$ref
        outdir   <- file.path(baseoutdir, hptype)
        mapLabel <- "mapIter"
        maptag   <- paste(mapLabel,
                          paste0(litArrows(c(iteration, hptype,
                                             readtype, self$getLrMapper())),
                                 collapse = " "))
        flog.info("  Map partitioned longreads of haplotype %s", hptype,
                  name = "info")

        pileup <- mapReads(
          maptag = maptag, reffile = reffile, readfile = readfile,
          allele = allele, readtype = readtype, opts = opts, refname = refname,
          optsname = optsname, force = force, maxDepth = maxDepth,
          outdir = outdir, minMapq = minMapq, threshold = threshold,
          minBaseQuality = minBaseQuality, clean = clean,
          minNucleotideDepth = minNucleotideDepth,
          includeDeletions = TRUE, includeInsertions = includeInsertions,
          callInsertions = callInsertions, mapFun = mapFun,
          distributeGaps = TRUE, refseq = refseq)

        # ## Construct consensus sequence
        flog.info("   Constructing consensus ...", name = "info")
        conseqName <- "consensus." %<<% sub(".bam", "",
                                            basename(pileup$bamfile))
        conseq      <- conseq(pileup, name = conseqName, type = "prob",
                              excludeGaps = TRUE,
                              gapSuppressionRatio = gapSuppressionRatio)
        seqpath     <- file.path(outdir, conseqName %<<% ".fa")

        Biostrings::writeXStringSet(
          Biostrings::DNAStringSet(gsub("[-+]", "", conseq)),
          seqpath
        )

        ## Initialize structure
        self$mapIter[[iterationC]][[hptype]] = structure(
          list(
            dir     = self$relPath(outdir),
            reads   = self$relPath(readfile),
            ref     = conseqName,
            bamfile = self$relPath(pileup$bamfile),
            pileup  = pileup,
            conseq  = conseq,
            seqpath = self$relPath(seqpath),
            params  = list(gapSuppressionRatio = gapSuppressionRatio),
            tag     = maptag
          ),
          class = c("mapIter", "list")
        )

      }
    }
    if (plot) {
      flog.info(" Plot MapIter summary", name = "info")
      ## Coverage and base frequency
      plotlist <- foreach(iteration = seq_len(self$getIterations())) %do% {
        self$plotmapIterSummary(thin = 0.1, width = 4, iteration = iteration,
                                drop.indels = TRUE)
      }
      p <- cowplot::plot_grid(plotlist = plotlist, nrow = self$getIterations())
      cowplot::save_plot(p, filename = self$absPath("plot.MapIter.pdf"),
                         base_width = 12*length(hptypes),
                         title = paste(self$getLocus(),
                                       self$getSampleId(), sep = "." ),
                         base_height = 3*self$getIterations())
      cowplot::save_plot(p, filename = self$absPath(".plots/plot.MapIter.svg"),
                         base_width = 12*length(hptypes),
                         base_height = 3*self$getIterations())

    }
    createIgvConfigs(self,map = "mapIter", open = "FALSE")

    invisible(self)
  })

#' @export
print.mapIter <- function(x, ...) {
  msg <- sprintf("An object of class '%s'.\\n", class(x)[1])
  bamf <- ifelse(is.null(x$bamfile), " no bamfile", basename(x$bamfile %||% ""))
  msg <- sprintf("%s [Dir] %s\\n", msg, x$dir)
  readmsg <- if (is.null(x$reads)) "" else sprintf("[Reads] %s\\n ",
                                                   basename(x$reads))
  refmsg <- if (is.null(x$ref)) "" else sprintf("[Reference] %s\\n ",
                                                basename(x$ref))
  cat(msg %<<% readmsg %<<% refmsg)
}

## Method: partitionShortReads ####
#' @export
partitionShortReads.DR2S <- function(x,
                                     opts = list(),
                                     force = FALSE,
                                     ...) {
  x$runPartitionShortReads(opts = opts,
                           force = force,
                           ...)
  invisible(x)
}
# TODO: look at arguments and make same
DR2S_$set("public", "runPartitionShortReads", function(opts = list(),
                                                       force = FALSE,
                                                       optsname = "",
                                                       minMapq = 0,
                                                       ...) {

  ## debug
  # opts = list()
  # force = FALSE
  # optsname = ""
  # minMapq = 0

  flog.info("Step 3: PartitionShortReads ...", name = "info")
  flog.info(" Partition shortreads based on initial mapping and " %<<%
              "longread clustering", name = "info")

  ## Check if reporting is already finished and exit safely
  if (.checkReportStatus(self)) return(invisible(self))

  ## Overide default arguments
  args <- self$getOpts("partitionSR")
  if (!is.null(args)) {
    env  <- environment()
    list2env(args, envir = env)
  }

  ## exit if no shortreads provided
  if (is.null(self$getConfig("shortreads"))) {
    flog.warn(" Cannot partition shortreads. No shortreads provided",
              name = "info")
    return(invisible(self))
  }

  ## Check if there is a shortread mapping from mapInit and use it.
  ## If not, map to the reference
  if (self$getPartSR()) {
    flog.info(" Found shortread mapping from MapInit", name = "info")
    bamfile <- self$absPath(self$mapInit$SR2$bamfile)
    ref <- self$mapInit$SR1$conseq
    refname <- names(ref)
  } else {
    mapfmt  <- "mapPartSR <%s> <%s> <%s> <%s>"
    maptag  <- sprintf(mapfmt, self$mapInit$SR1$ref, self$getSrdType(),
                       self$getSrMapper(), optstring(opts, optsname))

    flog.warn(" Found no shortread mapping from MapInit", name = "info")
    flog.info(" Map shortreads against provided reference", name = "info")

    ref <- self$getRefSeq()
    refname <- names(ref)
    reffile <- self$getRefPath()
    ## Run mapper
    flog.info("  Indexing ...", name = "info")
    samfile <- self$getSrMapFun(
      reffile  = reffile,
      readfile = self$getShortreads(),
      allele   = refname,
      readtype = self$getSrdType(),
      opts     = opts,
      refname  = "",
      optsname = optsname,
      force    = force,
      outdir   = self$getOutdir()
    )
    ## Run bam - sort - index pipeline
    flog.info("  Indexing ...", name = "info")
    bamfile <- .bamSortIndex(
      samfile = samfile,
      reffile = reffile,
      minMapq = minMapq,
      force = force,
      clean = TRUE
    )
  }
  hptypes <- self$getHapTypes()
  prtMat <- self$partition$mat
  seqs <- lapply(self$partition$hpl, function(x) .getSeqsFromMat(
    as.matrix(prtMat[x,])))

  mats <- lapply(seqs, function(x, names) {
    magrittr::set_colnames(createPWM(x), names)
  }, names = colnames(prtMat))

  # Run partitioning
  srpartition <- getSRPartitionScores(refname, bamfile, mats, cores = "auto")

  ## Assign read to haplotype with highest probability,
  ## i.e. product over probabilities of each haplotype and choose max
  flog.info(" Get highest-scoring haplotype for each read", name = "info")
  srpartition$haplotypes <- scoreHighestSR(srpartition$srpartition,
                                           diffThreshold = 0.001)

  # Write fastqs
  foreach(hptype = hptypes ) %do% {
    srfilenames <- c()
    flog.info(" Write shortread fastq for haplotype %s", hptype, name = "info")
    fqs <- self$getShortreads()
    dontUseReads <- srpartition$haplotypes$read[
      !srpartition$haplotypes$read %in% dplyr::filter(
        srpartition$haplotypes, haplotype == hptype)$read]

    # write fastq's
    foreach(fq = fqs) %do% {
      srFastqHap <- self$absPath(
        file.path(
          self$mapIter$`0`[[hptype]]$dir,
          dot(c(strsplit1(basename(fq), "\\.")[1], hptype, "fastq.gz"))))
      .writePartFq(fq = fq, srFastqHap = srFastqHap,
                   dontUseReads = dontUseReads)
      srfilenames <- c(srfilenames, srFastqHap)
      self$srpartition[[hptype]]$srpartition <- srpartition
    }
    self$srpartition[[hptype]]$SR[[fq]] <- self$relPath(srfilenames)
  }
  return(invisible(self))
})

## Method: mapFinal ####
#' @export
mapFinal.DR2S <- function(x,
                          opts = list(),
                          minBaseQuality = 3,
                          minMapq = 50,
                          maxDepth = 1e5,
                          minNucleotideDepth = 3,
                          includeDeletions = TRUE,
                          includeInsertions = TRUE,
                          force = FALSE,
                          fullname = FALSE,
                          plot = TRUE,
                          createIgv = TRUE,
                          clip = FALSE) {
  x$runMapFinal(opts = opts,
                minBaseQuality = minBaseQuality,
                minMapq = minMapq,
                maxDepth = maxDepth,
                minNucleotideDepth = minNucleotideDepth,
                includeDeletions = includeDeletions,
                includeInsertions = includeInsertions,
                force = force,
                fullname = fullname,
                plot = plot,
                createIgv = createIgv,
                clip = clip)
  invisible(x)
}
DR2S_$set("public", "runMapFinal", function(opts = list(),
                                            minBaseQuality = 3,
                                            threshold = NULL,
                                            minMapq = 50,
                                            maxDepth = 1e5,
                                            minNucleotideDepth = 3,
                                            includeDeletions = TRUE,
                                            includeInsertions = TRUE,
                                            force = FALSE,
                                            fullname = FALSE,
                                            plot = TRUE,
                                            createIgv = TRUE,
                                            clip = FALSE) {

  ## debug
  # opts = list()
  # minBaseQuality = 3
  # minMapq = 50
  # maxDepth = 1e5
  # minNucleotideDepth = 3
  # includeDeletions = TRUE
  # includeInsertions = TRUE
  # force = FALSE
  # fullname = TRUE
  # plot = TRUE
  # clip = TRUE
  # self <- dr2s
  # library(futile.logger)
  # library(foreach)
  # library(rlang)

  flog.info("Step 4: mapFinal ...", name = "info")
  flog.info(" Map shortreads and longreads against refined consensus sequences",
            name = "info")

  ## Check if reporting is already finished and exit safely
  if (!force && .checkReportStatus(self)) return(invisible(self))

  ## Overide default arguments
  args <- self$getOpts("mapFinal")
  if (!is.null(args)) {
    env  <- environment()
    list2env(args, envir = env)
  }
  if (is.null(threshold)) {
    threshold <- self$getThreshold()
  }

  ## stop if no shortreads provided
  if (is.null(self$getConfig("shortreads"))) {
    flog.warn(" Cannot run mapFinal. No shortreads provided", name = "info")
    return(invisible(self))
  }

  reftag       <- "mapFinal"
  outdir       <- .dirCreateIfNotExists(self$absPath(reftag))
  lastIter     <- self$mapIter[[max(names(self$mapIter))]]
  hptypes      <- self$getHapTypes()
  readfilesLR  <- vapply(hptypes, function(x)
    self$absPath(lastIter[[x]]$reads), FUN.VALUE = character(1))
  reffiles     <- vapply(hptypes, function(x)
    self$absPath(lastIter[[x]]$seqpath), FUN.VALUE = character(1))
  refseqs      <- set_names(lapply(hptypes, function(x)
    lastIter[[x]]$conseq), hptypes)
  readfilesSR  <- set_names(lapply(hptypes, function(x)
    self$absPath(unlist(self$srpartition[[x]]$SR))), hptypes)

  self$mapFinal = structure(
    list(
      dir          = self$relPath(outdir),
      sreads       = lapply(readfilesSR,self$relPath),
      lreads       = self$relPath(readfilesLR),
      ref          = self$relPath(reffiles),
      bamfile      = list(),
      pileup       = list(),
      tag          = list(),
      seqpath      = list(),
      homopolymers = NULL
    ), class = c("mapFinal", "list")
  )

  ## Remap long reads to the same reference sequences as short reads
  for (hptype in hptypes) {
    flog.info(" Run mapFinal for haplotype %s", hptype, name = "info" )
    reffile    <- reffiles[[hptype]]
    refseq     <- refseqs[[hptype]]
    mapgroupLR <- "LR" %<<% hptype
    readtype   <- self$getLrdType()
    maptagLR   <- paste("mapFinal", mapgroupLR, readtype, self$getLrMapper(),
                        optstring(opts), sep = ".")
    readfileLR <- readfilesLR[[hptype]]
    flog.info("  Map longreads to consensus", name = "info")

    pileup <- mapReads(
      maptag = maptagLR, reffile = reffile, readfile = readfileLR,
      threshold = threshold, allele = mapgroupLR, readtype = readtype,
      opts = opts, outdir = outdir, minMapq = minMapq, refname = hptype,
      optsname = optstring(opts), minBaseQuality = minBaseQuality,
      maxDepth = maxDepth, minNucleotideDepth = minNucleotideDepth,
      force = force, includeDeletions = includeDeletions, clean = TRUE,
      includeInsertions = includeInsertions,  mapFun = self$getLrMapFun(),
      callInsertions = FALSE, distributeGaps = TRUE, refseq = refseq)

    self$mapFinal$bamfile[[mapgroupLR]] = pileup$bamfile
    if (createIgv)
      self$mapFinal$igv[[mapgroupLR]] <- createIgvJsFiles(
        reffile, pileup$bamfile, self$getOutdir(), sampleSize = 100,
      fragmentReads = TRUE)
    self$mapFinal$pileup[[mapgroupLR]] = pileup
    self$mapFinal$tag[[mapgroupLR]] = maptagLR

    # calc new consensus
    cseq <- conseq(pileup$consmat, name = "mapFinal" %<<% hptype,
                   type = "ambig", threshold = 0.2, excludeGaps = FALSE)
    self$mapFinal$seq[[hptype]] <- cseq

    ## Map short reads
    if (!is.null(readfilesSR[[hptype]])) {
      flog.info("  Map shortreads to consensus", name = "info")
      mapgroupSR <- "SR" %<<% hptype
      maptagSR   <- paste("mapFinal", mapgroupSR, self$getLrdType(),
                          self$getSrMapper(),
                          optstring(opts), sep = ".")

      readfiles <- readfilesSR[[hptype]]
      readtype <- self$getSrdType()

      pileup <- mapReads(
        maptag = maptagSR, reffile = reffile,  readfile = readfiles,
        threshold = threshold, allele = mapgroupSR, readtype = readtype,
        opts = opts, outdir = outdir, minMapq = minMapq,
        optsname = optstring(opts), minBaseQuality = minBaseQuality + 10,
        maxDepth = maxDepth, minNucleotideDepth = minNucleotideDepth,
        force = force, includeDeletions = includeDeletions, clean = TRUE,
        includeInsertions = includeInsertions,  mapFun = self$getSrMapFun(),
        callInsertions = TRUE, distributeGaps = TRUE, refseq = refseq)
      # calc new consensus
      cseq <- conseq(pileup$consmat, name = "mapFinal" %<<% hptype,
                     type = "ambig", threshold = 0.2, excludeGaps = TRUE)

      self$mapFinal$bamfile[[mapgroupSR]] <- self$relPath(pileup$bamfile)
      if (createIgv)
        self$mapFinal$igv[[mapgroupSR]] <- createIgvJsFiles(
          reffile, pileup$bamfile, self$getOutdir(), sampleSize = 100)
      self$mapFinal$pileup[[mapgroupSR]] = pileup
      self$mapFinal$tag[[mapgroupSR]] = maptagSR
      self$mapFinal$seq[[hptype]] <- cseq
    }

  }

  if (plot) {
    ## Coverage and base frequency
    if (!is.null(self$mapFinal$sreads$A)) {
      readtypes <- c("LR", "SR")
    } else {
      readtypes <- c("LR")
    }
    plotlist <- foreach(readtype = readtypes) %do% {
      self$plotmapFinalSummary(iteration = "final", readtype = readtype,
                               thin = 0.25, width = 20)
    }
    p <- cowplot::plot_grid(plotlist = plotlist, nrow = 2, labels = readtypes)
    cowplot::save_plot(p, filename = self$absPath("plot.MapFinal.pdf"),
                       base_width = 12*length(hptypes),
                       title     = paste(self$getLocus(), self$getSampleId(),
                                         sep = "." ),
                       base_height = 3*length(readtypes))
    cowplot::save_plot(p, filename = self$absPath(".plots/plot.MapFinal.svg"),
                       base_width = 12*length(hptypes),
                       base_height = 3*length(readtypes))
  }

  return(invisible(self))
})

#' @export
print.mapFinal <- function(x, ...) {
  msg  <- sprintf("An object of class '%s'.\n", class(x)[1])
  bamf <- comma(basename(unlist(x$bamfile) %||% ""))
  msg <- sprintf(
    "%s [Dir] %s\n [Longreads] %s\n [Shortreads] %s\n [References] %s\n
    [Bamfile] %s",
    msg, x$dir,
    comma(basename(unlist(x$lreads))),
    comma(basename(unlist(x$sreads))),
    comma(basename(unlist(x$ref))),
    bamf
  )
  cat(msg)
}

## Method: runPipeline ####
DR2S_$set("public", "runPipeline", function() {
  steps_ <- self$getPipeline()
  while (length(steps_) > 0) {
    step <- steps_[1]
    self$run_(step)
    steps_ <- steps_[-1]
  }
  self
})


#' ## Method:  polish ####
#'
#' #' @export
#' polish.DR2S <- function(x) {
#'   flog.info("Step 6: Infer problematic positions and consensus sequence ...",
#'   name = "info")
#'   Sys.sleep(1)
#'   x$polish()
#'   message("  Done!\n")
#'   invisible(x)
#' }
#'
#' DR2S_$set("public", "polish", function() {
#'   self <- polish(self)
#'   return(invisible(self))
#' }
#' ## Method: report  ####
#'
#' #' @export
#' mapFinal.DR2S <- function(x) {
#'   flog.info("Step 7: report consensus sequences and problematic positions",
#'   name = "info")
#'   Sys.sleep(1)
#'   x$report()
#'   message("  Done!\n")
#'   invisible(x)
#' }
#'
#' DR2S_$set("public", "report", function() {
#'   self <- report(self)
#'   return(invisible(self))
#' }
