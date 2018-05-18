# Method: MapInit ####
#' @export
mapInit.DR2S <- function(x,
                         opts = list(),
                         optsname = "",
                         partSR = TRUE,
                         pct = 100,
                         threshold = 0.20,
                         min_base_quality = 3,
                         min_mapq = 50,
                         max_depth = 1e4,
                         min_nucleotide_depth = 3,
                         include_deletions = TRUE,
                         include_insertions = TRUE,
                         microsatellite = FALSE,
                         force = FALSE,
                         fullname = TRUE,
                         filterScores = TRUE,
                         forceMapping = FALSE,
                         plot = TRUE) {
  x$runMapInit(opts = opts,
               optsname = optsname,
               partSR = partSR,
               pct = pct,
               threshold = threshold,
               min_base_quality = min_base_quality,
               min_mapq = min_mapq,
               max_depth = max_depth,
               min_nucleotide_depth = min_nucleotide_depth,
               include_deletions = include_deletions,
               include_insertions = include_insertions,
               microsatellite = microsatellite,
               force = force,
               fullname = fullname,
               filterScores = filterScores,
               forceMapping = forceMapping,
               plot = plot)
  invisible(x)
}

DR2S_$set("public", "runMapInit", function(opts = list(),
                                           optsname = "",
                                           partSR = TRUE,
                                           pct = 100,
                                           threshold = 0.20,
                                           min_base_quality = 3,
                                           min_mapq = 50,
                                           max_depth = 1e4,
                                           min_nucleotide_depth = 3,
                                           include_deletions = TRUE,
                                           include_insertions = TRUE,
                                           microsatellite = FALSE,
                                           force = FALSE,
                                           fullname = TRUE,
                                           filterScores = TRUE,
                                           forceMapping = FALSE,
                                           plot = TRUE) {

  flog.info("Step 0: mapInit ...", name = "info")

  # # debug
  # opts = list()
  # partSR = TRUE
  # optsname = ""
  # pct = 100
  # threshold = 0.20
  # min_base_quality = 3
  # min_mapq = 0
  # max_depth = 1e4
  # min_nucleotide_depth = 3
  # include_insertions = TRUE
  # include_deletions = TRUE
  # force = FALSE
  # fullname = TRUE
  # plot = TRUE
  # microsatellite = TRUE
  # forceMapping = FALSE
  # filterScores = FALSE
  # library(ggplot2)
  # library(foreach)
  # library(futile.logger)
  # library(cowplot)
  # self <-dr2s
  
  stopifnot(pct > 0 && pct <= 100)
  ## TODO recode_header useless so
  recode_header <- FALSE
  ## Overide default arguments
  args <- self$getOpts("mapInit")
  if (!is.null(args)) {
    env  <- environment()
    list2env(args, envir = env)
  }

  microsatellite  <- self$getMicrosatellite()
  partSR          <- self$getPartSR()
  forceMapping <- self$getForceMapping()
  filterScores    <- self$getFilterScores()
  outdir          <- .dirCreateIfNotExists(self$absPath("mapInit"))
  .dirCreateIfNotExists(path = file.path(self$absPath(".plots")))

  if (recode_header) {
    stopifnot(
      recode_fastq_header(self$getLongreads()),
      recode_fastq_header(self$getShortreads()[1]),
      recode_fastq_header(self$getShortreads()[2])
    )
  }

  if (partSR) {
    mapfmt  <- "mapInit1 <%s> <%s> <%s> <%s>"
    maptag  <- sprintf(mapfmt, self$getReference(), self$getSrdType(),
                       self$getSrMapper(), optstring(opts, optsname))
    readfile <- self$getShortreads()
    # if (threshold != self$getThreshold()) {
    #   self$setThreshold(threshold)
    # }

    flog.info(" Map shortreads to provided reference", name = "info")

    ## Fetch mapper
    map_fun <- self$getSrMapFun()
    ## Run mapper
    flog.info("  Mapping ...", name = "info")
    samfile <- map_fun(
      reffile  = self$getRefPath(),
      readfile = readfile,
      allele   = self$getReference(),
      readtype = self$getSrdType(),
      opts     = opts,
      refname  = "",
      optsname = optsname,
      force    = force,
      outdir   = outdir
    )

    if (filterScores) {
      flog.info(" Filter reads with low alignment scores", name = "info")
      ## Run bam - sort - index pipeline
      bamfile <- bam_sort_index(samfile, self$getRefPath(), pct / 100,
                                min_mapq, force = force, clean = TRUE)
      ## Filter Reads
      bam <- Rsamtools::scanBam(bamfile,
                                param = Rsamtools::ScanBamParam(
                                  tag = "AS",
                                  what = c("qname", "pos", "cigar")
                                ))[[1]]
      readfilter <- .filterReads(bam = bam, preserve_ref_ends = TRUE)
      .fileDeleteIfExists(bamfile)

      flog.info(" Write new shortread fastqs to file", name = "info")
      fqs <- self$getShortreads()
      fqdir <- .dirCreateIfNotExists(file.path(outdir,self$getSrdType()))
      # write fastq's
      readfile <- c()
      readfile <- foreach(fq = fqs, .combine = c) %do% {
        srFastqHap = file.path(fqdir, basename(fq))
        write_part_fq(fq = fq, srFastqHap = srFastqHap,
                      dontUseReads = readfilter)
        srFastqHap
      }
      # set new shortread directory
      self$setConfig("filteredShortreads", self$relPath(fqdir))

      flog.info(" Map filtered shortreads to provided reference", name = "info")
      ## Rerun mapper
      flog.info("  Mapping ...", name = "info")
      samfile <- map_fun(
        reffile  = self$getRefPath(),
        readfile = readfile,
        allele   = self$getReference(),
        readtype = self$getSrdType(),
        opts     = opts,
        refname  = "",
        optsname = optsname,
        force    = force,
        outdir   = outdir
      )
    }

    ## Run bam - sort - index pipeline
    flog.info("  Indexing ...", name = "info")
    bamfile <- bam_sort_index(
      samfile = samfile,
      reffile = self$getRefPath(),
      sample = pct / 100,
      min_mapq = min_mapq,
      force = force,
      clean = TRUE
    )

    ## Calculate pileup from graphmap produced SAM file
    flog.info("  Piling up ...", name = "info")
    pileup <- Pileup(
      bamfile,
      self$getThreshold(),
      max_depth,
      min_base_quality = min_base_quality,
      min_mapq = min_mapq,
      min_nucleotide_depth = min_nucleotide_depth,
      include_deletions = include_deletions,
      include_insertions = include_insertions
    )

    if (include_insertions && is.null(ins(pileup$consmat))) {
      pileup <- pileup_include_insertions(x = pileup, threshold = 0.1)
    }
    ## distribution of gaps not necessary for shortreads (need to check if also
    ## true for homopolymer regions > 15)
    # pileup$consmat <- .distributeGaps(pileup$consmat, removeError = FALSE)

    # debug
    # pileup$consmat[which(pileup$consmat[,6] > 15),]
    if (max(rowSums(pileup$consmat)) /
        quantile(rowSums(pileup$consmat), 0.75) > 5) {
      flog.warn(" Shortreads seem corrupted or the reference is bad!",
                name = "info")
      maxCov <- max(rowSums(pileup$consmat))
      q75Cov <- quantile(rowSums(pileup$consmat), .75)
      flog.warn(paste0("   Maximum of coverage %s / 75%% quantile %s: %s > 5.",
                       " No equal distribution of coverage!",
                       " Have a look at the mapInit plot"),
                maxCov, q75Cov, maxCov/q75Cov, name = "info")
      if (!forceMapping) {
        file <- self$absPath(paste0("plot.MapInit.SR.",
                                    sub("bam$", "pdf", usc(basename(bamfile)))))
        plt <- plot_pileup_coverage(
          x = pileup,
          thin = 0.25,
          width = 2,
          label = self$getMapTag("init", "SR"),
          drop.indels = TRUE
        )
        flog.error(paste(" Aborting. If you want to force processing set ",
                         "forceMapping = TRUE in DR2S object initialisation",
                         " "),
                   name = "info")
        suppressWarnings(ggsave(file, plt, width = 12, height = 10,
                                onefile = TRUE,
                                title = paste(self$getLocus(),
                                              self$getSampleId(), sep = ".")))

        stop("Shortreads probably of bad quality. Bad coverage distribution.
             Run with forceMapping = TRUE to force processing.")
      } else {
        flog.warn(" Continue. Be aware that resulsts may not be correct!!",
                  name = "info")
      }
    }

    # calc initial consensus
    flog.info(" Construct initial consensus from shortreads", name = "info")

    conseq <- conseq(pileup$consmat, name = "mapInit", type = "prob",
                     threshold = self$getThreshold(), force_exclude_gaps = TRUE)
    conseq_name <- paste0("Init.consensus.",
                          sub(".sam.gz", "", basename(samfile)))
    conseqpath  <- file.path(outdir, paste0(conseq_name, ".fa"))
    Biostrings::writeXStringSet(
      Biostrings::DNAStringSet(gsub("[-+]", "N", conseq)),
      conseqpath)

    if (microsatellite) {
      mapfmt  <- "mapInit1.2 <%s> <%s> <%s> <%s>"
      maptag  <- sprintf(mapfmt, conseq_name, self$getSrdType(),
                         self$getSrMapper(), optstring(opts, optsname))
      readfile = self$getShortreads()
      if (threshold != self$getThreshold()) {
        self$setThreshold(threshold)
      }

      flog.info(" Refine microsatellites or repeats by extending the reference",
                name = "info")
      flog.info(" Remap shortreads to initial consensus from shortreads",
                name = "info")

      ## Fetch mapper
      map_fun <- self$getSrMapFun()
      ## Run mapper
      flog.info("  Mapping ...", name = "info")
      samfile <- map_fun(
        reffile  = conseqpath,
        readfile = readfile,
        allele   = conseq_name,
        readtype = self$getSrdType(),
        opts     = opts,
        refname  = "",
        optsname = optsname,
        force    = force,
        outdir   = outdir
      )

      ## Run bam - sort - index pipeline
      flog.info("  Indexing ...", name = "info")
      bamfile <- bam_sort_index(
        samfile = samfile,
        reffile = conseqpath,
        sample = pct / 100,
        min_mapq = min_mapq,
        force = force,
        clean = TRUE
      )

      ## Calculate pileup from graphmap produced SAM file
      flog.info("  Piling up ...", name = "info")
      pileup <- Pileup(
        bamfile,
        self$getThreshold(),
        max_depth,
        min_base_quality = min_base_quality,
        min_mapq = min_mapq,
        min_nucleotide_depth = min_nucleotide_depth,
        include_deletions = include_deletions,
        include_insertions = include_insertions
      )

      if (include_insertions && is.null(ins(pileup$consmat))) {
        pileup <- pileup_include_insertions(pileup, threshold = 0.1)
      }

      # Infer initial consensus
      flog.info(paste0(" Construct second consensus from shortreads ",
                       "with refined repeats"), name = "info")
      conseq <- conseq(pileup$consmat, name = "mapInit1.2", type = "prob",
                       threshold = 0.2, force_exclude_gaps = TRUE)
      conseq_name <- paste0("Init.consensus.2",
                            sub(".sam.gz", "", basename(samfile)))
      conseqpath  <- file.path(outdir, paste0(conseq_name, ".fa"))
      Biostrings::writeXStringSet(
        Biostrings::DNAStringSet(gsub("[-+]", "N", conseq)),
        conseqpath)
    }

    mapInitSR1 = structure(
      list(
        reads   = self$relPath(self$getShortreads()),
        bamfile = self$relPath(bamfile),
        pileup  = pileup,
        tag     = maptag,
        conseq  = conseq,
        seqpath = self$relPath(conseqpath),
        ref     = conseq_name
      ),
      class  = c("mapInit", "list")
    )

    ## Second mapping to infer polymorphic positions
    ## from same reference as longreads
    mapfmt  <- "mapInit2 <%s> <%s> <%s> <%s>"
    maptag  <- sprintf(mapfmt, mapInitSR1$ref,self$getSrdType(),
                       self$getSrMapper(), optstring(opts, optsname))
    readfile <- self$getShortreads()

    flog.info(" Remap shortreads to consensus for SNP calling", name = "info")

    ## Fetch mapper
    map_fun <- self$getSrMapFun()
    ## Run mapper
    flog.info("  Mapping ...", name = "info")
    samfile <- map_fun(
      reffile  = self$absPath(mapInitSR1$seqpath),
      readfile = readfile,
      allele   = mapInitSR1$ref,
      readtype = self$getSrdType(),
      opts     = opts,
      refname  = "",
      optsname = optsname,
      force    = force,
      outdir   = outdir
    )

    ## Run bam - sort - index pipeline
    flog.info("  Indexing ...", name = "info")
    bamfile <- bam_sort_index(
      samfile = samfile,
      reffile = self$absPath(mapInitSR1$seqpath),
      sample = pct / 100,
      min_mapq = min_mapq,
      force = force,
      clean = TRUE
    )

    ## Calculate pileup from graphmap produced SAM file
    flog.info("  Piling up ...", name = "info")
    pileup <- Pileup(
      bamfile,
      self$getThreshold(),
      max_depth,
      min_base_quality = min_base_quality,
      min_mapq = min_mapq,
      min_nucleotide_depth = min_nucleotide_depth,
      include_deletions = include_deletions,
      include_insertions = include_insertions
    )

    mapInitSR2 = structure(
      list(
        reads   = self$relPath(self$getShortreads()),
        bamfile = self$relPath(bamfile),
        pileup  = pileup,
        tag     = maptag,
        conseq  = conseq,
        seqpath = self$relPath(conseqpath),
        ref     = conseq_name
      ),
      class  = c("mapInit", "list")
    )
  }

  if (exists("mapInitSR1")) {
    flog.info(" Map longreads to consensus for clustering",
              name = "info")
    reffile <- self$absPath(mapInitSR1$seqpath)
    refseq <- mapInitSR1$conseq
    allele <- mapInitSR1$ref
  } else {
    flog.info(" Map longreads to provided reference for clustering",
              name = "info")
    reffile <- self$getRefPath()
    refseq <- self$getRefSeq()
    allele <- self$getReference()
  }

  mapfmt  <- "mapInit <%s> <%s> <%s> <%s>"
  maptag  <- sprintf(mapfmt, allele, self$getLrdType(), self$getLrMapper(),
                     optstring(opts, optsname))
  if (threshold != self$getThreshold()) {
    self$setThreshold(threshold)
  }
  ## Fetch mapper
  map_fun <- self$getLrMapFun()
  ## Run mapper
  flog.info("  Mapping ...", name = "info")
  samfile <- map_fun(
    reffile  = reffile,
    readfile = self$getLongreads(),
    allele   = allele,
    readtype = self$getLrdType(),
    opts     = opts,
    refname  = "",
    optsname = optsname,
    force    = force,
    outdir   = outdir
  )

  ## Run bam - sort - index pipeline
  flog.info("  Indexing ...", name = "info")
  bamfile <- bam_sort_index(
    samfile = samfile,
    reffile = reffile,
    sample = pct / 100,
    min_mapq = min_mapq,
    force = force,
    clean = TRUE
  )

  ## Calculate pileup from graphmap produced SAM file
  flog.info("  Piling up ...", name = "info")
  pileup <- Pileup(
    bamfile,
    self$getThreshold(),
    max_depth,
    min_base_quality = min_base_quality,
    min_mapq = min_mapq,
    min_nucleotide_depth = min_nucleotide_depth,
    include_deletions = TRUE,
    include_insertions = FALSE
  )

  pileup$consmat <- .distributeGaps(pileup$consmat, bamfile, reference = refseq,
                                    removeError = TRUE)

  self$mapInit = structure(
    list(
      reads   = self$relPath(self$getLongreads()),
      bamfile = self$relPath(bamfile),
      pileup  = pileup,
      tag     = maptag,
      SR1     = NULL,
      SR2     = NULL
    ),
    class  = c("mapInit", "list")
  )

  if (partSR) {
    self$mapInit$SR1 <- mapInitSR1
    self$mapInit$SR2 <- mapInitSR2
  }

  run_igv(self, map = "mapInit", open_now = "FALSE")

  if (plot) {
    flog.info(" Plot MapInit summary ", name = "info")
    ## Coverage and frequency of minor alleles
    p <- self$plotmapInitSummary(
      thin = 0.25,
      width = 2
    )
    plotRows <- ifelse(partSR, 2, 1)
    cowplot::save_plot(self$absPath("plot.MapInit.pdf"), plot = p, 
              ncol = 1, nrow = plotRows,
              base_aspect_ratio = as.numeric(paste(5, plotRows, sep = ".")),
              title = paste(self$getLocus(), self$getSampleId(), sep = "." ))
    cowplot::save_plot(self$absPath(".plots/plot.MapInit.svg"), plot = p, 
              ncol = 1, nrow = plotRows,
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
                                    threshold = NULL,
                                    skip_gap_freq = 2/3,
                                    dist_alleles = NULL,
                                    noGapPartitioning = FALSE,
                                    select_alleles_by = "count",
                                    plot = TRUE,
                                    ...) {
  x$runPartitionLongReads(threshold = threshold,
                          skip_gap_freq = skip_gap_freq,
                          noGapPartitioning = noGapPartitioning,
                          select_alleles_by = select_alleles_by,
                          dist_alleles = dist_alleles,
                          plot = plot)
  x$runSplitLongReadsByHaplotype(plot = plot)
  x$runExtractLongReads()
  invisible(x)
}

DR2S_$set("public",
          "runPartitionLongReads",
          function(threshold = NULL,
                   skip_gap_freq = 2/3,
                   dist_alleles = NULL,
                   noGapPartitioning = FALSE,
                   select_alleles_by = "count",
                   plot = TRUE) {
  # debug
  # threshold = NULL
  # skip_gap_freq = 2/3
  # dist_alleles = NULL
  # noGapPartitioning = TRUE
  # select_alleles_by = "count"
  # plot = TRUE
  # self <- dr2s
  # library(futile.logger)

  flog.info("Step 1: PartitionLongReads ...", name = "info")
  flog.info(" Partition longreads into haplotypes", name = "info")

  ## Overide default arguments
  args <- self$getOpts("partition")
  if (!is.null(args)) {
    env  <- environment()
    list2env(args, envir = env)
  }
  if (is.null(threshold)) {
    threshold <- self$getThreshold()
  }
  if (is.null(dist_alleles)) {
    dist_alleles <- self$getDistAlleles()
  }
  assertthat::assert_that(
    self$hasPileup(),
    is.double(skip_gap_freq),
    is.double(threshold),
    assertthat::is.count(dist_alleles),
    is.logical(plot)
  )

  tag  <- self$getMapTag("init", "LR")

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

  if (noGapPartitioning) {
    flog.info(" Use only non-gap positions for clustering", name = "info")
    ppos <- ppos %>%
      dplyr::filter(a1 != "-" & a2 != "-")
  }

  ## Check if already finished because it is a homozygous sample
  if (NROW(ppos) == 0) {
    flog.warn(" No polymorphic positions for clustering! Only single allele?",
              name = "info")
    flog.info(" Entering polish and report pipeline", name = "info")
    return(invisible(finish_cn1(self)))
  }

  mat <- if (tryCatch(
    !is(self$partition, "PartList"),
    error = function(e)
      TRUE
  ) ||
  !(all(ppos$position %in% colnames(self$partition$mat)) &&
    all(colnames(self$partition$mat) %in% ppos$position))) {
    SNPmatrix(bamfile = self$absPath(self$mapInit$bamfile), refseq = refseq,
              polymorphic_positions = ppos)
  } else {
    self$partition$mat
  }

  flog.info(" Partition %s longreads over %s SNPs", NROW(mat), NCOL(mat), 
            name = "info")
  prt <- partition_reads(x = mat,
                         skip_gap_freq = skip_gap_freq,
                         deepSplit = 1,
                         threshold = threshold,
                         dist_alleles = dist_alleles,
                         sort_by = select_alleles_by)
  ## Set sample haplotypes
  self$setHapTypes(levels(as.factor(PRT(prt))))

  # Check if we have only one cluster and finish the pipeline if so
  if (length(self$getHapTypes()) == 1) {
    flog.warn(" Only one allele left!")
    flog.info(" Entering polish and report pipeline", name = "info")
    return(invisible(finish_cn1(self)))
  }

  browse_seqs(SQS(prt),
              file = self$absPath("partition.fa.html"),
              openURL = FALSE)

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

# self <- dr2s
DR2S_$set("public", "runSplitLongReadsByHaplotype", function(plot = TRUE) {

  flog.info(" Split partitioned longreads by score", name = "info")

  ## Check if reporting is already finished and exit safely
  if (.checkReportStatus(self)) return(invisible(self))
  stopifnot(self$hasPartition())

  ## Overide default arguments
  args <- self$getOpts("split")
  if (!is.null(args)) {
    env  <- environment()
    list2env(args, envir = env)
  }

  prt  <- partition(self$getPartition())
  haplotypes <- levels(prt$haplotype)

  # Set all limits to NULL
  self$setLimits(sapply(haplotypes, function(x) NULL))
  prts <-  lapply(haplotypes, function(x) prt[prt$haplotype == x,])
  names(prts) <- haplotypes
  scores <- lapply(prts, function(x) x$mcoef)

  lmts <- optimal_partition_limits(scores)
  plainlmts <- as.list(lmts$limits$c)
  names(plainlmts) <- lmts$limits$haplotype

  self$setLimits(plainlmts)
  self$partition$lmt <- lmts$plt

  # Get only reads within the limit
  reads <- lapply(names(self$getLimits()), function(x) {
    dplyr::filter(prt, haplotype == x, mcoef >= self$getLimits()[x])
  })
  names(reads) <- names(self$getLimits())
  for (hp in haplotypes) {
    flog.info("  %s: Using %s longreads with score > %.2f",
              hp, nrow(reads[[hp]]), plainlmts[hp], name = "info")
  }

  # results Structure
  resStruct <- lapply(names(reads), function(x) {
    tr <- reads[[x]]
    structure(
      tr$read,
      q = tr$mcoef,
      freq = NROW(tr)/NROW(dplyr::bind_rows(reads)),
      limit = self$getLimits()[[x]]
    )})
  names(resStruct) <- haplotypes
  self$partition$hpl = structure(resStruct, class = c("HapList", "list"))
  if (plot) {
    p <- self$plotPartitionSummary(label = tag, 
                                   limits = unlist(self$getLimits()))
    
    cowplot::save_plot(self$absPath("plot.Partition.pdf"), plot = p, 
              title = paste(self$getLocus(), self$getSampleId(), sep = "." ),
              base_height = 12, base_width = 10
              )
    cowplot::save_plot(self$absPath(".plots/plot.Partition.svg"), plot = p, 
              base_aspect_ratio = 1.2)

              
    outf  <- self$absPath("plot.Sequence")
    ppos <- SNP(self$getPartition())
    names(ppos) <- 1:length(ppos)
    pwm <- lapply(PWM(self$getPartition()), function(pwm) {
      pwm[pwm < 0.1] <- 0
      pwm
    })
    p <- self$plotSeqLogo(ppos, pwm)
    cowplot::save_plot(filename  = self$absPath("plot.Sequence.pdf"),
                    plot      = p,
                    base_width     = 0.4*length(ppos)+1.4,
                    base_height    = 2.5*length(pwm),
                    title     = paste(self$getLocus(), self$getSampleId(), 
                                      sep = "." ),
                    units     = "cm",
                    limitsize = FALSE)
    cowplot::save_plot(filename  = self$absPath(".plots/plot.Sequence.svg"),
                    plot      = p,
                    base_width     = 0.4*length(ppos)+1.4,
                    base_height    = 2.5*length(pwm),
                    units     = "cm",
                    limitsize = FALSE)

  }
  return(invisible(self))
})

#' @export
print.HapList <- function(x, ...) {
  msg <- sprintf("An object of class '%s'.\n", class(x)[1])
  for (haplotype in names(x)) {
    msg <- paste0(msg,
                  sprintf("%s: n %s; frequency %s; limit %s\n",
                          haplotype,
                          length(x[haplotype]),
                          round(attr(x[[haplotype]], "freq"), 3),
                          attr(x[[haplotype]], "limit")))
  }
  cat(msg)
}

#' @export
## ToDo: change this
summary.HapList <- function(object, ....) {
  data.frame(
    n = c(length(object$A), length(object$B)),
    freq = round(c(
      attr(object$A, "freq"), attr(object$B, "freq")
    ), 2),
    limits = c(attr(object$A, "limit"), attr(object$B, "limit")),
    row.names = c("A", "B")
  )
}

# debug
#self <- dr2s
DR2S_$set("public", "runExtractLongReads", function() {

  flog.info(" Extract haplotyped longreads", name = "info")

  ## Check if reporting is already finished and exit safely
  if (.checkReportStatus(self)) return(invisible(self))

  stopifnot(self$hasHapList())

  ## Overide default arguments
  args <- self$getOpts("extract")
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
    file <- paste("hap", hptype, self$getLrdType(), self$getLrMapper(),
                  paste0("lim", 100 * abs(attr(self$getHapList(hptype),
                                               "limit"))),
                  paste0("n", length(fq)),
                  "fastq", "gz", sep = ".")
    out <- .fileDeleteIfExists(file.path(dir, file))
    ShortRead::writeFastq(fq, out, compress = TRUE)
    self$mapIter[["0"]][[hptype]] = structure(
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

## Method: mapIter ####
#' @export
mapIter.DR2S <- function(x,
                         opts = list(),
                         iterations = 1,
                         pct = 100,
                         min_base_quality = 3,
                         min_mapq = 0,
                         max_depth = 1e4,
                         min_nucleotide_depth = 3,
                         include_deletions = TRUE,
                         include_insertions = TRUE,
                         gap_suppression_ratio = 2/5,
                         force = FALSE,
                         fullname = TRUE,
                         plot = TRUE) {
  x$runMapIter(opts = opts,
               iterations = iterations,
               pct = pct,
               min_base_quality = min_base_quality,
               min_mapq = min_mapq,
               max_depth = max_depth,
               min_nucleotide_depth = min_nucleotide_depth,
               include_deletions = include_deletions,
               include_insertions = include_insertions,
               gap_suppression_ratio = gap_suppression_ratio,
               force = force,
               fullname = fullname,
               plot = plot)
  invisible(x)
}

# self <- hla.fq
# self <- dpb1_3
DR2S_$set("public", "runMapIter", function(opts = list(),
                                           iterations = 1,
                                           pct = 100,
                                           min_base_quality = 3,
                                           min_mapq = 0,
                                           max_depth = 1e4,
                                           min_nucleotide_depth = 3,
                                           include_deletions = TRUE,
                                           include_insertions = TRUE,
                                           gap_suppression_ratio = 2/5,
                                           force = FALSE,
                                           fullname = TRUE,
                                           plot = TRUE) {
  # # debug
  # self <- dr2s
  # opts = list()
  # pct = 100
  # min_base_quality = 3
  # min_mapq = 0
  # max_depth = 1e4
  # min_nucleotide_depth = 3
  # include_deletions = TRUE
  # include_insertions = TRUE
  # gap_suppression_ratio = 2/5
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
  hptypes <- self$getHapTypes()
  iterations <- self$getIterations()

  ## Mapper
  map_fun <- self$getLrMapFun()

  # ### HACK: Turn off Distribute Gaps!!! REMOVE AFTER MIC
  # distGaps <- FALSE

  # Construct consensus from initial mapping with the clustered reads
  flog.info(" Construct consensus sequences using the mapInit reference", 
            name = "info")
  bamfile <- self$absPath(self$mapInit$bamfile)
  if (self$getPartSR()) {
    ref <- self$mapInit$SR1$conseq
  } else {
    ref <- self$getRefSeq()
  }
  mat <- msa_from_bam(bamfile, ref, paddingLetter = ".")

  #hptype <- "A"
  foreach(hptype = hptypes) %do% {
    flog.info("  Constructing a consensus for haplotype %s ...", 
              hptype, name = "info")
    readIds <- self$getHapList(hptype)
    cmat <- consmat(t(
      Biostrings::consensusMatrix(mat[readIds], as.prob = FALSE)[
        VALID_DNA(include = "indel"), ]
    ), freq = FALSE)
    conseq_name <- paste0("consensus.mapIter.0.", hptype)
    conseq <- conseq(cmat, name = conseq_name, type = "prob", 
                     exclude_gaps = FALSE)
    seqpath <- self$absPath(
      file.path(self$mapIter[["0"]][[hptype]]$dir, paste0(conseq_name, ".fa")))
    self$mapIter[["0"]][[hptype]]$ref     <- "mapIter0"
    self$mapIter[["0"]][[hptype]]$conseq  <- conseq
    self$mapIter[["0"]][[hptype]]$seqpath <- self$relPath(seqpath)
    self$mapIter[["0"]][[hptype]]$tag     <- "mapIter0"
    Biostrings::writeXStringSet(
      Biostrings::DNAStringSet(gsub("[-+]", "", conseq)),
      seqpath
    )
    NULL
  }

  # iteration <- 1
  for (iteration in 1:iterations) {
    flog.info(" Iteration %s of %s", iteration, iterations, name = "info")

    iterationC <- as.character(iteration)
    prevIteration <- self$mapIter[[as.character(as.numeric(iteration - 1))]]

    # hptype = "B"
    foreach(hptype = hptypes) %do% {
      reftag   <- prevIteration[[hptype]]$ref
      outdir   <- self$absPath(prevIteration[[hptype]]$dir)
      readpath <- self$absPath(prevIteration[[hptype]]$reads)
      refpath  <- self$absPath(prevIteration[[hptype]]$seqpath)
      refseq   <- prevIteration[[hptype]]$conseq

      optsname <- sprintf("%s", hptype)
      mapfmt  <- "mapIter <%s> <%s> <%s> <%s>"
      maptag  <- sprintf(mapfmt, iteration, hptype, self$getLrdType(),
                         self$getLrMapper(), optstring(opts, optsname))

      self$mapIter[[as.character(iteration)]][[hptype]] = structure(
        list(
          dir     = self$relPath(outdir),
          reads   = self$relPath(readpath),
          ref     = self$relPath(refpath),
          bamfile = list(),
          pileup  = list(),
          conseq  = list(),
          seqpath = list(),
          params  = list(gap_suppression_ratio = gap_suppression_ratio),
          tag     = list()
        ),
        class = c("mapIter", "list")
      )

      flog.info("  Map partitioned longreads of haplotype %s", hptype, 
                name = "info")

      ## Run mapper
      flog.info("   Mapping ...", name = "info")
      samfile <- map_fun(
        reffile  = refpath,
        readfile = readpath,
        allele   = paste0("mapIter", iteration),
        readtype = self$getLrdType(),
        opts     = opts,
        refname  = reftag,
        optsname = optsname,
        force    = force,
        outdir   = outdir
      )

      # samfile
      ## Run bam - sort - index pipeline
      flog.info("   Indexing ...", name = "info")
      bamfile <- bam_sort_index(
        samfile,
        refpath,
        pct / 100,
        min_mapq,
        force = force,
        clean = TRUE
      )
      self$mapIter[[iterationC]][[hptype]]$bamfile = self$relPath(bamfile)

      ## Calculate pileup from graphmap produced SAM file
      flog.info("   Piling up ...", name = "info")
      pileup <- Pileup(
        bamfile,
        self$getThreshold(),
        max_depth = max_depth,
        min_base_quality = min_base_quality,
        min_mapq = min_mapq,
        min_nucleotide_depth = min_nucleotide_depth,
        include_deletions = include_deletions,
        include_insertions = include_insertions
      )
      if (include_insertions && is.null(ins(pileup$consmat))) {
        pileup <- pileup_include_insertions(x = pileup, threshold = 0.2)
      }
      pileup$consmat <- .distributeGaps(mat = pileup$consmat,
                                        bamfile = bamfile,
                                        reference = refseq,
                                        removeError = TRUE)
      self$mapIter[[iterationC]][[hptype]]$pileup = pileup

      # ## Construct consensus sequence
      flog.info("   Constructing consensus ...", name = "info")
      conseq_name <- paste0("consensus.", sub(".sam.gz", "", basename(samfile)))
      conseq      <- conseq(pileup, name = conseq_name, type = "prob",
                            exclude_gaps = TRUE,
                            gap_suppression_ratio = gap_suppression_ratio)
      seqpath     <- file.path(outdir, paste0(conseq_name, ".fa"))
      self$mapIter[[iterationC]][[hptype]]$seqpath = self$relPath(seqpath)
      self$mapIter[[iterationC]][[hptype]]$conseq  = conseq
      self$mapIter[[iterationC]][[hptype]]$ref     = conseq_name
      Biostrings::writeXStringSet(
        Biostrings::DNAStringSet(gsub("[-+]", "", conseq)),
        seqpath
      )

      ## Set maptag
      self$mapIter[[iterationC]][[hptype]]$tag = maptag
    }

  }
  if (plot) {
    flog.info(" Plot MapIter summary", name = "info")
    ## Coverage and base frequency
    plotlist <- foreach(iteration = 1:self$getIterations()) %do% {
      self$plotmapIterSummary(thin = 0.1, width = 4, iteration = iteration,
                              drop.indels = TRUE)
    }
    p <- cowplot::plot_grid(plotlist = plotlist, nrow = self$getIterations())
    cowplot::save_plot(p, filename = self$absPath("plot.MapIter.pdf"),
              base_width = 24*length(hptypes),
              title     = paste(self$getLocus(), 
                                self$getSampleId(), sep = "." ),
              base_height = 6*self$getIterations())
    cowplot::save_plot(p, filename = self$absPath(".plots/plot.MapIter.svg"),
              base_width = 24*length(hptypes),
              base_height = 6*self$getIterations())
    
  }
  run_igv(self,map = "mapIter", open_now = "FALSE")

  invisible(self)
})

#' @export
print.mapIter <- function(x, ...) {
  msg <- sprintf("An object of class '%s'.\\n", class(x)[1])
  bamf <- ifelse(is.null(x$bamfile), " no bamfile", basename(x$bamfile %||% ""))
  msg <- sprintf(
    "%s [Dir] %s\\n [Reads] %s\\n [Reference] %s\\n [Bamfile] %s\\n",
    msg,
    x$dir,
    basename(x$reads),
    basename(x$ref),
    bamf
  )
  cat(msg)
}

## Method: partitionShortReads ####
#' @export
partitionShortReads.DR2S <- function(x,
                                     opts = list(),
                                     force = FALSE) {
  x$runPartitionShortReads(opts = opts,
                           force = force)
  invisible(x)
}

DR2S_$set("public", "runPartitionShortReads", function(opts = list(),
                                                       force = FALSE,
                                                       optsname = "",
                                                       pct = 100,
                                                       threshold = 0.20,
                                                       min_mapq = 0) {

  ## debug
  # opts = list()
  # force = FALSE
  # optsname = ""
  # pct = 100
  # threshold = 0.20
  # min_mapq = 0

  flog.info("Step 3: PartitionShortReads ...", name = "info")
  flog.info(paste0(" Partition shortreads based on initial mapping and ", 
                   "longread clustering"), name = "info")

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

    if (threshold != self$getThreshold()) {
      self$setThreshold(threshold)
    }

    flog.warn(" Found no shortread mapping from MapInit", name = "info")
    flog.info(" Map shortreads against provided reference", name = "info")

    ref <- self$getRefSeq()
    refname <- names(ref)
    reffile <- self$getRefPath()

    ## Fetch mapper
    map_fun <- self$getSrMapFun()
    ## Run mapper
    flog.info("  Indexing ...", name = "info")
    samfile <- map_fun(
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
    bamfile <- bam_sort_index(
      samfile = samfile,
      reffile = reffile,
      sample = pct / 100,
      min_mapq = min_mapq,
      force = force,
      clean = TRUE
    )
  }

  hptypes <- self$getHapTypes()
  prt_mat <- self$partition$mat
  seqs <- lapply(self$partition$hpl, function(x) get_seqs_from_mat(
    as.matrix(prt_mat[x,])))
  names(seqs) <- hptypes

  mats <- lapply(seqs, function(x) create_PWM(x))
  mats <- foreach(m = mats) %do% {
    colnames(m) <- colnames(prt_mat)
    m
  }
  names(mats) <- names(seqs)

  # Run partitioning
  srpartition <- get_SR_partition_scores(refname, bamfile, mats, cores = "auto")

  ## Assign read to haplotype with highest probability,
  ## i.e. product over probabilities of each haplotype and choose max
  flog.info(" Get highest-scoring haplotype for each read", name = "info")
  srpartition$haplotypes <- score_highest_SR(srpartition$srpartition,
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
      write_part_fq(fq = fq, srFastqHap = srFastqHap, 
                    dontUseReads = dontUseReads)
      srfilenames <- c(srfilenames, srFastqHap)
      self$srpartition[[hptype]]$srpartition <- srpartition
    }
    self$srpartition[[hptype]]$SR[[fq]] <- self$relPath(srfilenames)
  }
  return(invisible(self))
})

#' @export
print.partitionShortReads <- function(x, ...) {
  msg  <- sprintf("An object of class '%s'.\n", class(x)[1])
  bamf <- paste0(basename(unlist(x$bamfile) %||% ""), collapse = ", ")
  seqp <- paste0(basename(unlist(x$seqpath) %||% ""), collapse = ", ")

  msg <- sprintf(
    "%s [Dir] %s\n [Longreads] %s\n [Shortreads] %s\n [References] %s\n
    [Bamfile] %s\n [Seqpath] %s\n",
    msg, x$dir,
    paste0(basename(unlist(x$reads)), collapse = ", "),
    paste0(basename(x$sreads), collapse = ", "),
    paste0(basename(unlist(x$ref)), collapse = ", "),
    bamf, seqp
  )
  cat(msg)
}

## Method: mapFinal ####
#' @export
mapFinal.DR2S <- function(x,
                          opts = list(),
                          pct = 100,
                          min_base_quality = 3,
                          min_mapq = 50,
                          max_depth = 1e5,
                          min_nucleotide_depth = 3,
                          include_deletions = TRUE,
                          include_insertions = TRUE,
                          force = FALSE,
                          fullname = TRUE,
                          plot = TRUE,
                          clip = FALSE) {
  x$runMapFinal(opts = opts,
                pct = pct,
                min_base_quality = min_base_quality,
                min_mapq = min_mapq,
                max_depth = max_depth,
                min_nucleotide_depth = min_nucleotide_depth,
                include_deletions = include_deletions,
                include_insertions = include_insertions,
                fullname = fullname,
                plot = plot,
                clip = clip)
  invisible(x)
}

DR2S_$set("public", "runMapFinal", function(opts = list(),
                                            pct = 100,
                                            min_base_quality = 3,
                                            min_mapq = 50,
                                            max_depth = 1e5,
                                            min_nucleotide_depth = 3,
                                            include_deletions = TRUE,
                                            include_insertions = TRUE,
                                            force = FALSE,
                                            fullname = TRUE,
                                            plot = TRUE,
                                            clip = FALSE) {

  ## debug
  # self <- dpb1_3
  # opts = list()
  # pct = 100
  # min_base_quality = 3
  # min_mapq = 50
  # max_depth = 1e5
  # min_nucleotide_depth = 3
  # include_deletions = TRUE
  # include_insertions = TRUE
  # force = FALSE
  # fullname = TRUE
  # plot = TRUE
  # clip = FALSE
  # self <- dr2s

  flog.info("Step 4: mapFinal ...", name = "info")
  flog.info(" Map shortreads and longreads against refined consensus sequences", 
            name = "info")

  ## Check if reporting is already finished and exit safely
  if (.checkReportStatus(self)) return(invisible(self))

  ## Overide default arguments
  args <- self$getOpts("mapFinal")
  if (!is.null(args)) {
    env  <- environment()
    list2env(args, envir = env)
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
  readpathsLR  <- sapply(hptypes, function(x) 
    self$absPath(lastIter[[x]]$reads))
  names(readpathsLR) <- hptypes
  refpaths     <- sapply(hptypes, function(x) 
    self$absPath(lastIter[[x]]$seqpath))
  refseqs      <- sapply(hptypes, function(x) 
    lastIter[[x]]$conseq)
  names(refpaths) <- hptypes
  readpathsSR <- lapply(hptypes, function(x) 
    self$absPath(unlist(self$srpartition[[x]]$SR)))
  names(readpathsSR) <- hptypes

  self$mapFinal = structure(
    list(
      dir          = self$relPath(outdir),
      sreads       = lapply(readpathsSR,self$relPath),
      lreads       = self$relPath(readpathsLR),
      ref          = self$relPath(refpaths),
      bamfile      = list(),
      pileup       = list(),
      tag          = list(),
      seqpath      = list(),
      homopolymers = NULL
    ), class = c("mapFinal", "list")
  )


  ## Remap long reads to the same reference sequences as short reads
  # debug
  # hptype = "B"
  for (hptype in hptypes) {
    flog.info(" Run mapFinal for haplotype %s", hptype, name = "info" )
    refpath  <- refpaths[[hptype]]
    refseq  <- refseqs[[hptype]]
    mapgroupLR <- paste0("LR", hptype)
    maptagLR   <- paste("mapFinal", mapgroupLR, self$getLrdType(),
                        self$getLrMapper(),
                        optstring(opts), sep = ".")
    readpathLR <- readpathsLR[[hptype]]
    flog.info("  Map longreads to consensus", name = "info")
    ## Mapper
    map_fun <- self$getLrMapFun()
    ## Run mapper
    flog.info("  Mapping ...", name = "info")
    samfile <- map_fun(
      reffile  = refpath,
      readfile = readpathLR,
      allele   = mapgroupLR,
      readtype = self$getLrdType(),
      opts     = opts,
      refname  = hptype,
      optsname = optstring(opts),
      force    = force,
      outdir   = outdir
    )

    ## Run bam - sort - index pipeline
    flog.info("  Indexing ...", name = "info")
    bamfile <- bam_sort_index(
      samfile,
      refpath,
      pct / 100,
      min_mapq,
      force = force,
      clean = TRUE
    )
    self$mapFinal$bamfile[[mapgroupLR]] = bamfile

    ## Calculate pileup from graphmap produced SAM file
    flog.info("  Piling up ...", name = "info")
    pileup <- Pileup(
      bamfile,
      self$getThreshold(),
      max_depth = max_depth,
      min_base_quality = min_base_quality,
      min_mapq = min_mapq,
      min_nucleotide_depth = min_nucleotide_depth,
      include_deletions = include_deletions,
      include_insertions = TRUE
    )
    pileup$consmat <- .distributeGaps(pileup$consmat,
                                      bamfile,
                                      refseq,
                                      removeError = TRUE)

    self$mapFinal$pileup[[mapgroupLR]] = pileup

    ## Set maptag
    self$mapFinal$tag[[mapgroupLR]] = maptagLR
    # calc new consensus
    cseq <- conseq(pileup$consmat, name = paste0("mapFinal", hptype),
                   type = "ambig", threshold = 0.2, exclude_gaps = FALSE)
    self$mapFinal$seq[[hptype]] <- cseq


    ## Map short reads
    if (!is.null(readpathsSR[[hptype]])) {
      flog.info("  Map shortreads to consensus", name = "info")
      mapgroupSR <- paste0("SR", hptype)
      maptagSR   <- paste("mapFinal", mapgroupSR, self$getLrdType(),
                          self$getSrMapper(),
                          optstring(opts), sep = ".")

      readfiles <- readpathsSR[[hptype]]
      # debug:
      #readfiles <- self$getShortreads()
      ## Mapper
      map_fun <- self$getSrMapFun()

      ## Run mapper
      flog.info("  Mapping ...", name = "info")
      samfile <- map_fun(
        reffile  = refpath,
        readfile = readfiles,
        allele   = paste0(mapgroupSR, ".", self$getLrdType()),
        readtype = self$getSrdType(),
        opts     = opts,
        refname  = hptype,
        optsname = optstring(opts),
        force    = force,
        outdir   = outdir
      )

      if (clip) {
        flog.info("  Trimming softclips and polymorphic ends ...",
                  name = "info")
        ## Run bam - sort - index pipeline
        bamfile <- bam_sort_index(samfile, refpath, pct / 100, min_mapq,
                                  force = force, clean = TRUE)
        ## Trim softclips
        fq <- .trimSoftclippedEnds(bam = Rsamtools::scanBam(bamfile)[[1]],
                                    preserve_ref_ends = TRUE)
        ## Trim polymorphic ends
        fq <- trim_polymorphic_ends(fq)
        ## Write new shortread file to disc
        fqdir  <- .dirCreateIfNotExists(file.path(self$getOutdir(), 
                                                     "mapFinal"))
        fqfile <- paste("sread", hptype, self$getSrMapper(), "trimmed", "fastq",
                        "gz", sep = ".")
        fqout  <- .fileDeleteIfExists(file.path(fqdir, fqfile))
        ShortRead::writeFastq(fq, fqout, compress = TRUE)
        .fileDeleteIfExists(bamfile)
        ## Rerun mapper
        flog.info("   Mapping trimmed short reads against latest consensus ...",
                  name = "info")
        samfile <- map_fun(
          reffile  = refpath,
          readfile = fqout,
          allele   = paste0(mapgroupSR, ".", self$getLrdType()),
          readtype = self$getSrdType(),
          opts     = list(A = 1, B = 4, O = 2),
          refname  = hptype,
          optsname = optstring(opts),
          force    = force,
          outdir   = outdir
        )
        # cleanup
        .fileDeleteIfExists(fqout)
      }

      ## Run bam - sort - index pipeline
      flog.info("  Indexing ...", name = "info")
      bamfile <- bam_sort_index(
        samfile,
        refpath,
        pct / 100,
        min_mapq,
        force = force,
        clean = TRUE
      )
      self$mapFinal$bamfile[[mapgroupSR]] = self$relPath(bamfile)

      ## Calculate pileup from graphmap produced SAM file
      flog.info("  Piling up ...", name = "info")
      pileup <- Pileup(
        bamfile,
        self$getThreshold(),
        max_depth = max_depth,
        min_base_quality = min_base_quality + 10,
        min_mapq = min_mapq,
        min_nucleotide_depth = min_nucleotide_depth,
        include_deletions = include_deletions,
        include_insertions = TRUE
      )

      if (include_insertions && is.null(ins(pileup$consmat))) {
        pileup <- pileup_include_insertions(pileup, threshold = 0.2)
      }
      self$mapFinal$pileup[[mapgroupSR]] = pileup

      ## Set maptag
      self$mapFinal$tag[[mapgroupSR]] = maptagSR
      # calc new consensus
      cseq <- conseq(pileup$consmat, name = paste0("mapFinal", hptype),
                     type = "ambig", threshold = 0.2, exclude_gaps = TRUE)
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
      plotlist <- foreach (readtype = readtypes) %do% {
        self$plotmapFinalSummary(iteration = "final", readtype = readtype,
                                 thin = 0.25, width = 20)
      }
      p <- cowplot::plot_grid(plotlist = plotlist, nrow = 2, labels = readtypes)
      cowplot::save_plot(p, filename = self$absPath("plot.MapFinal.pdf"),
                base_width = 24*length(hptypes),
                title     = paste(self$getLocus(), self$getSampleId(), 
                                  sep = "." ),
                base_height = 6*length(readtypes))
      cowplot::save_plot(p, filename = self$absPath(".plots/plot.MapFinal.svg"),
                base_width = 24*length(hptypes),
                base_height = 6*length(readtypes))
  }
  run_igv(self, map = "mapFinal", open_now = FALSE)

  return(invisible(self))
})

#' @export
print.mapFinal <- function(x, ...) {
  msg  <- sprintf("An object of class '%s'.\n", class(x)[1])
  bamf <- paste0(basename(unlist(x$bamfile) %||% ""), collapse = ", ")
  msg <- sprintf(
    "%s [Dir] %s\n [Longreads] %s\n [Shortreads] %s\n [References] %s\n
    [Bamfile] %s",
    msg, x$dir,
    paste0(basename(unlist(x$lreads)), collapse = ", "),
    paste0(basename(unlist(x$sreads)), collapse = ", "),
    paste0(basename(unlist(x$ref)), collapse = ", "),
    bamf
  )
  cat(msg)
}

## Method: runPipeline ####
DR2S_$set("public", "runPipeline", function() {
  steps_ <- self$getPipeline()
  while (length(steps_) > 0) {
    step <- steps_[1]
    private$run_(step)
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
