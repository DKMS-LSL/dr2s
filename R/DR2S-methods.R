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
                         forceBadMapping = FALSE,
                         plot = TRUE) {
  flog.info("Step 0: Initial mapping ... ", name = "info")
  Sys.sleep(1)
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
               forceBadMapping = forceBadMapping,
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
                                           forceBadMapping = FALSE,
                                           plot = TRUE) {

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
  # forceBadMapping = FALSE
  # filterScores = FALSE
  # library(ggplot2)
  # library(foreach)
  # library(futile.logger)
  # self <- dr2s


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
  forceBadMapping <- self$getForceBadMapping()
  filterScores    <- self$getFilterScores()
  outdir       <- dir_create_if_not_exists(self$absPath("mapInit"))

  if (recode_header) {
    stopifnot(
      recode_fastq_header(self$getLongreads()),
      recode_fastq_header(self$getShortreads()[1]),
      recode_fastq_header(self$getShortreads()[2])
    )
  }
  if (partSR){
    mapfmt  <- "mapInit1 <%s> <%s> <%s> <%s>"
    maptag  <- sprintf(mapfmt, self$getReference(), self$getSrdType(),
                       self$getSrMapper(), optstring(opts, optsname))
    readfile <- self$getShortreads()
    if (threshold != self$getThreshold()) {
      self$setThreshold(threshold)
    }
    ## Fetch mapper
    map_fun <- self$getSrMapFun()
    ## Run mapper
    flog.info(" First mapping of shortreads to provided reference", name = "info")
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
      flog.info(" Filter reads with low alignment scores...", name = "info")
      ## Run bam - sort - index pipeline
      bamfile <- bam_sort_index(samfile, self$getRefPath(), pct / 100,
                                min_mapq, force = force, clean = TRUE)
      ## Filter Reads
      bam = Rsamtools::scanBam(bamfile,
                               param = Rsamtools::ScanBamParam(tag="AS",
                                                               what = c("qname",
                                                                        "pos",
                                                                        "cigar")
                                                               ))[[1]]
      readfilter <- filter_reads(bam = bam, preserve_ref_ends = TRUE)
      file_delete_if_exists(bamfile)

      flog.info(" Write new shortread fastqs to file ...", name = "info")
      fqs <- self$getShortreads()
      fqdir  <- dir_create_if_not_exists(file.path(outdir,self$getSrdType()))
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

      ## Rerun mapper
      flog.info("  Mapping filtered short reads against latest consensus ... ",
                name = "info")
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
      pileup <- pileup_include_insertions(pileup, threshold = 0.1)
    }
    ## distribution of gaps not necessary for shortreads (need to check if also
    ## true for homopolymer regions > 15)
    # pileup$consmat <- .distributeGaps(pileup$consmat, removeError = FALSE)

    # debug
    # pileup$consmat[which(pileup$consmat[,6] > 15),]
    if (max(rowSums(pileup$consmat))/quantile(rowSums(pileup$consmat), 0.75) > 5){
      flog.warn("  Shortreads seems corrupted or the reference is bad!",
                name = "info")
      maxCov <- max(rowSums(pileup$consmat))
      q75Cov <- quantile(rowSums(pileup$consmat), .75)
      flog.warn("   maximum of coverage %s / 75%% quantile %s: %s > 5.
                         No equal distribution of coverage!
                Have a look at the mapInit Plot",
                maxCov, q75Cov, maxCov/q75Cov, name = "info")
      if (!forceBadMapping) {
        gfile <- self$absPath(paste0("plot.MapInit.SR.",
                                  sub("bam$", "pdf", usc(basename(bamfile)))))
        plt <- plot_pileup_coverage(
          x = pileup,
          thin = 0.25,
          width = 2,
          label = self$getMapTag("init", "SR"),
          drop.indels = TRUE
        )
        flog.error("  Aborting. If you want to force processing set forceBadMapping = TRUE in DR2S object initialisation",
                   name = "info")
        suppressWarnings(ggsave(gfile, plt, width = 12, height = 10,
                                onefile = TRUE,
                                title = paste(self$getLocus(),
                                              self$getSampleId(), sep = ".")))

        stop("Shortreads probably of bad quality. Bad coverage distribution.
             Run with forceBadMapping = TRUE to force processing.")
      } else {
        flog.warn("  Continue. Be aware that resulsts may not be correct!!",
                  name = "info")
      }
    }

    # calc initial consensus
    flog.info("  Construct initial consensus from shortreads", name = "info")

    conseq <- conseq(pileup$consmat, "mapInit", "prob",
                     force_exclude_gaps = TRUE,
                     threshold = self$getThreshold())
    conseq_name <- paste0("Init.consensus.", sub(".sam.gz", "",
                                                 basename(samfile)))
    conseqpath     <- file.path(outdir, paste0(conseq_name, ".fa"))
    Biostrings::writeXStringSet(
      Biostrings::DNAStringSet(gsub("[-+]", "N", conseq)),
      conseqpath)

    if (microsatellite){
      mapfmt  <- "mapInit1.2 <%s> <%s> <%s> <%s>"
      maptag  <- sprintf(mapfmt, conseq_name, self$getSrdType(),
                         self$getSrMapper(), optstring(opts, optsname))
      readfile = self$getShortreads()
      if (threshold != self$getThreshold()) {
        self$setThreshold(threshold)
      }
      ## Fetch mapper
      map_fun <- self$getSrMapFun()
      ## Run mapper
      flog.info(" Refine microsatellites or repeats by extending reference ",
                name = "info")
      flog.info("  Second mapping of shortreads to reference from first mapping ",
                name = "info")
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
      flog.info("   Indexing ...", name = "info")
      bamfile <- bam_sort_index(
        samfile = samfile,
        reffile = conseqpath,
        sample = pct / 100,
        min_mapq = min_mapq,
        force = force,
        clean = TRUE
      )

      ## Calculate pileup from graphmap produced SAM file
      flog.info("   Piling up ...", name = "info")
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
      ## distribution of gaps not necessary for shortreads (need to check if also
      ## true for homopolymer regions > 15)
      # pileup$consmat <- .distributeGaps(pileup$consmat, removeError = FALSE)

      # Infer initial consensus
      flog.info("   Construct second consensus from shortreads with refined repeats ...", name = "info")
      conseq <- conseq(pileup$consmat, "mapInit1.2", "prob",
                       force_exclude_gaps = TRUE, threshold = 0.2)
      conseq_name <- paste0("Init.consensus.2", sub(".sam.gz", "",
                                                    basename(samfile)))
      conseqpath     <- file.path(outdir, paste0(conseq_name, ".fa"))
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

    ## Fetch mapper
    map_fun <- self$getSrMapFun()
    ## Run mapper
    flog.info(" Mapping shortreads against mapInit reference for calling SNPs ...", name = "info")
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

    ## check if calling insertions is necessary; Problem: insertions are not
    ## used for clustering in either way
    # if (include_insertions && is.null(ins(pileup$consmat))) {
    #   pileup <- pileup_include_insertions(pileup, threshold = 0.1)
    # }

    ## distribution of gaps seems not necessary for shortreads
    # pileup$consmat <- .distributeGaps(pileup$consmat, removeError = FALSE)

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

  flog.info(" Map longreads against MapInit reference for clustering ...",
            name = "info")
  if (exists("mapInitSR1")){
    reffile <- self$absPath(mapInitSR1$seqpath)
    refseq <- mapInitSR1$conseq
    allele <- mapInitSR1$ref
  } else {
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


  pileup$consmat <- .distributeGaps(pileup$consmat, bamfile, reference = refseq, removeError = TRUE)

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

  if (partSR){
    self$mapInit$SR1 <- mapInitSR1
    self$mapInit$SR2 <- mapInitSR2
  }
  run_igv(self,map = "mapInit", open_now = "FALSE")
  if (plot) {
    flog.info(" Plotting MapInit Summary ...", name = "info")
    ## Coverage and frequency of minor alleles
    gfile <- self$absPath(paste0("plot.MapInit.LR.",
                              sub("bam$", "pdf",
                                  usc(basename(self$mapInit$bamfile)))))
    pdf(file = gfile, width = 12, height = 10, onefile = TRUE,
        title = paste(self$getLocus(), self$getSampleId(), sep = "." ))
    self$plotmapInitSummary(
      thin = 0.25,
      width = 2,
      label = self$getMapTag("init", "LR"),
      readtype = "LR"
    )
    dev.off()

    if (partSR){
      gfile <- self$absPath(paste0("plot.MapInit.SR.",
                                sub("bam$", "pdf",
                                    usc(basename(self$mapInit$SR2$bamfile)))))
      pdf(file = gfile, width = 12, height = 10, onefile = TRUE,
          title = paste(self$getLocus(), self$getSampleId(), sep = "." ))
      self$plotmapInitSummary(
        thin = 0.25,
        width = 2,
        label = self$getMapTag("init", "SR"),
        readtype = "SR"
      )
      dev.off()
    }
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
                                      plot = TRUE,
                                      ...) {
  flog.info("Step 1: Partition longreads into haplotypes ...", name = "info")
  Sys.sleep(1)
  x$runPartitionLongReads(threshold = threshold,
                             skip_gap_freq = skip_gap_freq,
                             dist_alleles = dist_alleles,
                             plot = plot)
  message("  Done!\n")
  invisible(x)

}
DR2S_$set("public", "runPartitionLongReads", function(threshold = NULL,
                                                         skip_gap_freq = 2/3,
                                                         dist_alleles = NULL,
                                                         plot = TRUE) {

  # debug
  # skip_gap_freq = 2/3
  # dist_alleles = 3
  # threshold = NULL
  # plot = TRUE
  # self <- dr2s


  ## Overide default arguments
  args <- self$getOpts("partition")
  if (!is.null(args)) {
    env  <- environment()
    list2env(args, envir = env)
  }
  if (is.null(threshold)){
    threshold <- self$getThreshold()
  }
  if (is.null(dist_alleles)){
    dist_alleles <- self$getDistAlleles()
  }
  assertthat::assert_that(self$hasPileup(),
                          is.double(skip_gap_freq),
                          is.double(threshold),
                          assertthat::is.count(dist_alleles),
                          is.logical(plot))

  tag  <- self$getMapTag("init", "LR")

  ## Get the reference sequence
  if (!is.null(self$mapInit$SR1)){
    useSR <- TRUE
    refseq <- self$mapInit$SR1$conseq
    flog.info(" Generating SNP matrix from shortreads ...", name = "info")
  } else {
    useSR <- FALSE
    refseq <- self$getRefSeq()
    flog.info(" Generating SNP matrix from longreads ...", name = "info")
  }
  ppos <- self$polymorphicPositions(useSR = useSR)

  ## Check if already finished because it is a homozygous sample
  if (NROW(ppos) == 0){
    flog.warn("No polymorphic positions for clustering! Only single allele?",
              name = "info")
    flog.info("Entering polish and report pipeline", name = "info")
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
  flog.info(" Partitioning %s reads over %s SNPs ...",
            NROW(mat), NCOL(mat), name = "info")
  prt <- partition_reads(x = mat,
                         skip_gap_freq = skip_gap_freq,
                         deepSplit = 1,
                         threshold = threshold,
                         dist_alleles = dist_alleles)
  ## Set sample haplotypes
  self$setHapTypes(levels(as.factor(PRT(prt))))

  # Check if we have only one cluster and finish the pipeline if so
  if (length(self$getHapTypes()) == 1){
    flog.warn("Only one allele left!")
    flog.info("Entering polish and report pipeline", name = "info")
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

## Method: splitReadsByHaplotype ####

#' @export
splitReadsByHaplotype.DR2S <- function(x,
                                          limits,
                                          ...) {
  flog.info(" Split partitioned reads by score ...", name = "info")
  x$runSplitReadsByHaplotype(limits)
  invisible(x)
}

#self <- dpb1_3
DR2S_$set("public", "runSplitReadsByHaplotype", function(limits,
                                                      plot = TRUE){

  ## Check if reporting is already finished and exit safely
  if (check_report_status(self)) return(invisible(self))
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
  self.setLimits <- sapply(haplotypes, function(x) NULL)

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
    flog.info("  %s: Using %s longreads with score > %s ...",
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
    tag   <- self$getMapTag("init", "LR")
    bnm   <- basename(self$mapInit$bamfile)
    outf  <- self$absPath(paste0("plot.partition.", sub("bam$", "pdf", usc(bnm))))
    pdf(file = outf, width = 10, height = 12, onefile = TRUE,
        title = paste(self$getLocus(), self$getSampleId(), sep = "." ))
    self$plotPartitionSummary(label = tag, limits = unlist(self$getLimits()))
    # debug
    # label = tag
    # limits = unlist(self$getLimits())
    dev.off()
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

## Method: extractLongReads ####

#' @export
extractLongReads.DR2S <- function(x,
                              nreads = NULL,
                              replace = FALSE,
                              nalign = 40,
                              ...) {
  flog.info("Step 2: Extracting haplotyped reads ...", name = "info")
  Sys.sleep(1)
  x$runExtractLongReads(nreads = nreads, replace = replace, nalign = nalign)
  message("  Done!\n")
  invisible(x)
}
# debug
#self <- dr2s
DR2S_$set("public", "runExtractLongReads", function(nreads = NULL,
                                             replace = FALSE,
                                             nalign = 40) {

  # debug
  # nreads = NULL
  # replace = FALSE
  # nalign = 40

  ## Check if reporting is already finished and exit safely
  if (check_report_status(self)) return(invisible(self))

  stopifnot(self$hasHapList())

  ## Overide default arguments
  args <- self$getOpts("extract")
  if (!is.null(args)) {
    env  <- environment()
    list2env(args, envir = env)
  }

  ## override nreads in config if not NULL
  if (!is.null(nreads))
    self$setNreads(nreads)

  ## do this for each haptype
  hptypes <- self$getHapTypes()
  for (hptype in hptypes){
    dir <- dir_create_if_not_exists(normalizePath(file.path(self$getOutdir(),
                                                            "mapIter",
                                                            (hptype)),
                                                  mustWork = FALSE))

    qnames <- self$getHapList(hptype)

    # x = self$mapInit$bamfile
    # qnames = qnames
    # n = self$getNreads()
    # replace = replace

    fq  <- .extractFastq(
      x = self$absPath(self$mapInit$bamfile),
      qnames = qnames,
      n = self$getNreads(),
      replace = replace
    )
    if (!is.null(nreads)) {
      attr(self$partition$hpl[hptype], "index") <-
        which(qnames %in% as.character(ShortRead::id(fq)))
      attr(self$partition$hpl[hptype], "n") <- self$getNreads()
    }
    file <- paste("hap", hptype, self$getLrdType(), self$getLrMapper(),
                  paste0("lim", 100 * abs(attr(self$getHapList(hptype),
                                               "limit"))),
                  paste0("n", length(fq)),
                  "fastq", "gz", sep = ".")
    out <- file_delete_if_exists(file.path(dir, file))
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
                         pruning_cutoff = 0.75,
                         force = FALSE,
                         fullname = TRUE,
                         plot = TRUE) {
  flog.info("Step 3: Iteratively mapping partitioned longreads against consensus sequences from initial mapping using only partitioned reads ... ",
            name = "info")
  Sys.sleep(1)
  x$runMapIter(opts = opts, iterations = 1, pct = pct,
               min_base_quality = min_base_quality,
               min_mapq = min_mapq, max_depth = max_depth,
               min_nucleotide_depth = min_nucleotide_depth,
               include_deletions = include_deletions,
               include_insertions = include_insertions,
               pruning_cutoff = pruning_cutoff, force = force,
               fullname = fullname, plot = plot)
  message("  Done!\n")
  invisible(x)
}

# self <- hla.fq
#self <- dpb1_3
DR2S_$set("public", "runMapIter", function(opts = list(),
                                           iterations = 1,
                                           pct = 100,
                                           min_base_quality = 3,
                                           min_mapq = 0,
                                           max_depth = 1e4,
                                           min_nucleotide_depth = 3,
                                           include_deletions = TRUE,
                                           include_insertions = TRUE,
                                           pruning_cutoff = 0.75,
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
  # pruning_cutoff = 0.75
  # force = FALSE
  # fullname = TRUE
  # plot = TRUE
  # iterations = 2
  # ##

  ## Check if reporting is already finished and exit safely
  if (check_report_status(self)) return(invisible(self))

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
  flog.info(" Constructing a consensus from initial mapping ...", name = "info")
  bamfile <- self$absPath(self$mapInit$bamfile)
  if (self$getPartSR()){
    ref <- self$mapInit$SR1$conseq
  } else {
    ref <- self$getRefSeq()
  }
  mat <- msa_from_bam(bamfile, ref, paddingLetter = ".")

  foreach(hptype = hptypes) %do% {
    flog.info(paste("  Constructing a consensus from initial mapping for",
                    "haplotype %s ...", " "), hptype, name = "info")
    readIds <- self$getHapList(hptype)
    cmat <- consmat(t(Biostrings::consensusMatrix(mat[readIds],
                                        as.prob = FALSE)[VALID_DNA(
                                          include = "indel"),]), freq = FALSE)

    conseq_name <- paste0("consensus.mapIter.0.", hptype)
    conseq <- conseq(x = cmat, name = conseq_name, type = "prob",
                     force_exclude_gaps = FALSE, prune_matrix = FALSE,
                     cutoff = pruning_cutoff)

    seqpath     <- self$absPath(file.path(self$mapIter[["0"]][[hptype]]$dir,
                             paste0(conseq_name, ".fa")))
    self$mapIter[["0"]][[hptype]]$ref <- "mapIter0"
    self$mapIter[["0"]][[hptype]]$conseq <- conseq
    self$mapIter[["0"]][[hptype]]$seqpath <- self$relPath(seqpath)
    self$mapIter[["0"]][[hptype]]$tag <- "mapIter0"
    Biostrings::writeXStringSet(
      Biostrings::DNAStringSet(gsub("[-+]", "", conseq)),
      seqpath
    )
    conseq
  }

  # iteration <- 1
  for (iteration in 1:iterations) {
    flog.info(" Iteration %s ", iteration, " of ", iterations, name = "info")

    iterationC = as.character(iteration)
    prevIteration <- self$mapIter[[as.character(as.numeric(iteration-1))]]

    # hptype = "A"
    foreach (hptype = hptypes) %do% {
      reftag   <- prevIteration[[hptype]]$ref
      outdir   <- self$absPath(prevIteration[[hptype]]$dir)
      readpath <- self$absPath(prevIteration[[hptype]]$reads)
      refpath  <- self$absPath(prevIteration[[hptype]]$seqpath)
      refseq  <- prevIteration[[hptype]]$conseq

      ## try collapsing homopolymers; Did not seem to improve the mapping
      # if (iteration %% 2 == 1)
      #   refpath <- .collapseHomopolymers(refpath,4)

      nreads <- if (!is.null(attr(self$getHapList(hptype), "n"))) {
        attr(self$getHapList(hptype), "n")
      } else {
        length(self$getHapList(hptype))
      }
      optsname <- sprintf("%s [%s]", hptype, nreads)
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
          params  = list(pruning_cutoff = pruning_cutoff),
          tag     = list()
        ),
        class = c("mapIter", "list")
      )

      # self$mapIter[["0"]][[hptype]]$ref <- "multialign"
      # self$mapIter[["0"]][[hptype]]$conseq <- cons
      # self$mapIter[["0"]][[hptype]]$seqpath <- consout
      # self$mapIter[["0"]][[hptype]]$tag <- "multialign"

      ## Run mapper
      flog.info("  Mapping of partitioned longreads of haplotype %s ...",
                hptype, name = "info")
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
      pileup$consmat <- .distributeGaps(pileup$consmat,
                                        bamfile,
                                        refseq,
                                        removeError = TRUE)
      self$mapIter[[iterationC]][[hptype]]$pileup = pileup

      # ## Construct consensus sequence
      flog.info("   Constructing a consensus ...", name = "info")
      conseq_name <- paste0("consensus.", sub(".sam.gz", "", basename(samfile)))
      conseq      <- conseq(x = pileup, name = conseq_name, type = "prob",
                            exclude_gaps = FALSE, prune_matrix = FALSE,
                            cutoff = pruning_cutoff)
      seqpath     <- file.path(outdir, paste0(conseq_name, ".fa"))
      self$mapIter[[iterationC]][[hptype]]$seqpath = self$relPath(seqpath)
      self$mapIter[[iterationC]][[hptype]]$conseq = conseq
      self$mapIter[[iterationC]][[hptype]]$ref = conseq_name
      Biostrings::writeXStringSet(
        Biostrings::DNAStringSet(gsub("[-+]", "", conseq)),
        seqpath
      )

      ## Set maptag
      self$mapIter[[iterationC]][[hptype]]$tag = maptag
    }

    if (plot) {
      flog.info(" Plotting MapIter Summary ...", name = "info")
      ## Coverage and base frequency
      gfile <- self$absPath(paste("plot.MapIter", iteration, self$getLrdType(), self$getLrMapper(),
              "pdf", sep = ".")
      )
      pdf(file = gfile, width = 8 * length(hptypes), height = 8, onefile = TRUE,
          title = paste(self$getLocus(), self$getSampleId(), sep = "." ))
      self$plotmapIterSummary(thin = 0.1, width = 4, iteration = iteration,
                              drop.indels = TRUE)
      ## map1: 0.1; 4
      dev.off()
      ## Consensus sequence probability
      gfile <- self$absPath(paste("plot.MapIter.conseq", iteration, self$getLrdType(),
              self$getLrMapper(), "pdf", sep = ".")
      )
      pdf(file = gfile, width = 16, height = 4*length(hptypes), onefile = TRUE,
          title = paste(self$getLocus(), self$getSampleId(), sep = "." ))
      self$plotmapIterSummaryConseqProb(text_size = 1.75,
                                        iteration = iteration,
                                        point_size = 0.75,
                                        threshold = "auto")
      dev.off()
    }
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

## Method: partition Shortreads ####

#' @export
partitionShortReads.DR2S <- function(x,
                                     opts = list(),
                                     force = FALSE) {
  flog.info("Step 4: Partition shortreads based on initial mapping and longread clustering ... ", name = "info")
  Sys.sleep(1)
  x$runPartitionShortReads(opts = opts,
                           force = force)
  message("  Done!\n")
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

  ## Check if reporting is already finished and exit safely
  if (check_report_status(self)) return(invisible(self))

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
  if (self$getPartSR()){
    flog.info(" Found shortread mapping from MapInit ...", name = "info")
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
    ref <- self$getRefSeq()
    refname <- names(ref)
    reffile <- self$getRefPath()
    ## Fetch mapper
    map_fun <- self$getSrMapFun()
    ## Run mapper
    flog.warn(" Found no shortread mapping from MapInit ...", name = "info")
    flog.info(" Mapping short reads against provided reference ...",
              name = "info")
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
  srpartition <- get_SR_partition_scores(refname, bamfile, mats,
                                         cores = "auto")

  ## Assign read to haplotype with highest probability,
  ## i.e. product over probabilities of each haplotype and choose max
  flog.info(" Get highest scoring haplotype for each read", name = "info")
  srpartition$haplotypes <- score_highest_SR(srpartition$srpartition,
                                             diffThreshold = 0.001)

  # Write fastqs
  foreach(hptype = hptypes ) %do% {
    srfilenames <- c()
    flog.info(" Write shortread fastq for haplotype %s ...", hptype,
              name = "info")
    fqs <- self$getShortreads()
    dontUseReads <- srpartition$haplotypes$read[
      !srpartition$haplotypes$read %in% dplyr::filter(
        srpartition$haplotypes, haplotype == hptype)$read]

    # write fastq's
    foreach(fq = fqs) %do% {

      srFastqHap = self$absPath(file.path(self$mapIter$`0`[[hptype]]$dir,
                             paste0(c(strsplit(basename(fq), "\\.")[[1]][1],
                                      hptype, "fastq.gz"), collapse = ".")))
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
  flog.info("Step 5: Mapping short and long reads against refined consensus sequences ...",
            name = "info")
  Sys.sleep(1)
  x$runMapFinal(opts = opts, pct = pct, min_base_quality = min_base_quality,
                min_mapq = min_mapq, max_depth = max_depth,
                min_nucleotide_depth = min_nucleotide_depth,
                include_deletions = include_deletions,
                include_insertions = include_insertions,
                fullname = fullname, plot = plot, clip = clip)
  message("  Done!\n")
  invisible(x)
}

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

  ## Check if reporting is already finished and exit safely
  if (check_report_status(self)) return(invisible(self))

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
  outdir       <- dir_create_if_not_exists(self$absPath(reftag))
  lastIter     <- self$mapIter[[max(names(self$mapIter))]]
  hptypes      <- self$getHapTypes()
  readpathsLR  <- sapply(hptypes, function(x) self$absPath(lastIter[[x]]$reads))
  names(readpathsLR) <- hptypes
  refpaths     <- sapply(hptypes, function(x) self$absPath(lastIter[[x]]$seqpath))
  refseqs      <- sapply(hptypes, function(x) lastIter[[x]]$conseq)
  names(refpaths) <- hptypes
  readpathsSR <- lapply(hptypes, function(x) self$absPath(unlist(self$srpartition[[x]]$SR)))
  names(readpathsSR) <- hptypes

  self$mapFinal = structure(
    list(
      dir     = self$relPath(outdir),
      sreads  = lapply(readpathsSR,self$relPath),
      lreads  = self$relPath(readpathsLR),
      ref     = self$relPath(refpaths),
      bamfile = list(),
      pileup  = list(),
      tag     = list(),
      seqpath = list()
    ), class = c("mapFinal", "list")
  )


  ## Remap long reads to the same reference sequences as short reads
  # debug
  # hptype = "B"
  for (hptype in hptypes) {
    flog.info(" Run mapFinal for haplotype %s ...", hptype, name = "info" )
    refpath  <- refpaths[[hptype]]
    refseq  <- refseqs[[hptype]]

    mapgroupLR <- paste0("LR", hptype)
    maptagLR   <- paste("mapFinal", mapgroupLR, self$getLrdType(),
                        self$getLrMapper(),
                        optstring(opts), sep = ".")
    readpathLR <- readpathsLR[[hptype]]
    ## Mapper
    map_fun <- self$getLrMapFun()
    ## Run mapper
    flog.info("  Mapping long reads against latest consensus ...",
              name = "info")
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
    cseq <- conseq(pileup$consmat, paste0("mapFinal", hptype), "ambig",
                   exclude_gaps = FALSE, threshold = 0.2)
    self$mapFinal$seq[[hptype]] <- cseq


    ## Map short reads
    if (!is.null(readpathsSR[[hptype]])){
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
      flog.info("  Mapping short reads against latest consensus ...",
                name = "info")

      samfile <- map_fun(
        reffile  = refpath,
        readfile = readfiles,
        # if we run shortreads against both pacbio and nanopore data
        # this hack makes sure that we can distinguish the bam files ->
        # we get pacbio.illumina.bwamem.A...bam and nanopore.illumina.bwamem.A...bam
        allele   = paste0(mapgroupSR, ".", self$getLrdType()),
        readtype = self$getSrdType(),
        opts     = opts,
        refname  = hptype,
        optsname = optstring(opts),
        force    = force,
        outdir   = outdir
      )

      if (clip) {
        flog.info(" Trimming softclips and polymorphic ends ...", name = "info")
        ## Run bam - sort - index pipeline
        bamfile <- bam_sort_index(samfile, refpath, pct / 100, min_mapq,
                                  force = force, clean = TRUE)
        ## Trim softclips
        fq <- trim_softclipped_ends(bam = Rsamtools::scanBam(bamfile)[[1]],
                                    preserve_ref_ends = TRUE)
        ## Trim polymorphic ends
        fq <- trim_polymorphic_ends(fq)
        ## Write new shortread file to disc
        fqdir  <- dir_create_if_not_exists(file.path(self$getOutdir(), "final"))
        fqfile <- paste("sread", hptype, self$getSrMapper(), "trimmed", "fastq",
                        "gz", sep = ".")
        fqout  <- file_delete_if_exists(file.path(fqdir, fqfile))
        ShortRead::writeFastq(fq, fqout, compress = TRUE)
        file_delete_if_exists(bamfile)
        ## Rerun mapper
        flog.info("  Mapping trimmed short reads against latest consensus ... ",
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
        file_delete_if_exists(fqout)
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
      cseq <- conseq(pileup$consmat, paste0("mapFinal", hptype), "ambig",
                     exclude_gaps = FALSE, threshold = 0.2)
      self$mapFinal$seq[[hptype]] <- cseq
    }

  }


  if (plot) {
    flog.info(" Plotting final summary figures ...", name = "info")
    ## Coverage and base frequency

    if (!is.null(self$mapFinal$sreads$A)){
      for (readtype in c("LR", "SR")){
        #readtype <- "LR"
        gfile <- self$absPath(paste("plot.mapFinal", readtype,
                                 self$getLrdType(), self$getLrMapper(),
                                 "pdf", sep = "."))
        pdf(file = gfile, width = 8 * length(hptypes), height = 8,
            onefile = TRUE,  title = paste(self$getLocus(),
                                           self$getSampleId(), sep = "." ))
        self$plotmapFinalSummary(iteration = "final", readtype = readtype,
                                 thin = 0.25, width = 20)
        dev.off()
      }
    } else {
      gfile <- self$absPath(paste("plot.mapFinal.LR",
                                                 self$getLrdType(),
                                                 self$getLrMapper(), "pdf",
                                                 sep = "."))
      pdf(file = gfile, width = 8 * length(hptypes), height = 8,
          onefile = TRUE, title = paste(self$getLocus(), self$getSampleId(),
                                        sep = "." ))
      self$plotmapFinalSummary(iteration = "final", readtype = "LR",
                               thin = 0.25, width = 20)
      dev.off()
    }
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
#'   flog.info("Step 6: Infer problematic positions and consensus sequence ...", name = "info")
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
#'   flog.info("Step 7: report consensus sequences and problematic positions ...", name = "info")
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
