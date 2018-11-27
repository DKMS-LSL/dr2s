#' @include zzz.R utils.R
NULL

normaliseOpts <- function(opts, pipeline = "LR") {
  pipeline <- match.arg(pipeline, c("SR", "LR"))
  ## camelCasify option names of passed-in options
  opts <- lapply(opts, function(o) setNames(o, .camelCasify(names(o))))
  ## set up defaults according to the pipeline type
  opts0 <- list()
  ##
  ## set mapInit defaults
  ##
  opts0$mapInit <- compact(list(
    ## include deletions in pileup.
    includeDeletions = TRUE,
    ## include insertions in pileup.
    includeInsertions = TRUE,
    ## an insertion needs to be at frequency <callInsertionThreshold> for it
    ## to be included in the pileup.
    callInsertionThreshold = 1/5,
    ## don't filter for mapping quality unless specified.
    ## for <shortreads> we hardcode <minMapq = 50>
    minMapq = 0,
    ## Pick x top-scoring reads.
    topx = 0,
    ## if picking x top-scoring reads using a dynamic threshold:
    ## pickiness < 1: bias towards higher scores/less reads
    ## pickiness > 1: bias towards lower scores/more reads
    pickiness = 1,
    ## increase pickiness for the second iteration of LR mapping
    increasePickiness = 1,
    ## when picking x top-scoring reads the  minimum number of
    ## reads to pick if available.
    lowerLimit = 200,
    ## estimate the indel noise in a pileup and use this information to
    ## update the background model for PWM scoring
    updateBackgroundModel = FALSE,
    ## Subsample bam files for visualisation with IgvJs in the
    ## DR2S shiny app.
    createIgv = TRUE
  ))
  ##
  ## set partitionLongreads defaults
  ##
  opts0$partitionLongreads <- compact(list(
    ## Threshold to call a polymorphic position. Set to override global default.
    threshold = NULL,
    ## Don't partition based on gaps. Useful for samples with only few SNPs but
    ## with homopolymers. The falsely called gaps could mask the real variation.
    ## Set to override global default.
    distAlleles = NULL,
    ## The minumum frequency of the gap character required to call a gap position.
    skipGapFreq = 2/3,
    ## The number of distinct alleles in the sample. Set to override global default.
    noGapPartitioning = TRUE,
    ## Correlate polymorphic positions and cluster based on the absolute
    ## correlation coefficient. Extract positions from the cluster with the
    ## higher absolute mean correlation coefficient. This gets rid of positions
    ## that are not well distributed across the two alleles.
    restrictToCorrelatedPositions = FALSE,
    ## If more than <distAlleles> clusters are found select clusters based on:
    ## (1) "distance": The hamming distance of the resulting variant consensus
    ## sequences or (2) "count": Take the clusters with the most reads as the
    ## true alleles.
    selectAllelesBy = "distance",
    ## Minimum size of an allele cluster
    minClusterSize = 20,
    ## When selecting reads from allele clusters using a dynamic threshold:
    ## pickiness < 1: bias towards higher scores/less reads
    ## pickiness > 1: bias towards lower scores/more reads
    pickiness = 1,
    ## When selecting reads from allele clusters the minimum number of
    ## reads to pick if available.
    lowerLimit = 40
  ))
  ##
  ## set mapIter defaults
  ##
  opts0$mapIter <- compact(list(
    ## Number of <mapIter> iterations. How often are the
    ## clustered reads remapped to updated reference sequences.
    ## Set to override global default.
    iterations = NULL,
    ## Minimum occupancy (1 - fraction of gap) below which
    ## bases at insertion position are excluded from from consensus calling.
    columnOccupancy = 2/5,
    ## an insertion needs to be at frequency <callInsertionThreshold> for it
    ## to be included in the pileup.
    callInsertionThreshold = 1/5
  ))
  ##
  ## set partitionShortreads defaults
  ##
  if (pipeline == "SR") {
    opts0$partitionShortreads <- compact(list(

    ))
  }
  ##
  ## set mapFinal defaults
  ##
  opts0$mapFinal <- compact(list(
    ## include deletions in pileup.
    includeDeletions = TRUE,
    ## include insertions in pileup.
    includeInsertions = TRUE,
    ## an insertion needs to be at frequency <callInsertionThreshold> for it
    ## to be included in the pileup.
    callInsertionThreshold = 1/5,
    ## (for shortreads only) trim softclips and polymorphic ends of reads before
    ## the final mapping
    trimPolymorphicEnds = FALSE,
    ## Subsample bam files for visualisation with IgvJs in the
    ## DR2S shiny app.
    createIgv = TRUE
  ))
  ##
  ## set polish defaults
  ##
  opts0$polish <- compact(list(
    ## Threshold to call a polymorphic position. Set to override global default.
    threshold = NULL,
    ## Check the number of homopolymer counts. Compare the resulting sequence
    ## with the mode value and report differences.
    checkHpCount = TRUE,
    ## The minimal length of a homopolymer to be checked.
    hpCount = 10
  ))
  ##
  ## set report defaults
  ##
  opts0$report <- compact(list(
    ## Maximum number of sequence letters per line in pairwise alignment.
    blockWidth = 80,
    ## Suppress remapping of reads against final consensus.
    noRemap = FALSE,
    ## Subsample bam files for visualisation with IgvJs in the
    ## DR2S shiny app.
    createIgv = TRUE
  ))
  opts0 <- compact(opts0)
  ## update default with config settings
  opts1 <- utils::modifyList(compact(opts0), opts, keep.null = FALSE)
  opts1
}

.capitalise <- function(x) {
  capped <- grep("^[A-Z]", x, invert = TRUE)
  substr(x[capped], 1, 1) <- toupper(substr(x[capped], 1, 1))
  x
}

.camelCasify <- function(x, sep = "_") {
  xs <- strsplit(x, sep, fixed = TRUE)
  vapply(xs, function(x) {
    paste0(x[1L], paste0(.capitalise(x[-1L]), collapse = ""))
  }, FUN.VALUE = character(1))
}





