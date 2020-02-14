#' @include zzz.R utils.R
NULL

normaliseOpts <- function(opts, pipeline = "LR") {
  if (has_attr(opts, "valid")) {
    return(opts)
  }
  pipeline <- match.arg(pipeline, c("SR", "LR"))
  ## camelCasify option names of passed-in options
  if (!is.null(opts))
    opts <- setNames(opts, .camelCasify(names(opts)))
  else
    opts <- list()
  ## set up defaults according to the pipeline type
  opts0 <- list()
  ##
  ## mapInit() defaults ####
  ##
  opts0$mapInit <- compact(list(
    ##
    ## <includeDeletions>: include deletions in pileup.
    includeDeletions = TRUE,
    ##
    ## <includeInsertions>: include insertions in pileup.
    includeInsertions = TRUE,
    ##
    ## <callInsertionThreshold>: if <includeInsertions == TRUE>, an insertion
    ## needs to be at frequency <callInsertionThreshold> for it to be included
    ## in the pileup.
    callInsertionThreshold = 0.18,
    ##
    ## <microsatellite>: if <pipeline == "SR">, perform a second mapping of
    ## shortreads to the inferred reference. Set to TRUE if you suspect
    ## microsatellites or repetitive regions in your sequence. This extends
    ## the reference to a maximum length and enables a better mapping.
    microsatellite = FALSE,
    ##
    ## <forceMapping>: set to TRUE if you want to force processing of "bad"
    ## shortreads where the distribution of coverage is heavily unequal.
    ## Aborts the program if maximum coverage > 75 % quantile * 5.
    forceMapping = FALSE,
    ##
    ## <minMapq>: don't filter longreads for mapping quality unless specified.
    ## NOTE: for shortreads <minMapq = 50> is hardcoded.
    minMapq = 0,
    ##
    ## <topx>: pick the x top-scoring reads. Set to an integer value to pick
    ## a fixed number of reads. Set to "auto" to use a dynamically determined
    ## number of reads to be selected.
    topx = FALSE,
    ##
    ## <pickiness>: if <topx == "auto">: <pickiness < 1>: bias towards higher
    ## scores/less reads; <pickiness > 1>: bias towards lower scores/more reads
    pickiness = 0.8,
    ##
    ## <increasePickiness>: if <topx == "auto">: increase pickiness for the
    ## second iteration of the LR mapping
    increasePickiness = 1,
    ##
    ## <lowerLimit>: if <topx == "auto"> or <topx > 0>: the  minimum number
    ## of reads to pick if available.
    lowerLimit = 200,
    ##
    ## <updateBackgroundModel>: estimate the indel noise in a pileup and use
    ## this information to update the background model for PWM scoring
    updateBackgroundModel = TRUE,
    ##
    ## <createIgv>: subsample bam files for visualisation with IgvJs in the
    ## DR2S shiny app.
    createIgv = TRUE,
    ##
    ## <plot>: generate diagnostic plots.
    plot = TRUE
  ))
  ##
  ## partitionLongreads() defaults ####
  ##
  opts0$partitionLongreads <- compact(list(
    ##
    ## <threshold>: the threshold to call a polymorphic position. A minority
    ## nucleotide frequency below this threshold is considered noise rather
    ## than a valid polymorphism.
    threshold = 0.2,
    ##
    ## <distAlleles>: the expected number of distinct alleles in the sample.
    ## This should be 2 for heterozygous samples, 1 for homozygous samples,
    ## and may be set to >2 for some KIR loci.
    distAlleles = 2,
    ##
    ## <skipGapFreq>: the minumum frequency of the gap character required to
    ## call a gap position.
    skipGapFreq = 0.8,
    ##
    ## <noGapPartitioning>: don't partition based on gaps. Useful for samples
    ## with only few SNPs but with homopolymers. The falsely called gaps could
    ## mask the real variation. Set to override global default.
    noGapPartitioning = TRUE,
    ##
    ## <selectCorrelatedPositions>: correlate polymorphic positions and cluster
    ## based on the absolute correlation coefficient. Extract positions from the
    ## cluster with the higher absolute mean correlation coefficient. This gets
    ## rid of positions that are not well distributed across the two alleles.
    selectCorrelatedPositions = FALSE,
    ##
    ## <measureOfAssociation>: if <selectCorrelatedPositions> == TRUE, use
    ## <measureOfAssociation> ("cramer.V" or "spearman") to determine linkage
    ## between all polymorphic positions.
    measureOfAssociation = "cramer.V",
    ##
    ## <selectByColSum>: use <selectByColSum> to decide which SNP cluster to
    ## use based on overal sums. Otherwise the intra-cluster values are used
    selectByColSum = FALSE,
    ##
    ## <proportionOfOverlap>: we perform an equivalence test on clusters of
    ## polymorhic positions: Calculate the lower 1-sigma bound of the
    ## high-association cluster i. Calculate the upper 1-sigma bound of the
    ## low-association cluster j. Reject the clusters, if this bounds overlap
    ## by more than <proportionOfOverlap> of the average distance (dij) between
    ## clusters.
    proportionOfOverlap = 0.1,
    ##
    ## <minimumExpectedDifference>: by how much do we expect two clusters to
    ## minimally differ in mean Cramér's V. BIC-informed model-based clustering
    ## tends to split rather than lump and this is a heuristical attempt to
    ## forestall this.
    minimumExpectedDifference = 0.06,
    ##
    ## If more than <distAlleles> clusters are found select clusters based on:
    ## (1) "distance": The hamming distance of the resulting variant consensus
    ## sequences or (2) "count": Take the clusters with the most reads as the
    ## true alleles.
    selectAllelesBy = "distance",
    ##
    ## Minimum size of an allele cluster
    minClusterSize = 20,
    ##
    ## When selecting reads from allele clusters using a dynamic threshold:
    ## pickiness < 1: bias towards higher scores/less reads
    ## pickiness > 1: bias towards lower scores/more reads
    pickiness = 0.2,
    ##
    ## When selecting reads from allele clusters the minimum number of
    ## reads to pick if available.
    lowerLimit = 40,
    ##
    ## Generate diagnostic plots.
    plot = TRUE
  ))
  ##
  ## mapIter() defaults ####
  ##
  opts0$mapIter <- compact(list(
    ##
    ## Number of <mapIter> iterations. How often are the
    ## clustered reads remapped to updated reference sequences.
    iterations = 2,
    ##
    ## Minimum occupancy (1 - fraction of gap) below which
    ## bases at insertion position are excluded from from consensus calling.
    columnOccupancy = 0.2,
    ##
    ## an insertion needs to be at frequency <callInsertionThreshold> for it
    ## to be included in the pileup.
    callInsertionThreshold = 0.18,
    ##
    ## <createIgv>: subsample bam files for visualisation with IgvJs in the
    ## DR2S shiny app.
    createIgv = TRUE,
    ##
    ## Generate diagnostic plots.
    plot = TRUE
  ))
  ##
  ## partitionShortreads() defaults ####
  ##
  if (pipeline == "SR") {
    opts0$partitionShortreads <- compact(list(

    ))
  }
  ##
  ## mapFinal() defaults ####
  ##
  opts0$mapFinal <- compact(list(
    ##
    ## include deletions in pileup.
    includeDeletions = TRUE,
    ##
    ## include insertions in pileup.
    includeInsertions = TRUE,
    ##
    ## an insertion needs to be at frequency <callInsertionThreshold> for it
    ## to be included in the pileup.
    callInsertionThreshold = 0.18,
    ##
    ## (for shortreads only) trim softclips and polymorphic ends of reads before
    ## the final mapping
    trimPolymorphicEnds = FALSE,
    ##
    ## Subsample bam files for visualisation with IgvJs in the
    ## DR2S shiny app.
    createIgv = TRUE,
    ##
    ## Generate diagnostic plots.
    plot = TRUE
  ))
  ##
  ## report() defaults ####
  ##
  opts0$report <- compact(list(
    ##
    ## Suppress remapping of reads against final consensus.
    remap = TRUE,
    ##
    ## Threshold to call a polymorphic position. Set to override global default.
    threshold = NULL,
    ##
    ## Check the number of homopolymer counts. Compare the resulting sequence
    ## with the mode value and report differences.
    checkHpCount = TRUE,
    ##
    ## The minimal length of a homopolymer to be checked.
    hpCount = 10,
    ##
    ## Maximum number of sequence letters per line in pairwise alignment.
    blockWidth = 80,
    ##
    ## Subsample bam files for visualisation with IgvJs in the
    ## DR2S shiny app.
    createIgv = TRUE
  ))
  opts0 <- compact(opts0)
  ## update default with config settings
  opts1 <- utils::modifyList(compact(opts0), opts, keep.null = FALSE)
  validateOpts(opts = opts1)
}

MANDATORY_OPTS <- function() {
  c("mapInit", "partitionLongreads", "mapIter", "mapFinal", "report")
}

validateOpts <- function(opts) {
  fields <- MANDATORY_OPTS()
  assert_that(all(fields %in% names(opts)),
              msg = paste0("Missing opts for <",
                           comma(fields[!fields %in% names(opts)]),
                           "> in config"))
  ##
  ## mapInit() asserts ####
  ##
  assert_that(
    is.flag(opts$mapInit$includeDeletions),
    is.flag(opts$mapInit$includeInsertions),
    is.flag(opts$mapInit$microsatellite),
    is.flag(opts$mapInit$forceMapping),
    is.flag(opts$mapInit$updateBackgroundModel),
    is.number(opts$mapInit$minMapq),
    is.number(opts$mapInit$pickiness),
    is.number(opts$mapInit$increasePickiness),
    is.number(opts$mapInit$lowerLimit),
    is.flag(opts$mapInit$createIgv),
    is.flag(opts$mapInit$plot)
  )
  assert_that(
    is.number(opts$mapInit$callInsertionThreshold),
    opts$mapInit$callInsertionThreshold >= 0,
    opts$mapInit$callInsertionThreshold <= 1,
    msg = "<callInsertionThreshold> in mapInit() is not a number between 0 and 1")
  assert_that(
    opts$mapInit$topx == "auto" || is.number(opts$mapInit$topx) || is.flag(opts$mapInit$topx),
    msg = "<topx> in mapInit() is not \"auto\" or a number >= 0")
  ##
  ## partitionLonreads() asserts ####
  ##
  assert_that(
    is.flag(opts$partitionLongreads$noGapPartitioning),
    is.flag(opts$partitionLongreads$selectCorrelatedPositions),
    is.number(opts$partitionLongreads$minClusterSize),
    is.number(opts$partitionLongreads$pickiness),
    is.number(opts$partitionLongreads$lowerLimit),
    is.flag(opts$partitionLongreads$plot)
  )
  assert_that(
    is.number(opts$partitionLongreads$threshold),
    opts$partitionLongreads$threshold >= 0,
    opts$partitionLongreads$threshold <= 1,
    msg = "<threshold> in partitionLongreads() is not a number between 0 and 1")
  assert_that(
    is.count(opts$partitionLongreads$distAlleles),
    opts$partitionLongreads$distAlleles > 0,
    opts$partitionLongreads$distAlleles < 5,
    msg = "<distAlleles> in partitionLongreads() is not a count between 1 and 4")
  assert_that(
    is.number(opts$partitionLongreads$skipGapFreq),
    opts$partitionLongreads$skipGapFreq >= 0,
    opts$partitionLongreads$skipGapFreq <= 1,
    msg = "<skipGapFreq> in partitionLongreads() is not a number between 0 and 1")
  assert_that(
    is.string(opts$partitionLongreads$selectAllelesBy),
    opts$partitionLongreads$selectAllelesBy %in% c("count", "distance"),
    msg = "<selectAllelesBy> in partitionLongreads() is not 'count' nor 'distance'"
  )
  assert_that(
    is.string(opts$partitionLongreads$measureOfAssociation),
    opts$partitionLongreads$measureOfAssociation %in% c("cramer.V", "spearman"),
    msg = "<measureOfAssociation> in partitionLongreads() is not 'cramer.V' nor 'spearman'"
  )
  assert_that(
    is.logical(opts$partitionLongreads$selectByColSum),
    msg = "<selectByColSum> in partitionLongreads() is not logical"
  )
  assert_that(
    is.number(opts$partitionLongreads$proportionOfOverlap),
    opts$partitionLongreads$proportionOfOverlap >= 0,
    opts$partitionLongreads$proportionOfOverlap <= 1,
    msg = "<proportionOfOverlap> in partitionLongreads() is not a number between 0 and 1")
  assert_that(
    is.number(opts$partitionLongreads$minimumExpectedDifference),
    opts$partitionLongreads$minimumExpectedDifference >= 0,
    opts$partitionLongreads$minimumExpectedDifference <= 1,
    msg = "<minimumExpectedDifference> in partitionLongreads() is not a number between 0 and 1")
  ##
  ## mapIter() asserts ####
  ##
  assert_that(
    is.flag(opts$mapIter$plot)
  )
  assert_that(
    is.count(opts$mapIter$iterations),
    opts$mapIter$iterations > 0,
    opts$mapIter$iterations < 10,
    msg = "<iterations> in mapIter() is not a count between 1 and 9"
  )
  assert_that(
    is.number(opts$mapIter$columnOccupancy),
    opts$mapIter$columnOccupancy >= 0,
    opts$mapIter$columnOccupancy <= 1,
    msg = "<columnOccupancy> in mapIter() is not a number between 0 and 1")
  assert_that(
    is.number(opts$mapIter$callInsertionThreshold),
    opts$mapIter$callInsertionThreshold >= 0,
    opts$mapIter$callInsertionThreshold <= 1,
    msg = "<callInsertionThreshold> in mapIter() is not a number between 0 and 1")
  ##
  ## mapFinal() asserts ####
  ##
  assert_that(
    is.flag(opts$mapFinal$includeDeletions),
    is.flag(opts$mapFinal$includeInsertions),
    is.flag(opts$mapFinal$trimPolymorphicEnds),
    is.flag(opts$mapFinal$createIgv),
    is.flag(opts$mapFinal$plot)
  )
  assert_that(
    is.number(opts$mapFinal$callInsertionThreshold),
    opts$mapFinal$callInsertionThreshold >= 0,
    opts$mapFinal$callInsertionThreshold <= 1,
    msg = "<callInsertionThreshold> in mapFinal() is not a number between 0 and 1")
  ##
  ## report() asserts ####
  ##
  assert_that(
    is.count(opts$report$blockWidth),
    is.flag(opts$report$remap),
    is.flag(opts$report$createIgv),
    is.flag(opts$report$checkHpCount),
    is.count(opts$report$hpCount)
  )
  assert_that(
    is.null(opts$report$threshold) || (
      is.number(opts$threshold$threshold) && opts$report$threshold >= 0 && opts$report$threshold <= 1
    ),
    msg = "<threshold> in report()) is not NULL or a number between 0 and 1")

  attr(opts, "valid") <- TRUE
  return(opts)
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

.normaliseUtrLength <- function(utrLength) {
  loci <- names(utrLength)
  loci <- vapply(loci, .normaliseLocus, FUN.VALUE = character(1))
  names(utrLength) <- loci
  utrLength
}
