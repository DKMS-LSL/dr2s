context("Initialization")

outdir <- tempfile()
params <- list(
  sample = "ID14984936",
  locus = "A",
  longreads = list(type = "pacbio", dir = "Sequel"),
  shortreads = list(type = "illumina", dir = "Illumina"),
  datadir = system.file("inst/testData/", package = "DR2S"),
  outdir = tempfile(),
  reference = "01:01:01:01",
  srmapper = "bwamem",
  lrmapper = "minimap",
  iterations = 2,
  microsatellite = TRUE,
  partSR = TRUE,
  distAlleles = 2,
  filterScores = FALSE,
  threshold = 0.10,
  forceMapping = FALSE)

conf <- createDR2SConf(
  sample         = params[["sample"]],
  locus          = params[["locus"]],
  longreads      = params[["longreads"]],
  shortreads     = params[["shortreads"]],
  datadir        = params[["datadir"]],
  outdir         = params[["outdir"]],
  reference      = params[["reference"]],
  srmapper       = params[["srmapper"]],
  lrmapper       = params[["lrmapper"]],
  iterations     = params[["iterations"]],
  microsatellite = params[["microsatellite"]],
  partSR         = params[["partSR"]],
  distAlleles    = params[["distAlleles"]],
  filterScores   = params[["filterScores"]],
  threshold      = params[["threshold"]],
  forceMapping   = params[["forceMapping"]]
)

test_that("config is created from within R", {
  expect_is(conf, "DR2Sconf")
  expect_equal(conf$sample, params$sample)
  expect_equal(conf$locus, params$locus)
  expect_equal(conf$longreads, params$longreads)
  expect_equal(conf$shortreads, params$shortreads)
  expect_equal(conf$datadir, params$datadir)
  expect_equal(conf$outdir, params$outdir)
  expect_equal(conf$reference, params$reference)
  expect_equal(conf$srmapper, params$srmapper)
  expect_equal(conf$lrmapper, params$lrmapper)
  expect_equal(conf$iterations, params$iterations)
  expect_equal(conf$microsatellite, params$microsatellite)
  expect_equal(conf$partSR, params$partSR)
  expect_equal(conf$distAlleles, params$distAlleles)
  expect_equal(conf$filterScores, params$filterScores)
  expect_equal(conf$threshold, params$threshold)
  expect_equal(conf$forceMapping, params$forceMapping)
})
test_that("read from yaml", {
 TRUE 
}) 

test_that("validateDR2Sconf works correct", {
  expect_equal(conf, validateDR2SConf(conf))
  ## Check sampleID
  eConf <- conf
  eConf$sampleId <- NULL
  expect_error(validateDR2SConf(eConf))
  ## Check locus
  eConf <- conf
  eConf$locus <- ""
  expect_error(validateDR2SConf(eConf))
  eConf$locus <- 123
  expect_error(validateDR2SConf(eConf))
  eConf$locus <- "HLA-Z"
  expect_error(validateDR2SConf(eConf))
  eConf$locus <- "HLA-A"
  expect_equal(validateDR2SConf(eConf)$locus, "A")
  eConf$locus <- "3DL3"
  expect_equal(validateDR2SConf(eConf)$locus, "KIR3DL3")
  eConf$locus <- "KIR3DL8"
  expect_error(validateDR2SConf(eConf))
  ## Check longreads
  eConf <- conf
  expect_equal(validateDR2SConf(eConf)$longreads, 
               list(type = "pacbio", dir = "Sequel"))
  eConf$longreads <- "pacbio"
  expect_error(validateDR2SConf(eConf))
  eConf$longreads <- list(type = "illumina", dir = "Sequel")
  expect_error(validateDR2SConf(eConf))
  eConf$longreads <- list(type = "pacbio", dir = "Sequel1")
  expect_error(validateDR2SConf(eConf))
  eConf$longreads <- NULL
  expect_error(validateDR2SConf(eConf))
  ## Check shortreads
  eConf <- conf
  expect_equal(validateDR2SConf(eConf)$shortreads, 
               list(type = "illumina", dir = "Illumina"))
  conf0 <- createDR2SConf(
    sample         = params[["sample"]],
    locus          = params[["locus"]],
    longreads      = params[["longreads"]],
    datadir        = params[["datadir"]],
    outdir         = params[["outdir"]],
    reference      = params[["reference"]],
    srmapper       = params[["srmapper"]],
    lrmapper       = params[["lrmapper"]],
    iterations     = params[["iterations"]],
    microsatellite = params[["microsatellite"]],
    partSR         = params[["partSR"]],
    distAlleles    = params[["distAlleles"]],
    filterScores   = params[["filterScores"]],
    threshold      = params[["threshold"]],
    forceMapping   = params[["forceMapping"]]
  )
  expect_null(validateDR2SConf(conf0)$shortreads)
  eConf$shortreads <- ""
  expect_error(validateDR2SConf(eConf))
  eConf$shortreads <- list(type = "pacbio", dir = "Sequel")
  expect_error(validateDR2SConf(eConf))
  eConf$longreads <- list(type = "illumina", dir = "Illumina1")
  expect_error(validateDR2SConf(eConf))
  ## Test logicals
  eConf <- conf
  eConf$microsatellite <- ""
  expect_error(validateDR2SConf(eConf))
  eConf <- conf
  eConf$filterScores <- ""
  expect_error(validateDR2SConf(eConf))
  eConf <- conf
  eConf$partSR <- ""
  expect_error(validateDR2SConf(eConf))
  eConf <- conf
  eConf$forceMapping <- ""
  expect_error(validateDR2SConf(eConf))
  ## Test counts
  eConf <- conf
  eConf$threshold <- ""
  expect_error(validateDR2SConf(eConf))
  eConf <- conf
  eConf$iterations <- ""
  expect_error(validateDR2SConf(eConf))
  eConf <- conf
  eConf$distAlleles<- ""
  expect_error(validateDR2SConf(eConf))
  eConf$threshold
})

test_that("DR2S object is created", {
  initialiseDR2S(conf, createOutdir = TRUE)
})
