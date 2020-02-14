context("Initialization")

outdir <- tempfile()
params <- list(
  sample = "ID14984936",
  locus = "A",
  longreads = list(type = "pacbio", dir = "Sequel"),
  shortreads = list(type = "illumina", dir = "Illumina"),
  datadir = system.file("inst/testData/", package = "DR2S"),
  outdir = tempfile(),
  reference = "01:01:01:01")

conf <- createDR2SConf(
  sample         = params[["sample"]],
  locus          = params[["locus"]],
  longreads      = params[["longreads"]],
  shortreads     = params[["shortreads"]],
  datadir        = params[["datadir"]],
  outdir         = params[["outdir"]],
  reference      = params[["reference"]]
)

test_that("config is created from within R", {
  expect_is(conf, "DR2Sconf")
  expect_equal(conf$sampleId, params$sample)
  expect_equal(conf$locus, "HLA-A")
  expect_equal(conf$longreads, list(dir = "Sequel", type = "pacbio", mapper = "minimap"))
  expect_equal(conf$shortreads, list(dir = "Illumina", type = "illumina", mapper = "bwamem"))
  expect_equal(conf$datadir, params$datadir)
  expect_equal(conf$outdir, params$outdir)
  expect_equal(conf$reference, params$reference)
  expect_equal(conf$pipeline, "SR")
})

test_that("read from yaml", {
 TRUE
})

test_that("DR2S object is created", {
  initialiseDR2S(conf, createOutdir = TRUE)
})
