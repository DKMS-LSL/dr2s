context("Initialization")

outdir <- tempfile()
params <- list(
  sample = "sample3",
  locus = "A",
  longreads = list(type = "nanopore", dir = "nanopore_sampled", mapper = "minimap"),
  shortreads = list(type = "illumina", dir = "Illumina_sampled", mapper = "rsubread"),
  datadir = system.file("testData/", package = "DR2S"),
  outdir = tempfile(),
  reference = "HLA-A*01:01:01:01")
config <- createDR2SConf(
  sample         = params[["sample"]],
  locus          = params[["locus"]],
  longreads      = params[["longreads"]],
  shortreads     = params[["shortreads"]],
  datadir        = params[["datadir"]],
  outdir         = params[["outdir"]],
  reference      = params[["reference"]]
)
test_that("config is created from within R", {
  expect_is(config, "DR2Sconf")
  expect_equal(config$sampleId, params$sample)
  expect_equal(config$locus, "HLA-A")
  expect_equal(config$longreads, list(dir = "nanopore_sampled", type = "nanopore", mapper = "minimap"))
  expect_equal(config$shortreads, list(dir = "Illumina_sampled", type = "illumina", mapper = "rsubread"))
  expect_equal(
    normalizePath(config$datadir, mustWork = TRUE), 
    normalizePath(params$datadir, mustWork = TRUE))
  expect_equal(
    normalizePath(config$outdir, mustWork = FALSE), 
    normalizePath(params$outdir, mustWork = FALSE))
  expect_equal(config$reference, params$reference)
  expect_equal(config$pipeline, "SR")
})

## Init the config
aInit <- InitDR2S(config, createOutdir = TRUE)
test_that("DR2S object is created and of correct type", {
  expect_s3_class(aInit, "DR2S")
})

## Write the config to json
confFile <- paste0(tempfile(), ".json")
test_that("write config to json", {
  DR2S::writeDR2SConf(aInit, confFile, format = "json")
  expect_true(file.exists(confFile))
})

test_that("read config from json", {
  configYaml <- readDR2SConf(confFile)
  expect_equal(configYaml, config)
})

