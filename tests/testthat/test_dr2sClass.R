context("DR2S class")

test_that("microsatellite set and get works", {
  dr2s <- readRDS("dr2s.rds")
  testLogicalGetSet(dr2s, "Microsatellite")
  # dr2s <- InitDR2S(createDR2SConf())
})

test_that("partSR set and get works", {
  dr2s <- readRDS("dr2s.rds")
  # dr2s <- InitDR2S(createDR2SConf())
  testLogicalGetSet(dr2s, "PartSR")
})

test_that("forceMapping set and get works", {
  dr2s <- readRDS("dr2s.rds")
  # dr2s <- InitDR2S(createDR2SConf())
  testLogicalGetSet(dr2s, "ForceMapping")
})

test_that("filterScores set and get works", {
  dr2s <- readRDS("dr2s.rds")
  # dr2s <- InitDR2S(createDR2SConf())
  testLogicalGetSet(dr2s, "FilterScores")
})

test_that("getShortreads returns character", {
  dr2s <- readRDS("dr2s.rds")
  # dr2s <- InitDR2S(createDR2SConf())
  expect_true(is.character(dr2s$getShortreads()))
})

test_that("absPath gets correct paths", {
  dr2s <- readRDS("dr2s.rds")
  outdir <- dr2s$getOutdir()
  filename <- "a"
  expected <- setNames(file.path(outdir, filename), filename)
  expect_equal(dr2s$absPath(filename), expected)
  filename <- c("a", "b")
  expected <- setNames(file.path(outdir, filename), filename)
  expect_equal(dr2s$absPath(filename)[1], expected[1])
  expect_equal(dr2s$absPath(filename)[2], expected[2])
  filename <- 1
  expect_error(dr2s$absPath(filename))
  filename <- c(1,2)
  expect_error(dr2s$absPath(filename))
})



