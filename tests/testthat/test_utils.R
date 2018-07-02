context("Utils")

test_that(".dirCreateIfNotExists creates dirs as expected", {
  ## Test input
  testDir <- 11
  expect_error(.dirCreateIfNotExists(testDir))
  ## Test dir creation  and vectorization
  for (dirName in list("test", c("test1", "test2"))) {
    testDir <- file.path(tempdir(), dirName)
    expect_false(any(dir.exists(testDir)))
    expect_equal(unname(.dirCreateIfNotExists(testDir)), testDir)
    expect_true(all(dir.exists(testDir)))
    expect_equal(unname(.dirCreateIfNotExists(testDir)), testDir)
    unlink(testDir, recursive = TRUE)
    expect_false(any(dir.exists(testDir)))
  }
})