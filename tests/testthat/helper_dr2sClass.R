testLogicalGetSet <- function(x, funValue) {
  # dr2s <- InitDR2S(createDR2SConf())
  #funValue <- "Microsatellite"
  getter <- paste0("get", funValue)
  setter <- paste0("set", funValue)
  getterFun <- dr2s[[getter]]
  setterFun <- dr2s[[setter]]
  expect_is(getterFun(), "logical")
  expect_error(setterFun("x"))
  expect_error(setterFun("1"))
  expect_error(setterFun(c()))
  setterFun(TRUE)
  expect_equal(getterFun(), TRUE)
  setterFun(FALSE)
  expect_equal(getterFun(), FALSE)
}

testCharacterGetSet <- function(x, funValue) {
  # dr2s <- InitDR2S(createDR2SConf())
  #funValue <- "Microsatellite"
  getter <- paste0("get", funValue)
  setter <- paste0("set", funValue)
  getterFun <- dr2s[[getter]]
  setterFun <- dr2s[[setter]]
  expect_is(getterFun(), "character")
  expect_error(setterFun(TRUE))
  expect_error(setterFun(1))
  setterFun("TestPath")
  expect_equal(getterFun(), "TestPath")
  setterFun("")
  expect_equal(getterFun(), "")
}