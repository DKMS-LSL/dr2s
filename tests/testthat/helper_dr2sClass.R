# dr2sRds <- system.file(package = "DR2S", "test/DR2S.pacbio.minimap.rds")
# dr2s <- readRDS(dr2sRds)
# testLogicalGetSet <- function(x, funValue) {
#   # dr2s <- InitDR2S(createDR2SConf())
#   #funValue <- "Microsatellite"
#   getter <- "get" %<<% funValue
#   setter <- "set" %<<% funValue
#   getterFun <- dr2s[[getter]]
#   setterFun <- dr2s[[setter]]
#   expect_is(getterFun(), "logical")
#   expect_error(setterFun("x"))
#   expect_error(setterFun("1"))
#   expect_error(setterFun(c()))
#   setterFun(TRUE)
#   expect_equal(getterFun(), TRUE)
#   setterFun(FALSE)
#   expect_equal(getterFun(), FALSE)
# }
# testLogicalGet <- function(x, funValue) {
#   # dr2s <- InitDR2S(createDR2SConf())
#   #funValue <- "Microsatellite"
#   getter <- "get" %<<% funValue
#   getterFun <- dr2s[[getter]]
#   expect_is(getterFun(), "logical")
# }
# 
# testCharacterGetSet <- function(x, funValue) {
#   # dr2s <- InitDR2S(createDR2SConf())
#   #funValue <- "Microsatellite"
#   getter <- "get" %<<% funValue
#   setter <- "set" %<<% funValue
#   getterFun <- dr2s[[getter]]
#   setterFun <- dr2s[[setter]]
#   expect_is(getterFun(), "character")
#   expect_error(setterFun(TRUE))
#   expect_error(setterFun(1))
#   setterFun("TestPath")
#   expect_equal(getterFun(), "TestPath")
#   setterFun("")
#   expect_equal(getterFun(), "")
# }
