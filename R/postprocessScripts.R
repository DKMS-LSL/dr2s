## Make Rscript files for direct postprocessing
##
.writeReportCheckedConsensus <- function(path, libpath = ".pplib") {
  bashFile <- file.path(path, "run_reportCheckedConsensus.sh")
  rFile <- file.path(libpath, "reportCheckedConsensus.R")
  script <-
'library(DR2S)
## get the scripts dir for change wd
args <- commandArgs(trailingOnly = FALSE)
fileArgName <- "--file="
scriptName <- sub(fileArgName, "", args[grep(fileArgName, args)])
scriptBaseName <- dirname(scriptName)
setwd(scriptBaseName)

x <- readDR2S("..")
tryCatch(reportCheckedConsensus(x),
         error = function(e) {
           system(paste(
             shQuote("notify-send"),
             shQuote("-u"), shQuote("critical"),
             shQuote("Report failed! Ambiguous positions found"),
             shQuote("Run in terminal or look at the log to see whats wrong")))
         })
'

  write("#!/usr/bin/env bash\nRscript " %<<% rFile, bashFile)
  write(script, file.path(path, rFile))
  Sys.chmod(bashFile, mode = "775")
}

.writeCheckConsensus <- function(path, libpath = ".pplib") {
  bashFile <- file.path(path, "run_checkConsensus.sh")
  rFile <- file.path(libpath, "checkConsensus.R")
  script <-
'#!/usr/bin/env Rscript
library(DR2S)
## get the scripts dir for change wd
args <- commandArgs(trailingOnly = FALSE)
fileArgName <- "--file="
scriptName <- sub(fileArgName, "", args[grep(fileArgName, args)])
scriptBaseName <- dirname(scriptName)
setwd(scriptBaseName)

x <- readDR2S("..")
tryCatch(checkAlignmentFile(x),
         error = function(e) {
           system(paste(
             shQuote("notify-send"),
             shQuote("-u"), shQuote("critical"),
             shQuote(e),
             shQuote("Run in terminal to see whats wrong")))
         })
'
  write(script, file.path(path, rFile))
  write("#!/usr/bin/env bash\nRscript " %<<% rFile, bashFile)
  Sys.chmod(bashFile, mode = "775")
}

.writeRemapAlignments <- function(path, haptypes, libpath = ".pplib") {
  .writeScript <- function(hptype, path) {
    bashFile <- file.path(path, "run_remap" %<<% hptype %<<% ".sh")
    rFile    <- file.path(libpath, "remap" %<<% hptype %<<% ".R")
    script   <- sprintf(
'#!/usr/bin/env Rscript
library(DR2S)
## get the scripts dir for change wd
args <- commandArgs(trailingOnly = FALSE)
fileArgName <- "--file="
scriptName <- sub(fileArgName, "", args[grep(fileArgName, args)])
scriptBaseName <- dirname(scriptName)
setwd(scriptBaseName)

x <- readDR2S("..")
tryCatch({remapAlignment(x, "%s")
         },
         error = function(e) {
           system(paste(
             shQuote("notify-send"),
             shQuote("-u"), shQuote("critical"),
             shQuote(e),
             shQuote("Run in terminal to see whats wrong")))
         })
', hptype)
    write(script, file.path(path, rFile))
    write("#!/usr/bin/env bash\nRscript " %<<% rFile, bashFile)
    Sys.chmod(bashFile, mode = "775")
  }
  invisible(lapply(haptypes, function(hp) .writeScript(hp, path)))
}
