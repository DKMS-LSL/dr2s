## Make Rscript files for direct postprocessing
##
writeReportCheckedConsensus <- function(path, libpath = ".pplib") {
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

x <- read_dr2s("..")
tryCatch(report_checked_consensus(x),
         error = function(e) {
           system(paste(
             shQuote("notify-send"),
             shQuote("-u"), shQuote("critical"),
             shQuote("Report failed! Ambiguous positions found"),
             shQuote("Run in terminal or look at the log to see whats wrong")))
         })
'

  write(paste0("#!/usr/bin/env bash\nRscript ", rFile), bashFile)
  write(script, file.path(path, rFile))
  Sys.chmod(bashFile, mode = "775")
}

writeCheckConsensus <- function(path, libpath = ".pplib") {
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

x <- read_dr2s("..")
tryCatch(check_alignment_file(x),
         error = function(e) {
           system(paste(
             shQuote("notify-send"),
             shQuote("-u"), shQuote("critical"),
             shQuote(e),
             shQuote("Run in terminal to see whats wrong")))
         })
'
  write(script, file.path(path, rFile))
  write(paste0("#!/usr/bin/env bash\nRscript ", rFile), bashFile)
  Sys.chmod(bashFile, mode = "775")
}

writeRefineAlignments <- function(path, haptypes, libpath = ".pplib") {
  writeScript <- function(hptype, path) {
    bashFile <- file.path(path, paste0("run_remap", hptype, ".sh"))
    rFile    <- file.path(libpath, paste0("remap", hptype, ".R"))
    script   <- paste0(
'#!/usr/bin/env Rscript
library(DR2S)
## get the scripts dir for change wd
args <- commandArgs(trailingOnly = FALSE)
fileArgName <- "--file="
scriptName <- sub(fileArgName, "", args[grep(fileArgName, args)])
scriptBaseName <- dirname(scriptName)
setwd(scriptBaseName)

x <- read_dr2s("..")
tryCatch({refineAlignment(x, "', hptype, '")
         },
         error = function(e) {
           system(paste(
             shQuote("notify-send"),
             shQuote("-u"), shQuote("critical"),
             shQuote(e),
             shQuote("Run in terminal to see whats wrong")))
         })
')
    write(script, file.path(path, rFile))
    write(paste0("#!/usr/bin/env bash\nRscript ", rFile), bashFile)
    Sys.chmod(bashFile, mode = "775")
  }
  invisible(lapply(haptypes, function(hp) writeScript(hp, path)))
}
