## Make Rscript files for direct postprocessing
##
writeReportCheckedConsensus <- function(path) {
  file <- file.path(path, "reportCheckedConsensus.R")
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
tryCatch(report_checked_consensus(x),
         error = function(e) {
           system(paste(
             shQuote("notify-send"),
             shQuote("-u"), shQuote("critical"),
             shQuote("Report failed! Ambiguous positions found"),
             shQuote("Run in terminal or look at the log to see whats wrong")))
         })
'
  write(script, file)
  Sys.chmod(file, mode = "775")
}

writeCheckConsensus <- function(path) {
  file <- file.path(path, "checkConsensus.R")
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
  write(script, file)
}


writePlotDiagnosticAlignment <- function(path) {
  file <- file.path(path, "plotDiagnosticAlignment.R")
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
tryCatch(plot_diagnostic_alignment(x),
         error = function(e) {
           system(paste(
             shQuote("notify-send"),
             shQuote("-u"), shQuote("critical"),
             shQuote(e),
             shQuote("Run in terminal to see whats wrong")))
         })
'
  write(script, file)
  Sys.chmod(file, mode = "775")
}

writeRefineAlignments <- function(path, haptypes) {
  writeScript <- function(hptype, path) {
    file <- file.path(path, paste0("refineAlignment", hptype, ".R"))
    script <- paste0(
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
         run_igv(x, map = "refine", open_now = FALSE)
         },
         error = function(e) {
           system(paste(
             shQuote("notify-send"),
             shQuote("-u"), shQuote("critical"),
             shQuote(e),
             shQuote("Run in terminal to see whats wrong")))
         })
')
    script
    write(script, file)
    Sys.chmod(file, mode = "775")
  }
  invisible(lapply(haptypes, function(hp) writeScript(hp, path)))
}

