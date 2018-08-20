#!/usr/bin/env Rscript
library(DR2S)
## get the scripts dir for change wd
args <- commandArgs(trailingOnly = FALSE)
fileArgName <- "--file="
scriptName <- sub(fileArgName, "", args[grep(fileArgName, args)])
scriptBaseName <- dirname(scriptName)
setwd(scriptBaseName)

x <- readDR2S("..")
tryCatch({refineAlignment(x, "A")
         },
         error = function(e) {
           system(paste(
             shQuote("notify-send"),
             shQuote("-u"), shQuote("critical"),
             shQuote(e),
             shQuote("Run in terminal to see whats wrong")))
         })

