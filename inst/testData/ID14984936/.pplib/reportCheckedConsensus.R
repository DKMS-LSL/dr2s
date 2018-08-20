library(DR2S)
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

