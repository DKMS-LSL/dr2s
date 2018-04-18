#!/usr/bin/env Rscript

#' @name dr2s
#'
#' @usage dr2s [-h] <config>
#'
#' @param config yaml run configuration
#'
#' @author Gerhard Schöfl <schoefl@dkms-lab.de>
#' @date 2017-01-05
#' @version 0.1

## Dependencies ####
suppressPackageStartupMessages(stopifnot(
  require("optparse", quietly = TRUE),
  require("yaml", quietly = TRUE),
  require("foreach", quietly = TRUE),
  suppressMessages(require("ipd.Hsapiens.db", quietly = TRUE)),
  suppressMessages(require("DR2S", quietly = TRUE))
))

## Options ####
option_list <- list()
oparser <- OptionParser(usage = "%prog [-h] <config>",
                        option_list,
                        epilogue = '')
arguments   <- parse_args(oparser, positional_arguments = TRUE)
config_file <- arguments$args
configs     <- expand_dr2s_conf(read_dr2s_conf(config_file))

configs[[1]]
# config_file <- "./tests/config1.yaml"
rs <- foreach(conf = configs) %do% {
  mapper <- DR2Smap(conf)
  cat("\nRunning\n", sep = "")
  print(mapper$getConfig())
  cat("\n")
  print(mapper)
  cat("\n")
  Sys.sleep(4)
  mapper$runPipeline()
}

quit(status = 0)
