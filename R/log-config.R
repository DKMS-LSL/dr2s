# config for logging DR2S

conf_log <- function(outdir, logName = "info"){
  stopifnot( requireNamespace("futile.logger", quietly = TRUE) )
  flog.appender(appender.tee(file.path(outdir, paste0("DR2S_run.", logName), logName)))
}
