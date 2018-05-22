# config for logging DR2S
.confLog <- function(outdir, logName = "info"){
  stopifnot( requireNamespace("futile.logger", quietly = TRUE) )
  futile.logger::flog.appender(futile.logger::appender.tee(
    file.path(outdir, paste0("DR2S_run.", logName))),
    logName)
}
