.onLoad <- function(libname, pkgname) {
  if (.Platform$OS.type != "windows") {
    ## load python funs
    stopifnot(requireNamespace("rPython", quietly = TRUE))
    pyutils <- system.file("python", "pyutils.py", package = "DR2S", mustWork = TRUE)
    rPython::python.exec(readLines(pyutils))
  }
}

#' Initialize DR2S commandline scripts.
#'
#' @param bin path to bin folder
#'
#' @export
setup <- function(bin = "~/bin") {
  if (!file.exists(bin)) {
    dir.create(bin)
  }
  scripts <-
    Sys.glob(paste0(system.file("scripts", package = "DR2S", mustWork = TRUE), "/*"))
  dest <-
    file.path(
      bin,
      vapply(strsplit(basename(scripts), ".", fixed = TRUE), `[`, 1, FUN.VALUE = "")
    )
  foreach(s = scripts, d = dest) %do% {
    if (!file.exists(d)) {
      file.symlink(s, d)
      NULL
    } else if (Sys.readlink(d) == s) {
      NULL
    } else {
      tmp <- c("A utility ", basename(d), " already exists in ", bin,
               ".\nConsider removing/renaming before running setup()")
      message(tmp)
      NULL
    }
  }
  invisible(NULL)
}

