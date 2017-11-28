run_bwamem <- function(reffile,
                       readfile,
                       allele,
                       readtype,
                       opts = list(),
                       refname = "",
                       optsname = "",
                       force = FALSE,
                       outdir,
                       ...) {
  if (missing(allele)) {
    allele <- ""
  }
  if (missing(readtype)) {
    readtype <- "illumina"
  }
  if (!is.null(optsname)) {
    optsname <- gsub("[[:punct:][:space:]]", "", optstring(opts, optsname))
  }
  opts <- merge_list(opts, list(t = parallel::detectCores()/2))

  # debug
  # reffile <- self$getRefPath()
  # readfile <- self$getShortreads()
  # refname <- "n"
  # outdir <- self$getOutdir()

  cmd <- generate_mapping_commands("bwamem", readtype, reffile, readfile,
                                   allele, opts, refname, optsname,
                                   outdir = outdir)
  ## Don't execute if file exists and force is false
  if (force || !file.exists(cmd$outfile)) {
    system(cmd$cmd)
  }
  ## cleanup
  file_delete_if_exists(paste0(cmd$reffile, ".amb"))
  file_delete_if_exists(paste0(cmd$reffile, ".ann"))
  file_delete_if_exists(paste0(cmd$reffile, ".bwt"))
  file_delete_if_exists(paste0(cmd$reffile, ".fai"))
  file_delete_if_exists(paste0(cmd$reffile, ".pac"))
  file_delete_if_exists(paste0(cmd$reffile, ".sa"))

  # debug
  # samfile <- cmd$outfile

  cmd$outfile
}

# Options: -x STR read type (pacbio, ont2d, intractg)
#          -t INT number of threads
bwamem_cmd <- function(paths, opts) {
  exepath <- normalizePath(Sys.which("bwa"), mustWork = TRUE)
  .idx <- sprintf("%s index -a is '%s'", exepath, paths$reffile)
  sprintf(
    "%s && %s mem %s '%s' %s | gzip -3 > '%s'",
    .idx,
    exepath,
    make_opts(opts),
    paths$reffile,
    paths$readfile,
    paths$outfile
  )
}
