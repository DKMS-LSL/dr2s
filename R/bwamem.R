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

make_opts <- function(opts) {
  opts[vapply(opts, isTRUE, FALSE)] <- ""
  gsub("\\s+", " ", paste0(sprintf("-%s %s", names(opts), opts), collapse = " "))
}

generate_mapping_commands <- function(mapper,
                                      readtype,
                                      reffile,
                                      readfile,
                                      allele,
                                      opts = list(),
                                      refname = "",
                                      optsname = "",
                                      outdir = "./output") {
  mapper   <- match.arg(mapper, c("graphmap", "bwamem"))
  readtype <- match.arg(readtype, c("pacbio", "nanopore", "illumina"))
  mapfun <- match.fun(paste0(mapper, "_cmd"))
  outdir <- dir_create_if_not_exists(outdir)

  reffile <- normalizePath(reffile, mustWork = TRUE)
  ref <- file.path(outdir, basename(reffile))
  if (reffile != ref) {
    stopifnot(file.copy(reffile, ref, overwrite = TRUE))
  }

  reads <- paste0(wrap(normalizePath(readfile, mustWork = TRUE), "'"), collapse = " ")

  if (mapper == "bwamem") {
    opts <- compact(merge_list(opts, list(
      x = switch(
        readtype,
        pacbio = "pacbio",
        nanopore = "ont2d",
        illumina = NULL
      ),
      ##  Penalty for introducing Softclipping. Our read quality is usually
      ## quite good and clipping results mostly from bad mapping.
      L = switch(
        readtype,
        pacbio = 5,
        nanopore = 5,
        illumina = 15
      ),
      ## Use only reads above this mapping quality. Maximum is 60
      T = switch(
        readtype,
        pacbio = 50,
        nanopore = 30,
        illumina = 50
      )
    )))
    zip <- ".gz"
  } else {
    zip <- ""
  }

  ## mappers need generally three files: reffile, readfile, outfile
  ## reffile and outfile format
  ## [prefix.]allele.readtype.mapper.[refname.][optsname.][suffix.]ext
  # workaround for these damn windows filename conventions
  allele_nm  <- gsub("[*]", "#", gsub("[:]", "_", paste0(allele, collapse = "~")))
  #allele_nm <- paste0(allele, collapse = "~")
  allelename <- sprintf("%s%s.%s.", allele_nm %+% ".", readtype, mapper)
  refname <- refname %+% "."
  optsname <- optsname %+% "."

  paths <- list(
    reffile  = ref,
    readfile = reads,
    outfile  = normalizePath(file.path(
      outdir, sprintf("%s%s%ssam%s", allelename, refname, optsname, zip)
    ), mustWork = FALSE)
  )

  list(
    cmd      = mapfun(paths, opts),
    reffile  = paths$reffile,
    readfile = paths$readfile,
    outfile  = paths$outfile
  )
}
