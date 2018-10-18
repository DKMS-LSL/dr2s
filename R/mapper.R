## Commands for minimap2
runminimap <- function(reffile,
                       readfile,
                       allele,
                       readtype,
                       opts = list(),
                       refname = "",
                       optsname = "",
                       force = FALSE,
                       outdir,
                       cores = "auto",
                       ...) {
  assert_that(.hasCommand("minimap2"))
  if (missing(allele)) {
    allele <- ""
  }
  if (missing(readtype)) {
    readtype <- "illumina"
  }
  if (!is.null(optsname)) {
    optsname <- gsub("[[:punct:][:space:]]", "", optstring(opts, optsname))
  }
  if (cores == "auto") {
    cores <- .getIdleCores()
  }
  assert_that(is.numeric(cores))

  opts <- .mergeList(opts, list(t = cores))

  cmd  <- Sys.which("minimap2")
  args <- .generateMappingCommands(mapper = "minimap", readtype,
                                   reffile, readfile, allele, opts,
                                   refname, optsname, outdir = outdir)
  ## Don't execute if file exists and force is false
  if (force || !file.exists(args$outfile)) {
    system2(cmd, args$args)
  }
  args$outfile
}

minimapCmd <- function(paths, opts) {
  sprintf(
    "-a %s '%s' %s | gzip -3 > '%s'",
    .makeOpts(opts),
    paths$reffile,
    paths$readfile,
    paths$outfile
  )
}

## Commands for bwamem
runbwamem <- function(reffile,
                      readfile,
                      allele,
                      readtype,
                      opts = list(),
                      refname = "",
                      optsname = "",
                      force = FALSE,
                      outdir,
                      cores = "auto",
                      ...) {
  assert_that(.hasCommand("bwa"))
  if (missing(allele)) {
    allele <- ""
  }
  if (missing(readtype)) {
    readtype <- "illumina"
  }
  if (!is.null(optsname)) {
    optsname <- gsub("[[:punct:][:space:]]", "", optstring(opts, optsname))
  }
  if (cores == "auto") {
    cores <- .getIdleCores()
  }
  assert_that(is.numeric(cores))
  opts <- .mergeList(opts, list(t = cores))

  # debug
  # reffile <- self$getRefPath()
  # readfile <- self$getShortreads()
  # refname <- "n"
  # outdir <- self$getOutdir()

  cmd <- normalizePath(Sys.which("bwa"), mustWork = TRUE)
  args <- .generateMappingCommands("bwamem", readtype, reffile, readfile,
                                   allele, opts, refname, optsname,
                                   outdir = outdir)
  ## Don't execute if file exists and force is false
  if (force || !file.exists(args$outfile)) {
    system2(cmd, args$args$idx)
    system2(cmd, args$args$map)
  }
  ## cleanup
  .fileDeleteIfExists(args$reffile %<<% ".amb")
  .fileDeleteIfExists(args$reffile %<<% ".ann")
  .fileDeleteIfExists(args$reffile %<<% ".bwt")
  .fileDeleteIfExists(args$reffile %<<% ".fai")
  .fileDeleteIfExists(args$reffile %<<% ".pac")
  .fileDeleteIfExists(args$reffile %<<% ".sa")

  args$outfile
}

# Options: -x STR read type (pacbio, ont2d, intractg)
#          -t INT number of threads
bwamemCmd <- function(paths, opts) {
  .idx <- sprintf("index -a is '%s'", paths$reffile)
  map <- sprintf(
    "mem %s '%s' %s | gzip -3 > '%s'",
    .makeOpts(opts),
    paths$reffile,
    paths$readfile,
    paths$outfile
  )
  list(idx = .idx,
       map = map)
}

.makeOpts <- function(opts) {
  opts[vapply(opts, isTRUE, FALSE)] <- ""
  gsub("\\s+", " ", paste0(sprintf("-%s %s", names(opts), opts),
                           collapse = " "))
}

.generateMappingCommands <- function(mapper,
                                     readtype,
                                     reffile,
                                     readfile,
                                     allele,
                                     opts = list(),
                                     refname = "",
                                     optsname = "",
                                     outdir = "./output") {
  mapper   <- match.arg(mapper, c("minimap", "bwamem"))
  readtype <- match.arg(readtype, c("pacbio", "nanopore", "illumina"))
  mapfun   <- match.fun(mapper %<<% "Cmd")
  outdir   <- .dirCreateIfNotExists(outdir)

  reffile <- normalizePath(reffile, mustWork = TRUE)
  ref <- file.path(outdir, basename(reffile))
  if (reffile != ref) {
    stopifnot(file.copy(reffile, ref, overwrite = TRUE))
  }

  reads <- paste0(wrap(normalizePath(readfile, mustWork = TRUE), "'"),
                  collapse = " ")

  if (mapper == "bwamem") {
    opts <- compact(.mergeList(opts, list(
      x = switch(
        readtype,
        pacbio = "pacbio",
        nanopore = "ont2d",
        illumina = NULL
      ),
      # ## Matching score
      c = switch(
        readtype,
        pacbio = NULL,
        nanopore = NULL,
        illumina = 4
      ),
      # # ## Mismatch penalty
      B = switch(
        readtype,
        pacbio = NULL,
        nanopore = NULL,
        illumina = 6
      ),
      r = switch(
        readtype,
        pacbio = NULL,
        nanopore = NULL,
        illumina = 0.1
      ),
      # ##  Penalty for gap opening
      O = switch(
        readtype,
        pacbio = NULL,
        nanopore = "6,6",
        illumina = "18,18"
      ),
      ##  Penalty for introducing Softclipping. Our read quality is usually
      ## quite good and clipping results mostly from bad mapping.
      L = switch(
        readtype,
        pacbio = 5,
        nanopore = 5,
        illumina = 60
      ),
      ## Use only reads above this mapping quality. Maximum is 60
      T = switch(
        readtype,
        pacbio = 50,
        nanopore = 30,
        illumina = 50
      ),
      w = switch(
        readtype,
        pacbio = NULL,
        nanopore = 30,
        illumina = 100
      )
    )))
    zip <- ".gz"
  } else if (mapper == "minimap") {
    opts <- compact(.mergeList(opts, list(
      x = switch(
        readtype,
        pacbio = "map-pb",
        nanopore = "map-ont",
        illumina = "sr"
      ),
      # ## Matching score
      k = switch(
        readtype,
        pacbio = 13,
        nanopore = NULL,
        illumina = 6
      ),
      # # ## Mismatch penalty
      w = switch(
        readtype,
        pacbio = 4,
        nanopore = NULL,
        illumina = 2
      ),
      N = switch(
        readtype,
        pacbio = 50,
        nanopore = NULL,
        illumina = 500
      ),
      A = switch(
        readtype,
        pacbio = NULL,
        nanopore = NULL,
        illumina = 2
      ),
      B = switch(
        readtype,
        pacbio = NULL,
        nanopore = NULL,
        illumina = 8
      ),
      O = switch(
        readtype,
        pacbio = NULL,
        nanopore = NULL,
        illumina = "16,32"
      ),
      E = switch(
        readtype,
        pacbio = NULL,
        nanopore = NULL,
        illumina = "4,2"
      ),
      r = switch(
        readtype,
        pacbio = NULL,
        nanopore = NULL,
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
  alleleNm  <- gsub("[*]", "_",
                    gsub("[:]", "_", paste0(allele, collapse = "~")))
  allelename <- sprintf("%s%s.%s.", alleleNm %<<% ".", readtype, mapper)
  refname <- if (nzchar(refname)) refname %<<% "." else ""
  optsname <- if (nzchar(optsname)) optsname %<<% "." else ""

  paths <- list(
    reffile  = ref,
    readfile = reads,
    outfile  = normalizePath(file.path(
      outdir, sprintf("%s%s%ssam%s", allelename, refname, optsname, zip)
    ), mustWork = FALSE)
  )

  list(
    args      = mapfun(paths, opts),
    reffile  = paths$reffile,
    readfile = paths$readfile,
    outfile  = paths$outfile
  )
}
