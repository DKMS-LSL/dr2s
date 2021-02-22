## Commands for minimap2
runminimap <- function(reffile,       ## <character>; file path to reference sequence
                       refname,       ## <character>; allele/reference identifier
                       readfile,      ## <character>; file path to read sequences
                       readtype,      ## <character>; "illumina", "pacbio", "nanopore"
                       outdir,        ## <character>; dir path to output directory
                       label = "",    ## <character>; prefix to output file
                       opts = list(), ## <named list>;
                       cores = "auto",
                       ...) {
  if (missing(refname)) {
    refname <- "ref"
  }
  if (missing(readtype)) {
    readtype <- "illumina"
  }
  if (cores == "auto") {
    cores <- .getIdleCores()
  }
  assert_that(is.numeric(cores))
  opts <- mergeList(opts, list(t = cores))

  cmd  <- Sys.which("minimap2")
  if (!nzchar(cmd))
    cmd <- system.file("extSrc/minimap2", package = "DR2S")
  cmd <- normalizePath(cmd, mustWork = TRUE)
  args <- .generateMappingCommands("minimap", reffile, refname,
                                   readfile, readtype, outdir,
                                   label, opts)
  ## Don't execute if file exists
  if (!file.exists(args$outfile)) {
    system2(cmd, args$args)
  } else {
    warning(sprintf("file <$s> already exists", args$outfile))
  }
  args$outfile
}

minimapCmd <- function(paths, opts) {
  opts
  sprintf(
    "--end-bonus 1 -a %s '%s' %s | gzip -3 > '%s'",
    .makeOpts(opts),
    paths$reffile,
    paths$readfile,
    paths$outfile
  )
}

## Commands for bwamem
runbwamem <- function(reffile,       ## <character>; file path to reference sequence
                      refname,       ## <character>; allele/reference identifier
                      readfile,      ## <character>; file path to read sequences
                      readtype,      ## <character>; "illumina", "pacbio", "nanopore"
                      outdir,        ## <character>; dir path to output directory
                      label = "",    ## <character>; prefix to output file
                      opts = list(), ## <named list>;
                      cores = "auto",
                      ...) {
  assert_that(.hasCommand("bwa"))
  if (missing(refname)) {
    refname <- "ref"
  }
  if (missing(readtype)) {
    readtype <- "illumina"
  }
  if (cores == "auto") {
    cores <- .getIdleCores()
  }
  assert_that(is.numeric(cores))
  opts <- mergeList(opts, list(t = cores))

  # debug
  # reffile <- self$getRefPath()
  # readfile <- self$getShortreads()
  # refname <- "n"
  # outdir <- self$getOutdir()

  cmd  <- normalizePath(Sys.which("bwa"), mustWork = TRUE)
  args <- .generateMappingCommands("bwamem", reffile, refname,
                                   readfile, readtype, outdir,
                                   label, opts)
  ## Don't execute if file exists
  if (!file.exists(args$outfile)) {
    system2(cmd, args$args$idx)
    system2(cmd, args$args$map)
  } else {
    warning(sprintf("file <$s> already exists", args$outfile))
  }
  ## cleanup
  .fileDeleteIfExists(args$reffile %<<% ".amb")
  .fileDeleteIfExists(args$reffile %<<% ".ann")
  .fileDeleteIfExists(args$reffile %<<% ".bwt")
  .fileDeleteIfExists(args$reffile %<<% ".fai")
  .fileDeleteIfExists(args$reffile %<<% ".pac")
  .fileDeleteIfExists(args$reffile %<<% ".sa")
  ##
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


#' Run rsubread
#' 
runrsubread <- function(reffile,       ## <character>; file path to reference sequence
                        refname,       ## <character>; allele/reference identifier
                        readfile,      ## <character>; file path to read sequences
                        readtype,      ## <character>; "illumina", "pacbio", "nanopore"
                        outdir,        ## <character>; dir path to output directory
                        label = "",    ## <character>; prefix to output file
                        opts = list(), ## <named list>;
                        ...) {
  require(Rsubread)
  
  if (missing(refname)) {
    refname <- "ref"
  }
  if (missing(readtype)) {
    readtype <- "illumina"
  }
  assert_that(readtype == "illumina")
  mapper <- "rsubread"
  optsname <- optstring(opts)
  
  filename <- sprintf("%s.%s.%s.", strip(refname, "_"), readtype, mapper)
  label    <- if (nzchar(label)) label %<<% "." else ""
  optsname <- if (nzchar(optsname)) optsname %<<% "." else ""
  outfile  <- normalizePath(file.path(
      outdir, sprintf("%s%s%ssam", label, filename, optsname)
    ), mustWork = FALSE)

  ## Build the index
  invisible(Rsubread::buildindex(basename = refname, reference = reffile))

  ## map reads
  if (length(readfile) == 2) {
    invisible(Rsubread::align(
      index = refname, 
      readfile1 = readfile[1], 
      readfile2 = readfile[2], 
      type = "dna", 
      output_file = outfile,
      output_format = "BAM", 
      unique = TRUE, 
      sortReadsByCoordinates = TRUE))
  } else {
    invisible(Rsubread::align(
      index = refname, 
      readfile1 = readfile, 
      type = "dna", 
      output_file = outfile,
      output_format = "BAM", 
      unique = TRUE, 
      sortReadsByCoordinates = TRUE))
  }
    
  ## cleanup
  .fileDeleteIfExists(reffile %<<% ".fai")
  .fileDeleteIfExists(outfile %<<% ".indel.vcf")
  .fileDeleteIfExists(outfile %<<% ".summary")
  .fileDeleteIfExists(outfile %<<% ".bai")
  ##
  outfile
}


.makeOpts <- function(opts) {
  opts[vapply(opts, isTRUE, FALSE)] <- ""
  gsub("\\s+", " ", paste0(sprintf("-%s %s", names(opts), opts), collapse = " "))
}

.generateMappingCommands <- function(mapper,
                                     reffile,
                                     refname,
                                     readfile,
                                     readtype,
                                     outdir = "./output",
                                     label,
                                     opts = list()) {
  mapper   <- match.arg(mapper, c("minimap", "bwamem"))
  readtype <- match.arg(readtype, c("pacbio", "nanopore", "illumina"))
  mapfun   <- match.fun(mapper %<<% "Cmd")
  outdir   <- .dirCreateIfNotExists(outdir)
  reffile  <- normalizePath(reffile, mustWork = TRUE)
  ref      <- file.path(outdir, basename(reffile))
  if (reffile != ref) {
    stopifnot(file.copy(reffile, ref, overwrite = TRUE))
  }
  optsname <- optstring(opts)
  reads <- paste0(wrap(normalizePath(readfile, mustWork = TRUE), "'"), collapse = " ")

  if (mapper == "bwamem") {
    opts <- compact(mergeList(opts, list(
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
    opts <- compact(mergeList(opts, list(
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
  ## reffile and outfile format:
  ##   [label.]refname.readtype.mapper.[optsname.][suffix.]ext
  # if refname is a HLA allele, e.g., HLA-A*01:01:01:01, make it compliant with
  # windos filesystem conventions
  filename <- sprintf("%s.%s.%s.", strip(refname, "_"), readtype, mapper)
  label    <- if (nzchar(label)) label %<<% "." else ""
  optsname <- if (nzchar(optsname)) optsname %<<% "." else ""
  paths <- list(
    reffile  = ref,
    readfile = reads,
    outfile  = normalizePath(file.path(
      outdir, sprintf("%s%s%ssam%s", label, filename, optsname, zip)
    ), mustWork = FALSE)
  )

  list(
    args     = mapfun(paths, opts),
    reffile  = paths$reffile,
    readfile = paths$readfile,
    outfile  = paths$outfile
  )
}
