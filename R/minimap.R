run_minimap <- function(reffile,
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
  opts <- merge_list(opts, list(t = parallel::detectCores()/4))

  # debug
  # reffile <- self$getRefPath()
  # readfile <- self$getShortreads()
  # refname <- "n"
  # outdir <- self$getOutdir()

  cmd <- generate_mapping_commands("minimap", readtype, reffile, readfile,
                                   allele, opts, refname, optsname,
                                   outdir = outdir)
  ## Don't execute if file exists and force is false
  if (force || !file.exists(cmd$outfile)) {
    system(cmd$cmd)
  }
  cmd$outfile
}

# Options: -x STR read type (map-pb, map-ont, sr)
#          -t INT number of threads
minimap_cmd <- function(paths, opts) {
  exepath <- normalizePath(Sys.which("minimap2"), mustWork = TRUE)
  sprintf(
    "%s -a %s '%s' %s | gzip -3 > '%s'",
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
  mapper   <- match.arg(mapper, c("minimap", "bwamem"))
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
      # ##  Penalty for introducing Softclipping. Our read quality is usually
      # ## quite good and clipping results mostly from bad mapping.
      # O = switch(
      #   readtype,
      #   pacbio = "6,6",
      #   nanopore = "6,6",
      #   illumina = "10,10"
      # ),
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
        illumina = 300
      )
    )))
    zip <- ".gz"
  } else if (mapper == "minimap") {
    opts <- compact(merge_list(opts, list(
      x = switch(
        readtype,
        pacbio = "map-pb",
        nanopore = "map-ont",
        illumina = "sr"
      ),
      # ## Matching score
      k = switch(
        readtype,
        pacbio = NULL,
        nanopore = NULL,
        illumina = 25
      ),
      # # ## Mismatch penalty
      w = switch(
        readtype,
        pacbio = NULL,
        nanopore = NULL,
        illumina = 5
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
        illumina = 16
      ),
      O = switch(
        readtype,
        pacbio = NULL,
        nanopore = NULL,
        illumina = "16,24"
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
        illumina = 20
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

    # opts <- compact(list(
    #   t = 6,
    #   x = switch(
    #     readtype,
    #     pacbio = "pacbio",
    #     nanopore = "ont2d",
    #     illumina = NULL
    #   ),
    #   # ## Matching score
    #   c = switch(
    #     readtype,
    #     pacbio = "1",
    #     nanopore = "1",
    #     illumina = 4
    #   ),
    #   # # ## Mismatch penalty
    #   B = switch(
    #     readtype,
    #     pacbio = "6,6",
    #     nanopore = "6,6",
    #     illumina = 6
    #   ),
    #
    #   # ##  Penalty for gap opening
    #   O = switch(
    #     readtype,
    #     pacbio = "6,6",
    #     nanopore = "6,6",
    #     illumina = "28,28"
    #   ),
    #   # ##  Penalty for introducing Softclipping. Our read quality is usually
    #   # ## quite good and clipping results mostly from bad mapping.
    #   r = switch(
    #     readtype,
    #     pacbio = "6,6",
    #     nanopore = "6,6",
    #     illumina = 0.5
    #   ),
    #   ##  Penalty for introducing Softclipping. Our read quality is usually
    #   ## quite good and clipping results mostly from bad mapping.
    #   L = switch(
    #     readtype,
    #     pacbio = 5,
    #     nanopore = 5,
    #     illumina = 60
    #   ),
    #   ## Use only reads above this mapping quality. Maximum is 60
    #   T = switch(
    #     readtype,
    #     pacbio = 50,
    #     nanopore = 30,
    #     illumina = 50
    #   ),
    #   w = switch(
    #     readtype,
    #     pacbio = 50,
    #     nanopore = 30,
    #     illumina = 300
    #   )
    # ))
