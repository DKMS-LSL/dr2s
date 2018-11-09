
# Class: Pileup -----------------------------------------------------------

## Internal constructor class "pileup"
Pileup_ <- function(...) {
  dots <- list(...)
  out  <- list()
  out$refpath  <- dots$refpath # <character>; path to reference fasta.
  out$bampath  <- dots$bampath # <character>; path to alignment file.
  out$param    <- dots$param   # <PileupParam>
  out$pileup   <- dots$pileup  # <tbl_df> with columns: "seqnames", "pos", "nucleotide", "count".
  out$consmat  <- consmat(dots$pileup, freq = FALSE) # <consmat>
  out$stats    <- dots$stats %||% list() # <named list>
  dots[c("refpath", "bampath", "param", "pileup", "stats")] <- NULL
  out$meta <- compact(mergeList(list(
      refname   = dots$refname,  # <character>; name of reference.
      readtype  = dots$readtype, # <character>; "illumina", "pacbio", "nanopore".
      reads     = dots$reads     # <character|NULL>; names of n top-scoring reads.
    ), dots))

  structure(out, class = c("pileup", "list"))
}

# Methods: Pileup ---------------------------------------------------------

#' @export
print.pileup <- function(x, asString = FALSE, ...) {
  params <- methods::slotNames(x$param)
  names(params) <- params
  values <- lapply(params, methods::slot, object = x$param)
  info <- paste(methods::slotNames(x$param), values, sep = ": ", collapse = "; ")
  msg <- if (asString) "" else
    sprintf("An object of class '%s'.\n", class(x)[1])
  msg <- sprintf("%s Bamfile: %s\n Reference: %s\n Readtype: %s\n %s\n",
                 msg, basename(bampath(x)), basename(refpath(x)), readtype(x),
                 paste0(strwrap(info, initial = "Params: ", exdent = 4), collapse = "\n"))
  if (asString)
    return(msg)
  else {
    cat(msg)
    print(consmat(x), n = 4)
  }
}

#' @export
bampath.pileup <- function(x, ...) {
  x$bampath
}

#' @export
refpath.pileup <- function(x, ...) {
  x$refpath
}

#' @export
stats.pileup <- function(x, name = NULL, ...) {
  if (!is.null(name)) {
    x$stats[[name]]
  } else {
    x$stats
  }
}

#' @export
meta.pileup <- function(x, name = NULL, ...) {
  if (!is.null(name)) {
    x$meta[[name]]
  } else {
    x$meta
  }
}

#' @export
consmat.pileup <- function(x, freq = FALSE, ...) {
  cm <- if (is(x$consmat, "consmat")) x$consmat else consmat(x$pileup)
  consmat(cm, freq = freq)
}

#' @export
`consmat<-.pileup` <- function(x, value) {
  assert_that(is(value, "consmat"))
  x$consmat <- value
  x
}

#' @export
refname.pileup <- function(x, ...) {
  meta(x, "readname")
}

#' @export
readtype.pileup <- function(x, ...) {
  meta(x, "readtype")
}

#' @rdname pileup
#' @export
reads <- function(x, ...) UseMethod("reads")
#' @export
reads.pileup <- function(x, ..) {
  meta(x, "reads")
}

#' @rdname pileup
#' @export
`reads<-` <- function(x, value) UseMethod("reads<-")
#' @export
`reads<-.pileup` <- function(x, value) {
  assert_that(is.character(value))
  x$meta$reads <- value
  x
}

# Class: MapList ----------------------------------------------------------

## Internal onstructor class "MapList"
MapList_ <- function(...) {
  dots <- list(...)
  out  <- list()
  out$readpath  <- dots$readpath # <character>; relative path(s) to reads.
  out$refpath   <- dots$refpath  # <character>; relative path to reference.
  out$bampath   <- dots$bampath  # <character>; relative path to alignment file.
  out$conspath  <- dots$conspath %||% "" # <character>; relative path to consensus sequence file.
  out$pileup    <- dots$pileup   # <pileup>; "pileup" object.
  out$stats     <- dots$stats %||% list() # <named list>
  dots[c("readpath", "refpath", "bampath", "conspath", "pileup", "stats")] <- NULL
  out$meta <- compact(mergeList( # <named list> that may hold additional fields
    list(
      maplabel = dots$maplabel,  # <character>; the label that indentifies the mapping.
      refname  = dots$refname,   # <character>; the label that indentifies the reference.
      mapper   = dots$mapper,    # <character>; the mapper used "minimap", "bwamem".
      opts     = dots$opts %||% list() # <named list>; non-default mapping options
    ), dots))

  structure(out, class = c("MapList", "list"))
}

# Methods: MapList --------------------------------------------------------

#' @export
print.MapList <- function(x, ...) {
  msg <- sprintf("An object of class '%s'\n", class(x)[1])
  msg <- sprintf(
    "%s [Tag]       %s\n [Reference] %s\n [Reads]     %s\n [Alignment] %s\n [Consensus] %s\n [Coverage]  %s\n",
    msg, tag(x),
    basename(refpath(x)),
    paste(basename(readpath(x)), collapse = "\n             "), # can have more than one readfile.
    basename(bampath(x)),
    basename(conspath(x) %||% ""),
    stats(x, "coverage")[["50%"]] %||% ""
  )
  cat(msg)
}

#' @export
tag.MapList <- function(x, ...) {
  # Descriptive tagline of the mapping
  # maplabel <refname> <readtype> <mapper> <opts>
  paste(meta(x, "maplabel"),
        paste0(litArrows(c(refname(x),
                           readtype(x),
                           meta(x, "mapper"),
                           optstring(meta(x, "opts")))),
               collapse = " "))
}

#' @export
readpath.MapList <- function(x, ...) {
  x$readpath
}

#' @export
bampath.MapList <- function(x, ...) {
  x$bampath
}

#' @export
refpath.MapList <- function(x, ...) {
  x$refpath
}

#' @export
conspath.MapList <- function(x, ...) {
  x$conspath
}

#' Access the consensus name.
#'
#' @param x An object containing consensus information.
#' @param ... Additional arguments.
#' @return A character vector.
#' @export
#' @examples
#' ###
consname <- function(x, ...) UseMethod("consname")
#' @export
consname.MapList <- function(x, ...) {
  if (is.null(x$conspath)) {
    return(NULL)
  }
  consname <- sub("\\.fa(s|sta)?$", "", basename(x$conspath))
  if (startsWith(consname, meta(x, "maplabel"))) {
    sub(paste0("^", meta(x, "maplabel"), "\\."), "", consname)
  } else consname
}

#' @export
refname.MapList <- function(x, ...) {
  meta(x, "refname")
}

#' @export
readtype.MapList <- function(x, ...) {
  readtype(x$pileup)
}

#' @export
stats.MapList <- function(x, name = NULL, ...) {
  if (!is.null(name)) {
    x$stats[[name]]
  } else {
    x$stats
  }
}

#' @export
meta.MapList <- function(x, name = NULL, ...) {
  if (!is.null(name)) {
    x$meta[[name]]
  } else {
    x$meta
  }
}

#' @export
consmat.MapList <- function(x, freq = FALSE, ...) {
  cm <- if (is(x$pileup$consmat, "consmat")) x$pileup$consmat else consmat(x$pileup$pileup)
  consmat(cm, freq = freq)
}


# Class: PartList ---------------------------------------------------------




