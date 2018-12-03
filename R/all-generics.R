

# constructor -------------------------------------------------------------


#' Create a \code{\link[=DR2S_]{DR2S}} configuration object.
#'
#' @section locus:
#' The locus can be specified generically, e.g. "A", "B", or "2DL1". In this
#' case a generic, locus-specific reference will initially be used.
#' Alternatively a locus can be specified more specifically by including
#' allele code information (if available), e.g. "DPB1*04:01:01:01". In this
#' case the specified allele will be used as initial reference.
#'
#' @section datadir:
#' A \code{datadir} must contain subdirectories (default: "pacbio" and
#' "illumina") containig long read FASTQs/FASTAs in the format
#' \code{SAMPLE_LOCUS_*.fast(q|a)(.gz)} and short read FASTQs/FASTAs in the
#' format \code{SAMPLE_LOCUS_*R[12]*.fast(q|a)(.gz)}, respectively.
#'
#' @section Outdir:
#' All output will be placed in directory hierarchy \code{OUTDIR/SAMPLE/LOCUS}
#'
#' @section Reference:
#' References can be specified as a path to a fasta file containing the
#' reference sequence.
#'
#' @param sample A unique sample identifier used to locate the long and short
#' read FASTQ/FASTA files.
#' @param locus The HLA or KIR locus (e.g., "A", "DPB1", or "2DL1"; see Note).
#' @param longreads Location, type, and mapper for longreads as a named list
#' with the fields \code{dir}, \code{type} ("pacbio" or "nanopore") and
#' \code{mapper} ("bwamem" or "minimap").
#' @param shortreads (optional) Location, type, and mapper for short reads
#' as a named list with the fields \code{dir}, \code{type} ("illumina") and
#' \code{mapper} ("bwamem" or "minimap").
#' @param datadir The data directory (See Note).
#' @param outdir The output directory (See Note).
#' @param reference (optional) Path to reference sequence (See Note).
#' @param details <named list> of sample metadata or \code{NULL}. Will be written
#' into the fasta header of the final sequences and stored in the config json
#' @param opts <named list> of arguments to the DR2S pipeline functions or
#'  \code{NULL}. Will be stored in the config json.
#' @param ... Further arguments.
#'
#' @return A \code{\link[DR2S]{DR2S}} config object
#' @family DR2S mapper functions
#' @export
#' @examples
#' \dontrun{
#' x <- InitDR2S(createDR2SConf(
#'   sample = "ID12300527",
#'   locus = "DPB1*04:01:01:01",
#'   datadir = "/path/to/data",
#'   outdir = "/path/to/output"
#'   )) %>%
#'   mapInit() %>%
#'   partitionLongreads() %>%
#'   mapIter() %>%
#'   partitionShortreads() %>%
#'   mapFinal() %>%
#'   polish() %>%
#'   report(blockWidth = 60)
#' }
createDR2SConf <- function(
  sample,
  locus,
  longreads  = list(dir = "pacbio", type = "pacbio", mapper = "minimap"),
  shortreads = list(dir = "illumina", type = "illumina", mapper = "bwamem"),
  datadir    = ".",
  outdir     = "./output",
  reference  = NULL,
  details    = NULL,
  opts       = NULL,
  ...) {
  UseMethod("DR2SConf")
}

#' Initialise a \code{\link[=DR2S_]{DR2S}} mapper.
#' @param config A DR2S configuration object. Either created by
#' \code{\link[=DR2S_]{createDR2SConfig}} or loaded by
#' \code{\link[=DR2S_]{readDR2SConf}}.
#' @param createOutdir Create the outdir if it doesn't exists.
#' @return A \code{\link[=DR2S_]{DR2S}} object.
#' @export
InitDR2S <- function(config, createOutdir = TRUE) {
  UseMethod("InitDR2S")
}


# mappers -----------------------------------------------------------------


#' Initial mapping.
#'
#' Short reads (Illumina) and long reads (PacBio or Nanopore) or long reads alone
#' are mapped to a single reference sequence.
#'
#' @param x A \code{\link[=DR2S_]{DR2S}} object.
#' @param opts List with options passed to the mapper.
#' @param ... Additional parameters passed to \code{\link[Rsamtools]{PileupParam}}.
#' @inheritParams Rsamtools::PileupParam
#' @return A \code{\link[=DR2S_]{DR2S}} object.
#' @family DR2S mapper functions
#' @export
#' @examples
#' \dontrun{
#' x <- DR2Smap(
#'   sample = "ID12300527",
#'   locus = "DPB1*04:01:01:01",
#'   datadir = "/path/to/data",
#'   outdir = "/path/to/output"
#'   ) %>%
#'   mapInit() %>%
#'   partitionLongreads() %>%
#'   mapIter() %>%
#'   partitionShortreads() %>%
#'   mapFinal() %>%
#'   polish() %>%
#'   report(blockWidth = 60)
#' }
mapInit <- function(x, opts = list(), ...) {
  UseMethod("mapInit")
}

#' Iterative mapping.
#'
#' Partitioned long reads are mapped against the MSA consensus sequence of
#' each haplotype. Indels are \strong{excluded} from pileup.
#' New reference consensus sequences for both alleles are inferred.
#'
#' @param x A \code{\link[=DR2S_]{DR2S}} object.
#' @param opts Named list with options passed to the mapper.
#' @param ... Additional parameters passed to \code{\link[Rsamtools]{PileupParam}}.
#' @inheritParams Rsamtools::PileupParam
#' @return A \code{\link[=DR2S_]{DR2S}} object.
#' @family DR2S mapper functions
#' @export
#' @examples
#' \dontrun{
#' x <- DR2Smap(
#'   sample = "ID12300527",
#'   locus = "DPB1*04:01:01:01",
#'   datadir = "/path/to/data",
#'   outdir = "/path/to/output"
#'   ) %>%
#'   mapInit() %>%
#'   partitionLongreads() %>%
#'   mapIter() %>%
#'   partitionShortreads() %>%
#'   mapFinal() %>%
#'   polish() %>%
#'   report(blockwidth = 60)
#' }
mapIter <- function(x, opts = list(), ...) {
  UseMethod("mapIter")
}

#' Final mapping (long reads and short reads).
#'
#' Partitioned long reads and partitioned or unpartitioned short reads are
#' mapped against the consensus sequences inferred at \code{mapIter}. Indels
#' are \strong{included} in pileup.
#'
#' @param x A \code{\link[=DR2S_]{DR2S}} object.
#' @param opts List with options passed to the mapper.
#' @param ... Additional parameters passed to \code{\link[Rsamtools]{PileupParam}}.
#' @inheritParams Rsamtools::PileupParam
#' @return A \code{\link[=DR2S_]{DR2S}} object.
#' @family DR2S mapper functions
#' @export
#' @examples
#' \dontrun{
#' x <- DR2Smap(
#'   sample = "ID12300527",
#'   locus = "DPB1*04:01:01:01",
#'   datadir = "/path/to/data",
#'   outdir = "/path/to/output"
#'   ) %>%
#'   mapInit() %>%
#'   partitionLongreads() %>%
#'   mapIter() %>%
#'   partitionShortreads() %>%
#'   mapFinal() %>%
#'   polish() %>%
#'   report(blockWidth = 60)
#' }
mapFinal <- function(x, opts = list(), ...) {
  UseMethod("mapFinal")
}


# partitioning ------------------------------------------------------------


#' Partition mapped long reads into haplotypes
#'
#' @param x A \code{\link[=DR2S_]{DR2S}} object.
#' @param plot Generate diagnostic plots.
#' @param ... Further arguments passed to methods.
#' @return A \code{\link[=DR2S_]{DR2S}} object.
#' @family DR2S mapper functions
#' @export
#' @examples
#' \dontrun{
#' x <- DR2Smap(
#'   sample = "ID12300527",
#'   locus = "DPB1*04:01:01:01",
#'   datadir = "/path/to/data",
#'   outdir = "/path/to/output"
#'   ) %>%
#'   mapInit() %>%
#'   partitionLongreads() %>%
#'   mapIter() %>%
#'   partitionShortreads() %>%
#'   mapFinal() %>%
#'   polish() %>%
#'   report(blockWidth = 60)
#' }
partitionLongreads <- function(x, plot = TRUE, ...) {
  UseMethod("partitionLongreads")
}

#' Assign short reads from mapInit to haplotypes
#'
#' @param x A \code{\link[=DR2S_]{DR2S}} object.
#' @param opts list with options passed to the mapper.
#' @param ... Further arguments passed to methods.
#' @return A \code{\link[=DR2S_]{DR2S}} object.
#' @family DR2S partition functions
#' @export
#' @examples
#' \dontrun{
#' x <- DR2Smap(
#'   sample = "ID12300527",
#'   locus = "DPB1*04:01:01:01",
#'   datadir = "/path/to/data",
#'   outdir = "/path/to/output"
#'   ) %>%
#'   clear() %>%
#'   mapInit() %>%
#'   partitionLongreads() %>%
#'   mapIter() %>%
#'   partitionShortreads() %>%
#'   mapFinal() %>%
#'   polish() %>%
#'   report(blockWidth = 60)
#' }
partitionShortreads <- function(x, opts = list(), ...) {
  UseMethod("partitionShortreads")
}


# polishing and reporting ------------------------------------------------


#' Polish final haplotype sequences.
#'
#' @param x A \code{\link[=DR2S_]{DR2S}} object.
#' @param ... Additional arguments passed to methods.
#' @return A \code{\link[=DR2S_]{DR2S}} object with the \code{consensus} field
#' populated.
#' @family DR2S mapper functions
#' @export
#' @examples
#' \dontrun{
#' x <- DR2Smap(
#'   sample = "ID12300527",
#'   locus = "DPB1*04:01:01:01",
#'   datadir = "/path/to/data",
#'   outdir = "/path/to/output"
#'   ) %>%
#'   mapInit() %>%
#'   partitionLongreads() %>%
#'   mapIter() %>%
#'   partitionShortreads() %>%
#'   mapFinal() %>%
#'   polish() %>%
#'   report(blockwidth = 60)
#' }
polish <- function(x, ...) {
  UseMethod("polish")
}


#' Report the final haplotype sequences.
#'
#' @param x A \code{\link[=DR2S_]{DR2S}} object.
#' @param whichMap Which mapping do we want to report. Will choose
#' \code{mapFinal} or \code{mapIter} in this order as available.
#' @param ... Additional arguments passed to methods.
#' @details
#' \code{report} will create a \code{report} directory containing the consensus
#' sequences \code{map\{1,2,3\}.\{A,B\}.\{readtype\}.\{mapper\}.unchecked.fa},
#' an html alignment file
#' \code{map\{1,2,3\}.aln.\{readtype\}.\{mapper\}.unchecked.html},
#' a pairwise alignment file
#' \code{map\{1,2,3\}.\{A,B\}.\{readtype\}.\{mapper\}.unchecked.pair},
#' and a tab-separated file
#' \code{problems.\{readtype\}.\{mapper\}.tsv} reporting
#' positions with potentially problematic results.
#'
#' To finalise the consensus sequences open the pairwise alignment file and save
#' it replacing \strong{unchecked} with \strong{checked} in the file name.
#' In the \strong{checked} file manual edits can be performed. After the
#' pairwise alignment is verified call \code{\link{reportCheckedConsensus}}.
#'
#' @return Side effect
#' @family DR2S mapper functions
#' @export
#' @examples
#' \dontrun{
#' x <- DR2Smap(
#'   sample = "ID12300527",
#'   locus = "DPB1*04:01:01:01",
#'   datadir = "/path/to/data",
#'   outdir = "/path/to/output"
#'   ) %>%
#'   mapInit() %>%
#'   partitionLongreads() %>%
#'   mapIter() %>%
#'   partitionShortreads() %>%
#'   mapFinal() %>%
#'   polish() %>%
#'   report(blockWidth = 60)
#' }
report <- function(x, whichMap, ...) {
  UseMethod("report")
}


# helpers -----------------------------------------------------------------


#' Clear the output directory of a DR2S object.
#'
#' @param x A \code{\link[=DR2S_]{DR2S}} object.
#' @param ... Further arguments passed to methods.
#'
#' @return A \code{\link[=DR2S_]{DR2S}} object.
#' @family DR2S mapper functions
#' @export
#' @examples
#' \dontrun{
#' x <- DR2Smap(
#'   sample = "ID12300527",
#'   locus = "DPB1*04:01:01:01",
#'   datadir = "/path/to/data",
#'   outdir = "/path/to/output"
#'   ) %>%
#'   clear() %>%
#'   mapInit() %>%
#'   partitionLongreads() %>%
#'   mapIter() %>%
#'   partitionShortreads() %>%
#'   mapFinal() %>%
#'   polish() %>%
#'   report(blockWidth = 60)
#' }
clear <- function(x, ...) {
  UseMethod("clear")
}

#' Save a DR2S object.
#'
#' @param x A \code{\link[=DR2S_]{DR2S}} object.
#' @param ... Further arguments passed to methods.
#'
#' @return A \code{\link[=DR2S_]{DR2S}} object.
#' @family DR2S mapper functions
#' @export
#' @examples
#' \dontrun{
#' x <- DR2Smap(
#'   sample = "ID12300527",
#'   locus = "DPB1*04:01:01:01",
#'   datadir = "/path/to/data",
#'   outdir = "/path/to/output"
#'   ) %>%
#'   clear() %>%
#'   mapInit() %>%
#'   partitionLongreads() %>%
#'   mapIter() %>%
#'   partitionShortreads() %>%
#'   mapFinal() %>%
#'   polish() %>%
#'   report(blockWidth = 60)
#' }
cache <- function(x, ...) {
  UseMethod("cache")
}

#' Access the name tag of an object.
#'
#' @param x An object containing name tags.
#' @param ... Additional arguments.
#' @return A character vector.
#' @export
#' @examples
#' ###
tag <- function(x, ...) {
  UseMethod("tag")
}

#' Access the path to a bam alignment file.
#'
#' @param x An object containing paths to bam files.
#' @param ... Additional arguments.
#' @return A character vector.
#' @export
#' @examples
#' ###
bampath <- function(x, ...) {
  UseMethod("bampath")
}

#' Access the path to a reference sequence file.
#'
#' @param x An object containing paths to reference files.
#' @param ... Additional arguments.
#' @return A character vector.
#' @export
#' @examples
#' ###
refpath <- function(x, ...) {
  UseMethod("refpath")
}

#' Access the path to a read sequence file.
#'
#' @param x An object containing paths to read files.
#' @param ... Additional arguments.
#' @return A character vector.
#' @export
#' @examples
#' ###
readpath <- function(x, ...) {
  UseMethod("readpath")
}

#' Access the path to a consensus sequence file.
#'
#' @param x An object containing paths to consensus files.
#' @param ... Additional arguments.
#' @return A character vector.
#' @export
#' @examples
#' ###
conspath <- function(x, ...) {
  UseMethod("conspath")
}

#' Access the reference name.
#'
#' @param x An object containing reference name information.
#' @param ... Additional arguments.
#' @return A character vector.
#' @export
#' @examples
#' ###
refname <- function(x, ...) {
  UseMethod("refname")
}

#' Access the readtype.
#'
#' @param x An object containing readtype information.
#' @param ... Additional arguments.
#' @return A character vector.
#' @export
#' @examples
#' ###
readtype <- function(x, ...) {
  UseMethod("readtype")
}

#' Access the stats of an object.
#'
#' @param x An object containing stats
#' @param name (Optional); slotname.
#' @param ... Additional arguments.
#' @return A named list.
#' @export
#' @examples
#' ###
stats <- function(x, name = NULL, ...) {
  UseMethod("stats")
}

#' Access the metadata of an object.
#'
#' @param x An object containing metadata
#' @param name (Optional); slotname.
#' @param ... Additional arguments.
#' @return A named list.
#' @export
#' @examples
#' ###
meta <- function(x, name = NULL, ...) {
  UseMethod("meta")
}
