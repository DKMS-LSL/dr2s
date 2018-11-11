

# constructor -------------------------------------------------------------


#' Create a \code{\link[=DR2S_]{DR2S}} configuration object.
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
#' References can be specified as allele codes or a path to a fasta file
#' containing the reference sequence.
#'
#' @param sample A unique sample identifier used to locate the long and short
#' read FASTQ/FASTA files.
#' @param locus The HLA or KIR locus (e.g., "A", "DPB1", or "2DL1")
#' @param longreads Location, type, and mapper for longreads as a named list
#' with the fields \code{dir}, \code{type} ("pacbio" or "nanopore") and
#' \code{mapper} ("bwamem" or "minimap").
#' @param shortreads (optional) Location, type, and mapper for short reads
#' as a named list with the fields \code{dir}, \code{type} ("illumina") and
#' \code{mapper} ("bwamem" or "minimap").
#' @param datadir The data directory (See Note).
#' @param outdir The output directory (See Note).
#' @param reference The reference allele (See Note).
#' @param threshold Threshold frequency for detecting polymorphisms.
#' @param iterations Number of iterations of the mapIter step.
#' @param microsatellite <\code{FALSE}>; Perform a second mapping of shortreads
#' to the inferred reference in mapInit. Set to TRUE if you suspect
#' microsatellites or other repetitive regions in your sequence. Usually extends
#' the reference to a maximum length and enables a better mapping.
#' @param distAlleles Number of different alleles in the sample. Should be 2
#' for heterozygous samples, 1 for homozygous samples and > 2 for some KIR loci.
#' @param filterScores use only reads passing a strict filtering step. TODO
#' @param forceMapping <\code{FALSE}>; set to TRUE if you want to force
#' processing of "bad" shortreads, i.e. when the distribution of coverage is
#' heavily unequal. Aborts the program if maximum coverage > 75 \% quantile * 5.
#' @param details Named list of sample metadata or \code{NULL}. Will be written
#' into the fasta header of the final sequences and stored in the config yaml.
#' @param ... Additional arguments.
#' @return A \code{\link[DR2S]{DR2S}} config object
#' @family DR2S mapper functions
#' @export
#' @importFrom S4Vectors metadata metadata<-
#' @examples
#' \dontrun{
#' x <- InitDR2S(createDR2SConf(
#'   sample = "ID12300527",
#'   locus = "DPB1",
#'   datadir = "/path/to/data",
#'   outdir = "/path/to/output",
#'   reference = "04:01:01:01"
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
  longreads      = list(dir = "pacbio", type = "pacbio", mapper = "minimap"),
  shortreads     = list(dir = "illumina", type = "illumina", mapper = "bwamem"),
  datadir        = ".",
  outdir         = "./output",
  reference      = NULL,
  threshold      = 0.20,
  iterations     = 1,
  microsatellite = FALSE,
  distAlleles    = 2,
  filterScores   = TRUE,
  forceMapping   = FALSE,
  details        = NULL,
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
#' @param opts Mapper options.
#' @param threshold Threshold to call a variant.
#' @param includeDeletions If \code{TRUE}, include deletions in pileup.
#' @param includeInsertions If \code{TRUE}, include insertions in pileup.
#' @param microsatellite If \code{TRUE} remap the shortreads again to expand
#' the reference to a maximum length. Works better for larger insertions and
#' repeats like microsatellites.
#' @param filterScores Apply a harsh filtering on the reads. Filters for mapping
#' quality and removes all softclipping reads. Usually not necessary.
#' @param forceMapping if \code{FALSE} the program throws an error in case the
#' coverage of shortreads is too different at different parts of the sequence.
#' @param topx Select the x best-scoring reads.
#' @param createIgv Subsample for looking with IgvJs in the shiny app.
#' @param force If \code{TRUE}, overwrite existing bam files.
#' @param plot Produce diagnostic plots.
#' @param ... Additional parameters passed to \code{\link[Rsamtools]{PileupParam}}.
#' @inheritParams Rsamtools::PileupParam
#' @return A \code{\link[=DR2S_]{DR2S}} object.
#' @family DR2S mapper functions
#' @export
#' @examples
#' \dontrun{
#' x <- DR2Smap(
#'   sample = "ID12300527",
#'   locus = "DPB1",
#'   datadir = "/path/to/data",
#'   outdir = "/path/to/output",
#'   reference = "04:01:01:01"
#'   ) %>%
#'   mapInit() %>%
#'   partitionLongreads() %>%
#'   mapIter() %>%
#'   partitionShortreads() %>%
#'   mapFinal() %>%
#'   polish() %>%
#'   report(blockWidth = 60)
#' }
mapInit <- function(x,
                    opts = list(),
                    threshold = NULL,
                    includeDeletions = TRUE,
                    includeInsertions = TRUE,
                    microsatellite = FALSE,
                    filterScores = TRUE,
                    forceMapping = FALSE,
                    topx = 0,
                    createIgv = TRUE,
                    force = FALSE,
                    plot = TRUE,
                    ...) {
  UseMethod("mapInit")
}

#' Iterative mapping.
#'
#' Partitioned long reads are mapped against the MSA consensus sequence of
#' each haplotype. Indels are \strong{excluded} from pileup.
#' New reference consensus sequences for both alleles are inferred.
#'
#' @param x A \code{\link[=DR2S_]{DR2S}} object.
#' @param opts Mapper options.
#' @param iterations Number of \code{mapIter} iterations. How often are the
#' clustered reads remapped to updated reference sequences.
#' @param columnOccupancy Minimum occupancy (1 - fraction of gap) below which
#' bases at insertion position are excluded from from consensus calling.
#' @param force If \code{TRUE}, overwrite existing bam file.
#' @param plot Plot diagnostics.
#' @param ... Additional parameters passed to \code{\link[Rsamtools]{PileupParam}}.
#' @inheritParams Rsamtools::PileupParam
#' @return A \code{\link[=DR2S_]{DR2S}} object.
#' @family DR2S mapper functions
#' @export
#' @examples
#' \dontrun{
#' x <- DR2Smap(
#'   sample = "ID12300527",
#'   locus = "DPB1",
#'   datadir = "/path/to/data",
#'   outdir = "/path/to/output",
#'   reference = "04:01:01:01"
#'   ) %>%
#'   mapInit() %>%
#'   partitionLongreads() %>%
#'   mapIter() %>%
#'   partitionShortreads() %>%
#'   mapFinal() %>%
#'   polish() %>%
#'   report(blockwidth = 60)
#' }
mapIter <- function(x,
                    opts = list(),
                    iterations = 1,
                    columnOccupancy = 0.4,
                    force = FALSE,
                    plot = TRUE,
                    ...) {
  UseMethod("mapIter")
}

#' Final mapping (long reads and short reads).
#'
#' Partitioned long reads and partitioned or unpartitioned short reads are
#' mapped against the consensus sequences inferred at \code{mapIter}. Indels
#' are \strong{included} in pileup.
#'
#' @param x A \code{\link[=DR2S_]{DR2S}} object.
#' @param opts Mapper options.
#' @param includeDeletions If \code{TRUE}, include deletions in pileup.
#' @param includeInsertions If \code{TRUE}, include insertions in pileup.
#' from bam files for
#'  variants ar polymorphic positions.
#' @param clip Clip homopolymeric read ends from shortreads.
#' @param force If \code{TRUE}, overwrite existing bam file.
#' @param createIgv Subsample for looking with IgvJs in the shiny app.
#' @param plot Plot diagnostics.
#' @param ... Additional parameters passed to \code{\link[Rsamtools]{PileupParam}}.
#' @inheritParams Rsamtools::PileupParam
#' @return A \code{\link[=DR2S_]{DR2S}} object.
#' @family DR2S mapper functions
#' @export
#' @examples
#' \dontrun{
#' x <- DR2Smap(
#'   sample = "ID12300527",
#'   locus = "DPB1",
#'   datadir = "/path/to/data",
#'   outdir = "/path/to/output",
#'   reference = "04:01:01:01"
#'   ) %>%
#'   mapInit() %>%
#'   partitionLongreads() %>%
#'   mapIter() %>%
#'   partitionShortreads() %>%
#'   mapFinal() %>%
#'   polish() %>%
#'   report(blockWidth = 60)
#' }
mapFinal <- function(x,
                     opts = list(),
                     includeDeletions = TRUE,
                     includeInsertions = TRUE,
                     clip = TRUE,
                     force = FALSE,
                     createIgv = TRUE,
                     plot = TRUE,
                     ...) {
  UseMethod("mapFinal")
}


# partitioning ------------------------------------------------------------


#' Partition mapped long reads into haplotypes
#',

#' @param x A \code{\link[=DR2S_]{DR2S}} object.
#' @param threshold The threshold when a SNP is a SNP.
#' @param skipGapFreq The gap frequenzy needed to call a gap position.
#' @param distAlleles The number of distinct alleles in the sample.
#' @param noGapPartitioning Don't partition based on gaps. Useful for samples
#' with only few SNPs but with homopolymers. The falsely called gaps could
#' mask the real variation.
#' @param selectAllelesBy If more than \code{distAlleles} clusters are found
#' select clusters based on: (1) "distance": The hamming distance of the
#' resulting variant consensus sequences or (2) "count": Take the clusters
#' with the most reads as the true alleles.
#' @param plot Plot
#' @param ... Further arguments passed to methods.
#' @return A \code{\link[=DR2S_]{DR2S}} object.
#' @family DR2S mapper functions
#' @export
#' @examples
#' \dontrun{
#' x <- DR2Smap(
#'   sample = "ID12300527",
#'   locus = "DPB1",
#'   datadir = "/path/to/data",
#'   outdir = "/path/to/output",
#'   reference = "04:01:01:01"
#'   ) %>%
#'   mapInit() %>%
#'   partitionLongreads() %>%
#'   mapIter() %>%
#'   partitionShortreads() %>%
#'   mapFinal() %>%
#'   polish() %>%
#'   report(blockWidth = 60)
#' }
partitionLongreads <- function(x,
                               threshold         = NULL,
                               skipGapFreq       = 2/3,
                               distAlleles       = NULL,
                               noGapPartitioning = FALSE,
                               selectAllelesBy   = "count",
                               plot              = TRUE,
                               ...) {
  UseMethod("partitionLongreads")
}

#' Assign short reads from mapInit to haplotypes
#'
#' @param x A \code{\link[=DR2S_]{DR2S}} object.
#' @param opts list with options passed to the mapper.
#' @param force force the creation of new fastq files if they already exist.
#' @param ... Further arguments passed to methods.
#' @return A \code{\link[=DR2S_]{DR2S}} object.
#' @family DR2S partition functions
#' @export
#' @examples
#' \dontrun{
#' x <- DR2Smap(
#'   sample = "ID12300527",
#'   locus = "DPB1",
#'   datadir = "/path/to/data",
#'   outdir = "/path/to/output",
#'   reference = "04:01:01:01"
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
partitionShortreads <- function(x,
                                opts = list(),
                                force = TRUE,
                                ...) {
  UseMethod("partitionShortreads")
}


# polishing and reporting ------------------------------------------------


#' Polish final haplotype sequences.
#'
#' @param x A \code{\link[=DR2S_]{DR2S}} object.
#' @param threshold When do we call a variant a variant.
#' @param checkHpCount Check the number of homopolymer counts in shortreads.
#' Compare the resulting sequence with the mode value and report differences.
#' @param cache Cache the updated \code{DR2S} object after assembling the
#' haplotypes.
#' @param ... Additional arguments passed to methods.
#' @return A \code{\link[=DR2S_]{DR2S}} object with the \code{consensus} field
#' populated.
#' @family DR2S mapper functions
#' @export
#' @examples
#' \dontrun{
#' x <- DR2Smap(
#'   sample = "ID12300527",
#'   locus = "DPB1",
#'   datadir = "/path/to/data",
#'   outdir = "/path/to/output",
#'   reference = "04:01:01:01"
#'   ) %>%
#'   mapInit() %>%
#'   partitionLongreads() %>%
#'   mapIter() %>%
#'   partitionShortreads() %>%
#'   mapFinal() %>%
#'   polish() %>%
#'   report(blockwidth = 60)
#' }
polish <- function(x,
                   threshold = x$getThreshold(),
                   checkHpCount = TRUE,
                   cache = TRUE,
                   ...) {
  UseMethod("polish")
}


#' Report the final haplotype sequences.
#'
#' @param x A \code{\link[=DR2S_]{DR2S}} object.
#' @param which Which mapping do we want to report. Will choose \code{mapFinal},
#' \code{mapIter} in this order if available.
#' @param blockWidth Maximum number of sequence letters per line in pairwise
#' alignment.
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
#'   locus = "DPB1",
#'   datadir = "/path/to/data",
#'   outdir = "/path/to/output",
#'   reference = "04:01:01:01"
#'   ) %>%
#'   mapInit() %>%
#'   partitionLongreads() %>%
#'   mapIter() %>%
#'   partitionShortreads() %>%
#'   mapFinal() %>%
#'   polish() %>%
#'   report(blockWidth = 60)
#' }
report <- function(x, which, blockWidth = 80, ...) {
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
#'   locus = "DPB1",
#'   datadir = "/path/to/data",
#'   outdir = "/path/to/output",
#'   reference = "04:01:01:01"
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
#'   locus = "DPB1",
#'   datadir = "/path/to/data",
#'   outdir = "/path/to/output",
#'   reference = "04:01:01:01"
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
