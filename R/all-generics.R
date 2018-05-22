

# constructor -------------------------------------------------------------


#' Constructor for \code{\link[=DR2S_]{DR2S}} mapper objects.
#'
#' Initialise a \code{\link[=DR2S_]{DR2S}} mapper.
#'
#' @section datadir:
#'
#' A \code{datadir} must contain arbitrarily named subdirectories
#' (default: "pacbio" and "illumina") with long read FASTQs in the format
#' \code{SAMPLE_LOCUS_*.fastq(.gz)} and short read FASTQs in the format
#' \code{SAMPLE_LOCUS_*R[12]*.fastq(.gz)}, respectively.
#'
#' @section Outdir:
#'
#' All output will be placed in directory hierarchy
#' \code{OUTDIR/SAMPLE/READTYPE/LOCUS}.
#'
#' @section Reference:
#'
#' References can be specified as allele codes or a path to a fasta file
#' containing the reference sequence. If \code{reference = NULL} a global generic reference
#' for a given locus will be used.
#'
#' @section Consensus:
#' Right now, the consensus method of multialign is removed. A de novo assembly of
#' longreads, maybe supported by shortreads might be implemented.
#' \dQuote{\bold{\code{mapping}}}: The reference is refined from the initial provided reference
#' during the mapInit step using short reads. Individual references for each found haplotype are
#' constructed from this reference using only longreads that are assigned to the haplotype.
#'
#' @param sample A unique sample identifier used to locate the long and short
#' read FASTQ files.
#' @param locus The HLA locus ("A", "B", "C", "DPB1", "DQB1", "DRB1") or KIR locus.
#' @param longreads The type and location of the long reads as a named list
#' with the fields \code{type} ("pacbio" or "nanopore") and \code{dir}.
#' Additional optional fields: \code{name} and \code{opts}.
#' @param shortreads The type and location of the short reads as a named list
#' with the fields \code{type} ("illumina") and \code{dir}. Short reads are optional.
#' \code{NULL} can be provided.
#' @param datadir The data directory (See Note).
#' @param outdir The output directory (See Note).
#' @param reference The reference allele(s).
#' @param consensus \dQuote{\code{mapping}} for now (See Note).
#' @param threshold Threshold frequency for polymorphisms.
#' @param iterations Number of iterations of the mapIter step.
#' @param fullname Truncate allele names.
#' @param partSR Use shortreads in the mapInit step for getting polymorphic positions and
#' a first reference.
#' @param microsatellite FALSE Perform a second mapping of shortreads to the inferred reference in
#' mapInit. Set to TRUE if you know you have repeats like in microsatellites. Usually extends the
#' reference to a maximum length and enables a better mapping.
#' @param forceMapping FALSE set to TRUE if you want to force processing of bad shortreads, i.e.
#' when the distribution of coverage is bad. Aborts the program if maximum coverage > 75 \% quantile * 5.
#' @return A \code{\link[=DR2S_]{DR2S}} object.
#' @family DR2S mapper functions
#' @export
#' @importFrom S4Vectors metadata metadata<-
#' @examples
#' \dontrun{
#' x <- DR2Smap(
#'   sample = "ID12300527",
#'   locus = "DPB1",
#'   datadir = "/path/to/data",
#'   outdir = "/path/to/output",
#'   reference = "04:01:01:01",
#'   consensus = "mapping"
#'   ) %>%
#'   clear() %>%
#'   mapInit() %>%
#'   partitionLongReads() %>%
#'   splitLongReads() %>%
#'   extractFastq() %>%
#'   mapIter() %>%
#'   partitionShortReads() %>%
#'   mapFinal() %>%
#'   polish() %>%
#'   report(blockWidth = 60)
#' }
DR2Smap <- function(sample,
                    locus,
                    longreads = list(type = "pacbio", dir = "pacbio"),
                    shortreads = list(type = "illumina", dir = "illumina"),
                    datadir = ".",
                    outdir = "./output",
                    reference = NULL,
                    consensus = "mapping",
                    threshold = 0.20,
                    iterations = 1,
                    partSR = TRUE,
                    fullname = TRUE,
                    ...) {
  UseMethod("DR2Smap")
}


# mappers -----------------------------------------------------------------


#' Initial mapping step (long reads).
#'
#' Long reads (PacBio or Nanopore) of a heterozygous sample are mapped to a
#' single reference sequence.
#'
#' @param x A \code{\link[=DR2S_]{DR2S}} object.
#' @param opts Mapper options.
#' @param refname Which allele do we use as primary reference; the
#' \strong{ref}erence allele or the \strong{alt}ernate allele.
#' @param optsname Additional text to describe the options used.
#' @param pct Percentage of templates to subsample (before sorting the bam file).
#' @param threshold Threshold to call a variant.
#' @param minBaseQuality Minimum \sQuote{QUAL} value for each nucleotide in an
#' alignment.
#' @param minMapq Minimum \sQuote{MAPQ} value for an alignment to be included
#' in pileup.
#' @param maxDepth Maximum number of overlapping alignments considered for
#' each position in the pileup.
#' @param minNucleotideDepth Minimum count of each nucleotide at a given
#' position required for that nucleotide to appear in the result.
#' @param includeDeletions If \code{TRUE}, include deletions in pileup.
#' @param includeInsertions If \code{TRUE}, include insertions in pileup.
#' @param force If \code{TRUE}, overwrite existing bam file.
#' @param fullname If \code{TRUE}, use the full name in reference fasta.
#' @param plot Produce diagnostic plots.
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
#'   reference = "04:01:01:01",
#'   consensus = "mapping"
#'   ) %>%
#'   clear() %>%
#'   mapInit() %>%
#'   partitionLongReads() %>%
#'   mapIter() %>%
#'   partitionShortReads() %>%
#'   mapFinal() %>%
#'   polish() %>%
#'   report(blockWidth = 60)
#' }
mapInit <- function(x,
                    opts = list(),
                    refname = "ref",
                    optsname = "",
                    pct = 100,
                    threshold = 0.20,
                    minBaseQuality = 7,
                    minMapq = 0,
                    maxDepth = 1e4,
                    minNucleotideDepth = 3,
                    includeDeletions = FALSE,
                    includeInsertions = FALSE,
                    force = FALSE,
                    fullname = TRUE,
                    plot = TRUE,
                    ...) {
  UseMethod("mapInit")
}

#' Iterative mapping step (long reads).
#'
#' Partitioned long reads are mapped against the consensus sequence from MSA of the best 40
#' sequences for each haplotype. Indels are \strong{excluded} from pileup. New consensus sequences
#' for both alleles are inferred.
#'
#' @param x A \code{\link[=DR2S_]{DR2S}} object.
#' @param opts Mapper options.
#' @param pct Percentage of templates to subsample (before sorting the bam file).
#' @param minBaseQuality Minimum \sQuote{QUAL} value for each nucleotide in an
#' alignment.
#' @param minMapq Minimum \sQuote{MAPQ} value for an alignment to be included
#' in pileup.
#' @param maxDepth Maximum number of overlapping alignments considered for
#' each position in the pileup.
#' @param minNucleotideDepth Minimum count of each nucleotide at a given
#' position required for that nucleotide to appear in the result.
#' @param includeDeletions If \code{TRUE}, include deletions in pileup.
#' @param includeInsertions If \code{TRUE}, include insertions in pileup.
#' @param gapSuppressionRatio The ratio of base/gap above which gaps at
#' insertion position are excluded from from consensus calling.
#' @param force If \code{TRUE}, overwrite existing bam file.
#' @param fullname If \code{TRUE}, use the full name in reference fasta.
#' @param plot Plot diagnostics.
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
#'   reference = "04:01:01:01",
#'   consensus = "mapping"
#'   ) %>%
#'   clear() %>%
#'   mapInit() %>%
#'   partitionLongReads() %>%
#'   mapIter() %>%
#'   partitionShortReads() %>%
#'   mapFinal() %>%
#'   polish() %>%
#'   report(blockwidth = 60)
#' }
mapIter <- function(x,
                    opts = list(),
                    pct = 100,
                    minBaseQuality = 7,
                    minMapq = 0,
                    maxDepth = 1e4,
                    minNucleotideDepth = 3,
                    includeDeletions = TRUE,
                    includeInsertions = TRUE,
                    gapSuppressionRatio = 2/5,
                    force = FALSE,
                    fullname = TRUE,
                    plot = TRUE,
                    ...) {
  UseMethod("mapIter")
}

#' Final mapping step (long reads and short reads).
#'
#' Partitioned long reads and partitioned or unpartitioned short reads are mapped against the
#' consensus sequences inferred at \code{mapIter}. Indels are \strong{included}
#' in pileup.
#'
#' @param x A \code{\link[=DR2S_]{DR2S}} object.
#' @param opts Mapper options.
#' @param pct Percentage of templates to subsample (before sorting the bam file).
#' @param minBaseQuality Minimum \sQuote{QUAL} value for each nucleotide in an
#' alignment.
#' @param minMapq Minimum \sQuote{MAPQ} value for an alignment to be included
#' in pileup.
#' @param maxDepth Maximum number of overlapping alignments considered for
#' each position in the pileup.
#' @param minNucleotideDepth Minimum count of each nucleotide at a given
#' position required for that nucleotide to appear in the result.
#' @param includeDeletions If \code{TRUE}, include deletions in pileup.
#' @param includeInsertions If \code{TRUE}, include insertions in pileup.
#' from bam files for variants ar polymorphic positions.
#' @param force If \code{TRUE}, overwrite existing bam file.
#' @param fullname If \code{TRUE}, use the full name in reference fasta.
#' @param plot Plot diagnostics.
#' @param clip Clip homopolymeric read ends from shortreads.
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
#'   reference = "04:01:01:01",
#'   consensus = "mapping"
#'   ) %>%
#'   clear() %>%
#'   mapInit() %>%
#'   partitionLongReads() %>%
#'   mapIter() %>%
#'   partitionShortReads() %>%
#'   mapFinal() %>%
#'   polish() %>%
#'   report(blockWidth = 60)
#' }
mapFinal <- function(x,
                     opts = list(),
                     pct = 100,
                     minBaseQuality = 7,
                     minMapq = 0,
                     maxDepth = 1e5,
                     minNucleotideDepth = 3,
                     includeDeletions = TRUE,
                     includeInsertions = TRUE,
                     force = FALSE,
                     fullname = TRUE,
                     plot = TRUE,
                     clip = TRUE,
                     ...) {
  UseMethod("mapFinal")
}


# partitioning ------------------------------------------------------------


#' Partition mapped long reads into haplotypes
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
#'   reference = "04:01:01:01",
#'   consensus = "mapping"
#'   ) %>%
#'   clear() %>%
#'   mapInit() %>%
#'   partitionLongReads() %>%
#'   mapIter() %>%
#'   partitionShortReads() %>%
#'   mapFinal() %>%
#'   polish() %>%
#'   report(blockWidth = 60)
#' }
partitionLongReads <- function(x,
                               skipGapFreq = 2/3,
                               distAlleles = NULL,
                               noGapPartitioning = FALSE,
                               selectAllelesBy = "count",
                               ...
                               ) {
  UseMethod("partitionLongReads")
}

#' Assign short reads from mapInit to haplotypes
#'
#' @param x A \code{\link[=DR2S_]{DR2S}} object.
#' @param force force the creation of new fastq files if they already exist
#' @param ... Further arguments passed to methods.
#'
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
#'   reference = "04:01:01:01",
#'   consensus = "mapping"
#'   ) %>%
#'   clear() %>%
#'   mapInit() %>%
#'   partitionLongReads() %>%
#'   mapIter() %>%
#'   partitionShortReads() %>%
#'   mapFinal() %>%
#'   polish() %>%
#'   report(blockWidth = 60)
#' }
partitionShortReads <- function(x,
                                force = TRUE,
                                ...) {
  UseMethod("partitionShortReads")
}

# polishing and reporting ------------------------------------------------

#' Polish final haplotype sequences.
#'
#' @param x A \code{\link[=DR2S_]{DR2S}} object.
#' @param threshold When do we call a variant a variant.
#' @param cache Cache the updated \code{DR2S} object after assembling the haplotypes.
#' @param ... Additional arguments passed to methods.
#'
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
#'   reference = "04:01:01:01",
#'   consensus = "mapping"
#'   ) %>%
#'   clear() %>%
#'   mapInit() %>%
#'   partitionLongReads() %>%
#'   mapIter() %>%
#'   partitionShortReads() %>%
#'   mapFinal() %>%
#'   polish() %>%
#'   report(blockwidth = 60)
#' }
polish <- function(x, threshold = x$getThreshold(), cache = TRUE, ...) {
  UseMethod("polish")
}


#' Report the final haplotype sequences.
#'
#' @param x A \code{\link[=DR2S_]{DR2S}} object.
#' @param which Which mapping do we want to report. Will choose \code{mapFinal},
#' \code{mapIter} in this order if available.
#' @param blockWidth Maximum number of sequence letters per line in pairwise alignment.
#' @param ... Additional arguments passed to methods.
#' @details
#' \code{report} will create a \code{report} directory containing the consensus
#' sequences \code{map\{1,2,3\}.\{A,B\}.\{readtype\}.\{mapper\}.unchecked.fa},
#' an html alignment file \code{map\{1,2,3\}.aln.\{readtype\}.\{mapper\}.unchecked.html},
#' a pairwise alignment file \code{map\{1,2,3\}.\{A,B\}.\{readtype\}.\{mapper\}.unchecked.pair},
#' and a tab-separated file \code{problems.\{readtype\}.\{mapper\}.tsv} reporting
#' positions with potentially problematic results.
#'
#' To finalise the consensus sequences open the pairwise alignment file and save
#' it replacing \strong{unchecked} with \strong{checked} in the file name.
#' In the \strong{checked} file manual edits can be performed. After the pairwise
#' alignment is verified call \code{\link{reportCheckedConsensus}}.
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
#'   reference = "04:01:01:01",
#'   consensus = "mapping"
#'   ) %>%
#'   clear() %>%
#'   mapInit() %>%
#'   partitionLongReads() %>%
#'   mapIter() %>%
#'   partitionShortReads() %>%
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
#'   reference = "04:01:01:01",
#'   consensus = "mapping"
#'   ) %>%
#'   clear() %>%
#'   mapInit() %>%
#'   partitionLongReads() %>%
#'   mapIter() %>%
#'   partitionShortReads() %>%
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
#'   reference = "04:01:01:01",
#'   consensus = "mapping"
#'   ) %>%
#'   clear() %>%
#'   mapInit() %>%
#'   partitionLongReads() %>%
#'   mapIter() %>%
#'   partitionShortReads() %>%
#'   mapFinal() %>%
#'   polish() %>%
#'   report(blockWidth = 60)
#' }
cache <- function(x, ...) {
  UseMethod("cache")
}
