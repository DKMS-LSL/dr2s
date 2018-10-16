

# constructor -------------------------------------------------------------


#' Constructor for configuration to create a \code{\link[=DR2S_]{DR2S}} object.
#'
#' @usage createDR2SConf(sample, locus, longreads = list(type = "pacbio", 
#' dir = "pacbio"), shortreads = list(type = "illumina", dir = "illumina"), 
#' datadir = ".", outdir = "./output", reference = NULL,
#' threshold = 0.20, iterations = 1, microsatellite = FALSE, distAlleles = 2, 
#' filterScores = TRUE, partSR = TRUE, forceMapping = FALSE, fullname = TRUE, 
#' details = NULL, ...)
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
#' \code{OUTDIR/SAMPLE}
#'
#' @section Reference:
#'
#' References can be specified as allele codes or a path to a fasta file
#' containing the reference sequence. 
#'
#' @param sample A unique sample identifier used to locate the long and short
#' read FASTQ files.
#' @param locus The HLA or KIR locus.
#' @param longreads The type and location of the long reads as a named list
#' with the fields \code{type} ("pacbio" or "nanopore") and \code{dir}.
#' Additional optional fields: \code{name} and \code{opts}.
#' @param shortreads The type and location of the short reads as a named list
#' with the fields \code{type} ("illumina") and \code{dir} (optional).
#' @param datadir The data directory (See Note).
#' @param outdir The output directory (See Note).
#' @param reference The reference allele(s).
#' @param threshold Threshold frequency for polymorphisms.
#' @param iterations Number of iterations of the mapIter step.
#' @param microsatellite FALSE Perform a second mapping of shortreads to the 
#' inferred reference in mapInit. Set to TRUE if you know you have repeats like 
#' in microsatellites. Usually extends the reference to a maximum length and 
#' enables a better mapping.
#' @param partSR Use shortreads in the mapInit step for getting polymorphic 
#' positions and a first reference.
#' @param forceMapping FALSE set to TRUE if you want to force processing of bad
#' shortreads, i.e. when the distribution of coverage is bad. Aborts the program
#' if maximum coverage > 75 \% quantile * 5.
#' @param filterScores use only reads passing a strict filtering step. TODO
#' @param distAlleles Number of different alleles in the sample. Should be 2
#' for heterozygous samples, 1 for homozygous samples and > 2 for some KIR loci.
#' @param fullname Truncate allele names.
#' @param details Metadata of a sample. Will be written to the fasta header of 
#' the final sequences and stored in the config dump.
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
#'   reference = "04:01:01:01",
#'   consensus = "mapping"
#'   )) %>%
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
createDR2SConf <- function(sample,
                    locus,
                    longreads      = list(type = "pacbio", dir = "pacbio"),
                    shortreads     = list(type = "illumina", dir = "illumina"),
                    datadir        = ".",
                    outdir         = "./output",
                    reference      = NULL,
                    threshold      = 0.20,
                    iterations     = 1,
                    microsatellite = FALSE,
                    distAlleles    = 2,
                    filterScores   = TRUE,
                    partSR         = TRUE,
                    forceMapping   = FALSE,
                    fullname       = TRUE,
                    details        = NULL,
                    ...) {
  UseMethod("DR2SConf")
}

#' Initialise a \code{\link[=DR2S_]{DR2S}} mapper.
#' @param config A DR2S configuration object. Either created by 
#' \code{\link[=DR2S_]{createDR2SConfig}} or loaded by 
#' \code{\link[=DR2S_]{readDR2SConf}}.
#' @param createOutdir Create the outdir if not exists.
#' @return A \code{\link[=DR2S_]{DR2S}} object.
#' @export
InitDR2S <- function(config,
                     createOutdir = TRUE) {
  UseMethod("InitDR2S")
}
# mappers -----------------------------------------------------------------


#' Initial mapping step (long reads).
#'
#' Long reads (PacBio or Nanopore) of a heterozygous sample are mapped to a
#' single reference sequence.
#'
#' @param x A \code{\link[=DR2S_]{DR2S}} object.
#' @param opts Mapper options.
#' @param optsname Additional text to describe the options used.
#' @param partSR If \code{TRUE} use shortreads to infer the polymorphic 
#' positions for clustering.
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
#' @param microsatellite If \code{TRUE} remap the shortreads again to expand
#' the reference to a maximum length. Works better for larger insertions and
#' repeats like microsatellites. 
#' @param force If \code{TRUE}, overwrite existing bam file.
#' @param fullname If \code{TRUE}, use the full name in reference fasta.
#' @param filterScores Apply a harsh filtering on the reads. Filters for mapping
#' quality and removes all softclipping reads. Usually not necessary.
#' @param forceMapping if \code{FALSE} the program throws an error in case the
#' coverage of shortreads is too different at different parts of the sequence. 
#' @param topx the best x reads. 
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
                    optsname = "",
                    partSR = TRUE,
                    threshold = NULL,
                    minBaseQuality = 3,
                    minMapq = 50,
                    maxDepth = 1e4,
                    minNucleotideDepth = 3,
                    includeDeletions = TRUE,
                    includeInsertions = TRUE,
                    microsatellite = FALSE,
                    force = FALSE,
                    fullname = TRUE,
                    filterScores = TRUE,
                    forceMapping = FALSE,
                    topx = 0,
                    plot = TRUE) {
  UseMethod("mapInit")
}

#' Iterative mapping step (long reads).
#'
#' Partitioned long reads are mapped against the consensus sequence from MSA of
#' the best 40 sequences for each haplotype. Indels are \strong{excluded} from 
#' pileup. New consensus sequences for both alleles are inferred.
#'
#' @param x A \code{\link[=DR2S_]{DR2S}} object.
#' @param opts Mapper options.
#' @param iterations Number of iterations. How often are the clustered reads 
#' remapped to the updated reference. 
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
                    iterations = 1,
                    minBaseQuality = 3,
                    minMapq = 0,
                    maxDepth = 1e4,
                    minNucleotideDepth = 3,
                    includeDeletions = TRUE,
                    includeInsertions = TRUE,
                    gapSuppressionRatio = 2/5,
                    force = FALSE,
                    fullname = TRUE,
                    plot = TRUE) {
  UseMethod("mapIter")
}

#' Final mapping step (long reads and short reads).
#'
#' Partitioned long reads and partitioned or unpartitioned short reads are 
#' mapped against the consensus sequences inferred at \code{mapIter}. Indels 
#' are \strong{included} in pileup.
#'
#' @param x A \code{\link[=DR2S_]{DR2S}} object.
#' @param opts Mapper options.
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
                     minBaseQuality = 7,
                     minMapq = 0,
                     maxDepth = 1e5,
                     minNucleotideDepth = 3,
                     includeDeletions = TRUE,
                     includeInsertions = TRUE,
                     force = FALSE,
                     fullname = TRUE,
                     plot = TRUE,
                     clip = TRUE) {
  UseMethod("mapFinal")
}


# partitioning ------------------------------------------------------------


#' Partition mapped long reads into haplotypes
#'
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
#' with the most reads as the alleles.
#' @param plot Plot
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
                               threshold         = NULL,
                               skipGapFreq       = 2/3,
                               distAlleles       = NULL,
                               noGapPartitioning = FALSE,
                               selectAllelesBy   = "count",
                               plot              = TRUE,
                               ...) {
  UseMethod("partitionLongReads")
}

#' Assign short reads from mapInit to haplotypes
#'
#' @param x A \code{\link[=DR2S_]{DR2S}} object.
#' @param force force the creation of new fastq files if they already exist.
#' @param opts list with options passed to the mapper.
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
                                opts = list(),
                                force = TRUE, 
                                ...) {
  UseMethod("partitionShortReads")
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
polish <- function(x, threshold = x$getThreshold(), checkHpCount = TRUE, 
                   cache = TRUE) {
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
