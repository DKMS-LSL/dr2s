

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
#' containing the reference sequence. Specifying both alleles as a character
#' vector is possible but only takes effect if \code{consensus} is set to
#' \dQuote{\code{mapping}}. If \code{reference = NULL} a global generic reference
#' for a given locus will be used.
#'
#' @section Consensus:
#' Two methods to generate initial consensus sequences from long reads are
#' available: \dQuote{\bold{\code{multialign}}}: After partitioning long reads into
#' haplotypes up to 30 "best" reads from each haplotype group are selected for
#' multiple sequence alignment and a preliminary consensus sequences construction.
#' All reads from each haplotype group are subsequently mapped against their
#' respective preliminary consensus sequences and these mappings are used for
#' constructing final long read consensus sequences.
#' \dQuote{\bold{\code{mapping}}}: Reads from each haplotype group are mapped
#' against the orignal reference sequence and these mappings are used for
#' constructing final long read consensus sequences. Note, that the second
#' approach is prone to introducing alignment artifacts if there is a lot of
#' mismatch between the reference used and the actual allele sequences.
#'
#' @param sample A unique sample identifier used to locate the long and short
#' read FASTQ files.
#' @param locus The HLA locus ("A", "B", "C", "DPB1", "DQB1", "DRB1")
#' @param longreads The type and location of the long reads as a named list
#' with the fields \code{type} ("pacbio" or "nanopore") and \code{dir}.
#' Additional optional fields: \code{name} and \code{opts}.
#' @param shortreads The type and location of the short reads as a named list
#' with the fields \code{type} ("illumina") and \code{dir}. Short reads are optional.
#' \code{NULL} can be provided.
#' @param datadir The data directory (See Note).
#' @param outdir The output directory (See Note).
#' @param reference (Optional) The reference allele(s). One or both alleles
#' specified as a character vector (See Note).
#' @param consensus How will the inital allele consensus sequences be generated
#' from the long reads. One of \dQuote{\code{multialign}} or \dQuote{\code{mapping}}
#' (See Note).
#' @param threshold Threshold frequency for polymorphisms.
#' @param fullname Truncate allele names.
#'
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
#'   consensus = "multialign"
#'   ) %>%
#'   clear() %>%
#'   map0() %>%
#'   partition_haplotypes() %>%
#'   split_reads_by_haplotype() %>%
#'   extract_fastq() %>%
#'   map1() %>%
#'   map2() %>%
#'   map3() %>%
#'   cache() %>%
#'   polish() %>%
#'   report(block_width = 60)
#' }
DR2Smap <- function(sample,
                    locus,
                    longreads = list(type = "pacbio", dir = "pacbio"),
                    shortreads = list(type = "illumina", dir = "illumina"),
                    datadir = ".",
                    outdir = "./output",
                    reference = NULL,
                    consensus = "multialign",
                    threshold = 0.20,
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
#' @param min_base_quality Minimum \sQuote{QUAL} value for each nucleotide in an
#' alignment.
#' @param min_mapq Minimum \sQuote{MAPQ} value for an alignment to be included
#' in pileup.
#' @param max_depth Maximum number of overlapping alignments considered for
#' each position in the pileup.
#' @param min_nucleotide_depth Minimum count of each nucleotide at a given
#' position required for that nucleotide to appear in the result.
#' @param include_deletions If \code{TRUE}, include deletions in pileup.
#' @param include_insertions If \code{TRUE}, include insertions in pileup.
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
#'   consensus = "multialign"
#'   ) %>%
#'   clear() %>%
#'   map0() %>%
#'   partition_haplotypes() %>%
#'   split_reads_by_haplotype(limitA = 0.5, limitB = 0.25) %>%
#'   extract_fastq() %>%
#'   map1() %>%
#'   map2() %>%
#'   map3() %>%
#'   cache() %>%
#'   polish() %>%
#'   report(block_width = 60)
#' }
map0 <- function(x,
                 opts = list(),
                 refname = "ref",
                 optsname = "",
                 pct = 100,
                 threshold = 0.20,
                 min_base_quality = 7,
                 min_mapq = 0,
                 max_depth = 1e4,
                 min_nucleotide_depth = 3,
                 include_deletions = FALSE,
                 include_insertions = FALSE,
                 force = FALSE,
                 fullname = TRUE,
                 plot = TRUE,
                 ...) {
  UseMethod("map0")
}

#' Second mapping step (long reads).
#'
#' Partitioned long reads are mapped against the reference and alternate
#' sequences. Indels are \strong{excluded} from pileup. New consensus sequences
#' for both alleles are inferred.
#'
#' @param x A \code{\link[=DR2S_]{DR2S}} object.
#' @param opts Mapper options.
#' @param pct Percentage of templates to subsample (before sorting the bam file).
#' @param min_base_quality Minimum \sQuote{QUAL} value for each nucleotide in an
#' alignment.
#' @param min_mapq Minimum \sQuote{MAPQ} value for an alignment to be included
#' in pileup.
#' @param max_depth Maximum number of overlapping alignments considered for
#' each position in the pileup.
#' @param min_nucleotide_depth Minimum count of each nucleotide at a given
#' position required for that nucleotide to appear in the result.
#' @param include_deletions If \code{TRUE}, include deletions in pileup.
#' @param include_insertions If \code{TRUE}, include insertions in pileup.
#' @param pruning_cutoff TODO
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
#'   consensus = "multialign"
#'   ) %>%
#'   clear() %>%
#'   map0() %>%
#'   partition_haplotypes() %>%
#'   split_reads_by_haplotype(limitA = 0.5, limitB = 0.25) %>%
#'   extract_fastq() %>%
#'   map1() %>%
#'   map2() %>%
#'   map3() %>%
#'   cache() %>%
#'   polish() %>%
#'   report(block_width = 60)
#' }
map1 <- function(x,
                 opts = list(),
                 pct = 100,
                 min_base_quality = 7,
                 min_mapq = 0,
                 max_depth = 1e4,
                 min_nucleotide_depth = 3,
                 include_deletions = TRUE,
                 include_insertions = TRUE,
                 pruning_cutoff = 0.75,
                 force = FALSE,
                 fullname = TRUE,
                 plot = TRUE,
                 ...) {
  UseMethod("map1")
}

#' Third mapping step (long reads).
#'
#' Partitioned long reads are mapped against the consensus sequences inferred
#' at \code{map1}. Indels are \strong{included} in pileup. New consensus
#' sequences for both alleles are inferred.
#'
#' @param x A \code{\link[=DR2S_]{DR2S}} object.
#' @param opts Mapper options.
#' @param pct Percentage of templates to subsample (before sorting the bam file).
#' @param min_base_quality Minimum \sQuote{QUAL} value for each nucleotide in an
#' alignment.
#' @param min_mapq Minimum \sQuote{MAPQ} value for an alignment to be included
#' in pileup.
#' @param max_depth Maximum number of overlapping alignments considered for
#' each position in the pileup.
#' @param min_nucleotide_depth Minimum count of each nucleotide at a given
#' position required for that nucleotide to appear in the result.
#' @param include_deletions If \code{TRUE}, include deletions in pileup.
#' @param include_insertions If \code{TRUE}, include insertions in pileup.
#' @param pruning_cutoff TODO
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
#'   consensus = "multialign"
#'   ) %>%
#'   clear() %>%
#'   map0() %>%
#'   partition_haplotypes() %>%
#'   split_reads_by_haplotype(limitA = 0.5, limitB = 0.25) %>%
#'   extract_fastq() %>%
#'   map1() %>%
#'   map2() %>%
#'   map3() %>%
#'   cache() %>%
#'   polish() %>%
#'   report(block_width = 60)
#' }
map2 <- function(x,
                 opts = list(),
                 pct = 100,
                 min_base_quality = 7,
                 min_mapq = 0,
                 max_depth = 1e4,
                 min_nucleotide_depth = 3,
                 include_deletions = TRUE,
                 include_insertions = TRUE,
                 pruning_cutoff = 0.75,
                 force = FALSE,
                 fullname = TRUE,
                 plot = TRUE,
                 ...) {
  UseMethod("map2")
}

#' Fourth mapping step (long reads and short reads).
#'
#' Partitioned long reads and unpartitioned short reads are mapped against the
#' consensus sequences inferred at \code{map2}. Indels are \strong{included}
#' in pileup.
#'
#' @param x A \code{\link[=DR2S_]{DR2S}} object.
#' @param opts Mapper options.
#' @param pct Percentage of templates to subsample (before sorting the bam file).
#' @param min_base_quality Minimum \sQuote{QUAL} value for each nucleotide in an
#' alignment.
#' @param min_mapq Minimum \sQuote{MAPQ} value for an alignment to be included
#' in pileup.
#' @param max_depth Maximum number of overlapping alignments considered for
#' each position in the pileup.
#' @param min_nucleotide_depth Minimum count of each nucleotide at a given
#' position required for that nucleotide to appear in the result.
#' @param include_deletions If \code{TRUE}, include deletions in pileup.
#' @param include_insertions If \code{TRUE}, include insertions in pileup.
#' @param include_read_ids If \code{TRUE}, read identifieres are extracted
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
#'   consensus = "multialign"
#'   ) %>%
#'   clear() %>%
#'   map0() %>%
#'   partition_haplotypes() %>%
#'   split_reads_by_haplotype(limitA = 0.5, limitB = 0.25) %>%
#'   extract_fastq() %>%
#'   map1() %>%
#'   map2() %>%
#'   map3() %>%
#'   cache() %>%
#'   polish() %>%
#'   report(block_width = 60)
#' }
map3 <- function(x,
                 opts = list(),
                 pct = 100,
                 min_base_quality = 7,
                 min_mapq = 0,
                 max_depth = 1e5,
                 min_nucleotide_depth = 3,
                 include_deletions = TRUE,
                 include_insertions = TRUE,
                 include_read_ids = TRUE,
                 force = FALSE,
                 fullname = TRUE,
                 plot = TRUE,
                 clip = TRUE,
                 ...) {
  UseMethod("map3")
}


# partitioning ------------------------------------------------------------


#' Partition mapped long reads into haplotypes
#'
#' @param x A \code{\link[=DR2S_]{DR2S}} object.
#' @param max_depth Maximum number of alignments considered for partitioning.
#' @param shuffle Randomly shuffle polymorphic positions.
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
#'   consensus = "multialign"
#'   ) %>%
#'   clear() %>%
#'   map0() %>%
#'   partition_haplotypes() %>%
#'   split_reads_by_haplotype(limitA = 0.5, limitB = 0.25) %>%
#'   extract_fastq() %>%
#'   map1() %>%
#'   map2() %>%
#'   map3() %>%
#'   cache() %>%
#'   polish() %>%
#'   report(block_width = 60)
#' }
partition_haplotypes <- function(x,
                                 max_depth = 1e4,
                                 plot = TRUE,
                                 ...) {
  UseMethod("partition_haplotypes")
}

#' Partition mapped long reads into haplotypes
#'
#' @param x A \code{\link[=DR2S_]{DR2S}} object.
#' @param limits List of Cutoff values for haplotypes.
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
#'   consensus = "multialign"
#'   ) %>%
#'   clear() %>%
#'   map0() %>%
#'   partition_haplotypes() %>%
#'   split_reads_by_haplotype() %>%
#'   extract_fastq() %>%
#'   map1() %>%
#'   map2() %>%
#'   map3() %>%
#'   cache() %>%
#'   polish() %>%
#'   report(block_width = 60)
#' }
split_reads_by_haplotype <- function(x,
                                     limits,
                                     ...) {
  UseMethod("split_reads_by_haplotype")
}

#' Extract FASTQs for partitioned reads
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
#'   consensus = "multialign"
#'   ) %>%
#'   clear() %>%
#'   map0() %>%
#'   partition_haplotypes() %>%
#'   split_reads_by_haplotype(limitA = 0.5, limitB = 0.25) %>%
#'   extract_fastq() %>%
#'   map1() %>%
#'   map2() %>%
#'   map3() %>%
#'   cache() %>%
#'   polish() %>%
#'   report(block_width = 60)
#' }
extract_fastq <- function(x, ...) {
  UseMethod("extract_fastq")
}

# polishing and reporting ------------------------------------------------

#' Polish final haplotype sequences.
#'
#' @param x A \code{\link[=DR2S_]{DR2S}} object.
#' @param threshold When do we call a variant a variant.
#' @param lower_limit When do we report a variant as uncertain.
#' @param cache Cache the updated \code{DR2S} object after assembling the haplotypes.
#' @param ... Additional arguments passed to methods.
#'
#' @return A \code{\link[=DR2S_]{DR2S}} object with the \code{consensus} field
#' populated.
#' @family DR2S mapper functions
#' @export
#' @examples
#' ## Run the full DR2S pipeline
#' \dontrun{
#' x <- DR2Smap(
#'   sample = "ID12300527",
#'   locus = "DPB1",
#'   datadir = "/path/to/data",
#'   outdir = "/path/to/output",
#'   reference = "04:01:01:01",
#'   consensus = "multialign"
#'   ) %>%
#'   clear() %>%
#'   map0() %>%
#'   partition_haplotypes() %>%
#'   split_reads_by_haplotype(limitA = 0.5, limitB = 0.25) %>%
#'   extract_fastq() %>%
#'   map1() %>%
#'   map2() %>%
#'   map3() %>%
#'   cache() %>%
#'   polish() %>%
#'   report(block_width = 60)
#' }
polish <- function(x, threshold = x$getThreshold(), lower_limit = 0.6, cache = TRUE, ...) {
  UseMethod("polish")
}


#' Report the final haplotype sequences.
#'
#' @param x A \code{\link[=DR2S_]{DR2S}} object.
#' @param which Which mapping do we want to report. Will choose \code{map3},
#' \code{map2}, \code{map1} in this order if available.
#' @param block_width Maximum number of sequence letters per line in pairwise alignment.
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
#' alignment is verified call \code{\link{report_checked_consensus}}.
#'
#' @return Side effect
#' @family DR2S mapper functions
#' @export
#' @examples
#' ## Run the full DR2S pipeline
#' \dontrun{
#' x <- DR2Smap(
#'   sample = "ID12300527",
#'   locus = "DPB1",
#'   datadir = "/path/to/data",
#'   outdir = "/path/to/output",
#'   reference = "04:01:01:01",
#'   consensus = "multialign"
#'   ) %>%
#'   clear() %>%
#'   map0() %>%
#'   partition_haplotypes() %>%
#'   split_reads_by_haplotype() %>%
#'   extract_fastq() %>%
#'   map1() %>%
#'   map2() %>%
#'   map3() %>%
#'   cache() %>%
#'   polish() %>%
#'   report(block_width = 60)
#' }
report <- function(x, which, block_width = 80, ...) {
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
#'   consensus = "multialign"
#'   ) %>%
#'   clear() %>%
#'   map0() %>%
#'   partition_haplotypes() %>%
#'   split_reads_by_haplotype(limitA = 0.5, limitB = 0.25) %>%
#'   extract_fastq() %>%
#'   map1() %>%
#'   map2() %>%
#'   map3() %>%
#'   cache() %>%
#'   polish() %>%
#'   report(block_width = 60)
#'   }
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
#'   consensus = "multialign"
#'   ) %>%
#'   clear() %>%
#'   map0() %>%
#'   partition_haplotypes() %>%
#'   split_reads_by_haplotype(limitA = 0.5, limitB = 0.25) %>%
#'   extract_fastq() %>%
#'   map1() %>%
#'   map2() %>%
#'   map3() %>%
#'   cache() %>%
#'   polish() %>%
#'   report(block_width = 60)
#'   }
cache <- function(x, ...) {
  UseMethod("cache")
}
