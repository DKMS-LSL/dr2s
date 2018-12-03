
<!-- README.md is generated from README.Rmd. Please edit that file -->

# DR2S - dual redundant reference sequencing

An R package designed to facilitate generating reliable, full-length
phase-defined reference sequences for novel HLA and KIR alleles.

**Note, that this package is still maturing. There’s no guarantee yet
that the API or the under-the-hood workings of the package won’t change
substantially**

## Workflow

### Input and output

`DR2S` is designed to integrate longread HLA and KIR data (e.g., PacBio
or Oxford Nanopore sequences) and shortread shotgun data (Illumina). It
is also possible to run `DR2S` in “longread-only mode”, but don’t expect
reference-quality consensus sequences that way.

As input, we expect longread and, optionally, shortread FASTQ files to
be placed in separate subdirectories within a working directory and to
follow the naming convention `SAMPLEID_LOCUS_.*.fastq(.gz)?`.

`SAMPLEID` can be any arbitrary unique identification code and `LOCUS`
should be one of `A`, `B`, `C`, `DQB1`, `DRB1`, or `DPB1` for HLA, or
one of `2DL1`, `2DL2`, `2DL3`, `2DL4`, `2DL5A`, `2DL5B`, `2DP1`, `2DS1`,
`2DS2`, `2DS3`, `2DS4`, `2DS5`, `3DL1`, `3DL2`, `3DL3`, `3DP1`, or
`3DS1` for KIR.

All output is placed in a directory tree `output\LOCUS\SAMPLEID\`.

An example:

``` 
~/dr2s_data
   |
   +-- pacbio
   |     |-- ID12912701_DPB1_lbc23.fastq.gz
   |
   |
   +-- illumina
   |     |-- ID12912701_DPB1_S23_L001_R1_001.fastq.gz
   |     |-- ID12912701_DPB1_S23_L001_R2_001.fastq.gz
   |
   +-- output
         |
         +-- DPB1
               |
               +-- ID12912701
                 
```

### Usage

Once all input files are put in place a DR2S analysis is started with a
call to the functions `createDR2SConf()` and `InitDR2S()`:

``` r
## a minimal example:
x <- InitDR2S(createDR2SConf(
  sample = "ID12912701",
  locus = "DPB1",
  longreads = list(dir = "pacbio", type = "pacbio", mapper = "minimap"),
  shortreads = list(dir = "illumina", type = "illumina", mapper = "bwamem"),
  datadir = "~/dr2s_data",
  outdir = "~/dr2s_data/output"
))
```

This call generates an `R6` object of class `DR2S` that encapsulates all
data and methods for all subsequent analysis steps.

  - `sample`: A unique sample identifier. The FASTQs associated with a
    sample need need to be prefixed with this identifier.
  - `locus`: One of the allowed HLA and KIR loci above. If allele
    information for a sample is available it can be specified as, e.g.
    `DPB1*04:02:01:01`. In this case this allele will be used as a
    reference against which an initial mapping of the longreads is
    performed. If this information is not given a generic locus-specific
    reference is used. NOTE: generic references are not yet implemented
    for KIR.
  - `longreads`: The location, type, and mapper for longreads as a named
    list with the fields `dir`, `type` (“pacbio” or “nanopore”) and
    `mapper` (“bwamem” or “minimap”).
  - `shortreads`: (optional) The location, type, and mapper for
    shortreads as a named list with the fields `dir`, `type`
    (“illumina”) and `mapper` (“bwamem” or “minimap”).
  - `datadir`: The data directory (see above).
  - `outdir`: The output directory (see above).
  - `reference`: (optional) Path to a fasta file containing the
    reference sequence.
  - `details`: (optional) Named list of sample metadata. These data will
    be included in the fasta headers of the final sequences and stored
    in the config file.
  - `opts`: (optional) Named list of arguments to the DR2S pipeline
    steps. They will be stored in the config file.

An analysis proceeds in a number of steps that can be chained together
using the pipe `%>%`:

``` r
x %>% 
  mapInit() %>%
  partitionLongreads() %>%
  mapIter() %>%
  partitionShortreads() %>%
  mapFinal() %>%
  polish() %>%
  report()
```

Alternatively, the complete pipeline can be run in one go:

``` r
x$runPipeline()
```

The individual steps perform the following analyses:

  - `mapInit()`:
    1.  Map the shortreads against an initial reference. Construct a
        consensus.
    2.  Perform an optional second mapping of shortreads to the
        consensus sequence from the previous step. This expands the
        reference and may be necessary if there are extensive repeat
        structures like microsatellites in your sequence. Construct a
        consensus.
    3.  A final mapping of shortreads to the consensus from the previous
        step.
    4.  Map the longreads to the consensus from the previous step..
  - `partitionLongreads()`:
    1.  Construct a SNP matrix from the longreads at polymorphic
        positions in shortreads.
    2.  Perform hierarchical clustering,
    3.  Chimera detection and
    4.  Assign a haplotype and score to each read.
    5.  Pick longreads that best represent the allele haplotypes based
        on the *haplotype score* from the previous step.
    6.  Extract the chosen longreads from the alignment file and write
        the data as FASTQs into subdirectories for each haplotype.
  - `mapIter()`:
    1.  Construct a consensus sequence for each haplotype from the
        initial mapping while using only reads that are assigned to a
        haplotype.
    2.  Iteratively refine the long read consensus sequences by
        remapping long reads per haplotype to the latest haplotype
        consensus.
  - `partitionShortreads()`:
    1.  Partition the shortreads into haplotypes using the inferred
        longread haplotypes and initial mapping.
  - `mapFinal()`:
    1.  Map the short reads against the refined long read consensus
        sequences.
  - `polish()`: 1.Disambiguate polymorphic positions in the short read
    mapping. Check and report inconsistencies, insertions, deletions
    along the way. This should not be that much with okayish seuence
    data
  - `report()`:
    1.  Report the finalised shortread-based consensus sequences as
        FASTA files. Provide a tsv file with suspicious positions that
        may warrant manual inspection. Report the alignment of all
        haplotypes.

Throughout this process a number of diagnostic plots are produced and
placed in the output directory for later inspection.

While this process works remarkably well, there are situations where
alignment artefacts or plain bad luck may introduce errors in the final
consensus sequences. You should **never** accept the result as ground
truth without some manual and visual consistency checks\!

`DR2S` provides some facilities to aid checking and signing-off of
finalised consensus sequences.

A typical post-processing workflow may look as follows:

``` r
plot_diagnostic_alignment(x)
run_igv(x, 3000)
check_alignment_file(x)
refineAlignment(x, "A")
run_igv(x, "refine")
report_checked_consensus(x)
```

  - **plot\_diagnostic\_alignment:** Displays an alignment of three
    preliminary long-read-based consensus sequences and the final
    short-read-based consensus sequence for all alleles in your browser.
    This can be used to spot inconsistencies between the long-read and
    the short-read evidence.

  - **run\_igv:** Opens an instance of the IGV Genome Browser for each
    haplotype at a specified position (one for each allele) displaying
    both the long read and short read data for manual inspection.

  - **check\_alignment\_file :** Opens a pairwise or multiple alignment
    of the final consensus sequences in your text editor. Use this to
    perform any manual edits on the consensus sequences. Editor options
    are: “subl”, “gvim” and “gedit”. Defaults to systems standard
    editor.

  - **report\_checked\_consensus:** Export final consensus sequences
    from the edited pairwise or multiple alignment as FASTAs into a
    separate subdirectory `./checked` in the output directory.

`DR2S` creates bash scripts for the convenient access to important
postprocessing functions:

  - **run\_checkConsensus.sh** Runs the `check_alignment_file` command.
  - **runIGV\_mapInit.sh** Opens an IGV instance of the initial mapping.
  - **runIGV\_mapIter.sh** Opens an IGV instance of the results of the
    mapIter step.
  - **runIGV\_mapFinal.sh** Opens an IGV instance of the final mapping
  - **run\_remap\[X\].sh** Remap the reads of a haplotype to the manuall
    curated sequence to look if it is finally correct. This command is
    available for all found haplotypes
  - **run\_reportCheckedConsensus.sh** report the manually checked
    consensus and state that its finished and can be used.

### Only longreads

If you want to find allele sequences only based on longreads you need to
skip the partitioning of shortreads. The longread partitioning in other
steps is automatically done if `partSR` is set to `FALSE` and no
shortreads are found in the datadir.

``` r
x %>% 
  mapInit() %>%
  partition_haplotypes() %>%
  split_reads_by_haplotype() %>%
  extract_fastq() %>%
  mapIter() %>%
  mapFinal() %>%
  polish() %>%
  report()
```

## Installation

This package is only available via gitlab now. It depends on a local
installation of `samtools`, `bwa` (\>= 0.7.11), python and the pysam
library, and a C++11 compliant compiler. An alternative mapper with
better results is `minimap2`, a successor of `bwa` For installation over
git via SSH, you need to have installed libssh2-1 and libssh2-1-dev
prior to install git2r/devtools or reinstall the package afterwards. You
do also need to have set up your ssh key in gitlab. Otherwise you need
to give your gitlab password in plain.

``` r
install.packages("devtools")  # if not already installed
library("devtools")
install_git("git@srvddgit01.labor.local:rlib/IPDdata.git")
install_git("git@srvddgit01.labor.local:rlib/hlatools.git")
install_git("git@srvddgit01.labor.local:sklas/gDR2S.git")
```
