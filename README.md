
<!-- README.md is generated from README.Rmd. Please edit that file -->

# DR2S - dual redundant reference sequencing

An R package designed to facilitate generating reliable, full-length
phase-defined reference sequences for novel HLA and KIR alleles.

## Installation

This package is only available on GitHub for now. It depends on a local
installation of `samtools`, `bwa` (&gt;= 0.7.11) and `minimap2`. Some
used R packages have additional system dependencies. Bash:

``` bash
sudo apt-get update
sudo apt-get install \
  build-essentials \ # for samtools, bwa 
  gcc \              # for samtools, bwa 
  autoconf \         # for samtools, bwa 
  libxml2-dev \
  libssl2-dev \
  libz-dev \
  libbz2-dev \
  liblzma-dev \
  libncurses5-dev
```

R:

``` r
install.packages("devtools")  # if not already installed
devtools::install_github("DKMS-LSL/DR2S")
```

A Docker image is also provided for convenience at docker hub. This can
be loaded and used with the following command:

``` bash
docker pull dkmslsl/testrepo
docker run --rm -p 8788:8787 -e PASSWORD=<yourpassword> -it  dkmslsl/testrepo
```

An rstudio server session can be accessed in your browser at
localhost:8788 (The default port 8787 is mapped to the host port 8788 in
case you have, like me, rstudio server already running). Login to
rstudio using the username “rstudio” and your set password. A detailed
introduction to use DR2S with example data can be found at the Vignette
“DR2S”.

## Workflow

### Input and output

`DR2S` is designed to integrate long-read HLA and KIR data (e.g., PacBio
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

All output is placed in a directory tree under the configured output
directory `LOCUS\SAMPLEID\`.

An example:

    ~/dr2s_data
       |
       +-- pacbio
       |     |-- ID123_DPB1_lbc23.fastq.gz
       |
       |
       +-- illumina
       |     |-- ID123_DPB1_S23_L001_R1_001.fastq.gz
       |     |-- ID123_DPB1_S23_L001_R2_001.fastq.gz
       |
       +-- output
             |
             +-- DQB1
                   |
                   +-- ID123
                          | --...
                          | --...
                          | --...

### Usage

Once all input files are put in place a DR2S analysis is started with a
call to the functions `createDR2SConf()` or `readDR2SConf()` and
`InitDR2S()`:

``` r
## a minimal example:
x <- InitDR2S(
  createDR2SConf(
    sample = "ID123",
    locus = "DPB1",
    longreads = list(dir = "Sequel", type = "pacbio"),
    shortreads = list(dir = "Illumina", type = "illumina"),
    datadir = "~/dr2s_data",
    outdir = "~/dr2s_data/output"
))
```

#### Arguments

-   `sample`: A unique sample identifier. The FASTQs associated with a
    sample need need to be prefixed with this identifier.
-   `locus`: One of the allowed HLA and KIR loci above. If allele
    information for a sample is available it can be specified as,
    e.g. `DPB1*04:02:01:01`. In this case this allele will be used as a
    reference against which an initial mapping of the longreads is
    performed. If this information is not given a generic locus-specific
    reference is used. NOTE: generic references are not yet implemented
    for KIR.
-   `longreads`: The location, type, and mapper for longreads as a named
    list with the fields `dir`, `type` (“pacbio” or “nanopore”) and
    `mapper` (“bwamem” or “minimap”).
-   `shortreads`: (optional) The location, type, and mapper for
    shortreads as a named list with the fields `dir`, `type`
    (“illumina”) and `mapper` (“bwamem” or “minimap”).
-   `datadir`: The data directory (see above).
-   `outdir`: The output directory (see above).
-   `reference`: (optional) Path to a fasta file containing the
    reference sequence.
-   `details`: (optional) Named list of sample metadata. These data will
    be included in the fasta headers of the final sequences and stored
    in the config file.
-   `opts`: (optional) Named list of arguments to the DR2S pipeline
    steps. They will be stored in the config file. See below for a
    detailed descriptions of options that control the DR2S pipeline.

This call generates an `R6` object of class `DR2S` that encapsulates all
data and methods for all subsequent analysis steps.

An alternative approach is to use yaml or json config files. An example
config file is provided in the toy example at `URL to github project`:

``` r
configFile <- "~/dr2s_data/dr2s_config.yaml"
x <- InitDR2S(
  readDR2SConf(configFile)
)
```

#### Pipeline

An analysis proceeds in a number of steps that can be chained together
using the pipe `%>%`:

``` r
x %>% 
  mapInit() %>%
  partitionLongreads() %>%
  mapIter() %>%
  partitionShortreads() %>%
  mapFinal() %>%
  report() %>% 
  cache()
```

Alternatively, the complete pipeline can be run in one go. There are two
available pipelines, one for the standard run (SR), and another for only
long-read data (LR). Which pipeline to run can be configured in the
config file.

``` r
x$runPipeline()
```

The individual steps perform the following analyses:

-   `mapInit()`:
    1.  Map the shortreads (SR) against an initial reference. Construct
        a tentative consensus.
    2.  Perform a second SR mapping to the consensus sequence from the
        previous step. This expands the reference and may be necessary
        if there are extensive repeat structures like microsatellites in
        your sequence. Construct a tentative consensus.
    3.  A final SR mapping to the consensus from the previous step.
    4.  A longread (LR) mapping to the consensus from the previous step.
-   `partitionLongreads()`:
    1.  Infer polymorphic positions from the SR mapping performed in the
        previous step.
    2.  Construct a SNP matrix from the LR mapping at the polymorphic
        positions infer from the SR mapping.
    3.  Perform hierarchical clustering the LR SNPs.
    4.  Attempt to detect chimeric reads.
    5.  Assign a haplotype and a *haplotype\_membership\_coefficient* to
        each read.
    6.  Pick longreads that best represent the allele haplotypes based
        on the *haplotype\_membership\_coefficient* from the previous
        step.
    7.  Extract the selected longreads from the alignment file and
        output the data as FASTQs into subdirectories for each
        haplotype.
-   `mapIter()`:
    1.  Construct a consensus sequence for each haplotype from the
        initial LR mapping whilst using only reads that have been
        assigned to that haplotype.
    2.  Iteratively refine the longread consensus sequences by remapping
        longreads per haplotype to the latest haplotype consensus.
-   `partitionShortreads()`:
    1.  Partition the shortreads into haplotypes using the inferred
        longread haplotypes and the initial SR mapping.
-   `mapFinal()`:
    1.  Map the shortreads against the refined longread consensus
        sequences.
-   `report()`:
    1.  Report the finalised shortread-based consensus sequences as
        FASTA files. Provide a tsv file with suspicious positions that
        may warrant manual inspection. Report the alignment of all
        haplotypes.

Throughout this process a number of diagnostic plots are produced and
placed in the output directory for later inspection.

While this process works remarkably well, there are situations where
alignment artefacts or plain bad luck may introduce errors in the final
consensus sequences. You should **never** accept the result as ground
truth without some manual and visual consistency checks!

`DR2S` provides some facilities to aid checking and signing-off of
finalised consensus sequences.

#### Postprocessing

A typical post-processing workflow may look as follows:

``` r
run_igv(x, 3000)
check_alignment_file(x)
refineAlignment(x, "A")
run_igv(x, "refine")
report_checked_consensus(x)
```

-   **run\_igv:** Opens an instance of the IGV Genome Browser for each
    haplotype at a specified position (one for each allele) displaying
    both the long read and short read data for manual inspection.

-   **check\_alignment\_file :** Opens a pairwise or multiple alignment
    of the final consensus sequences in your text editor. Use this to
    perform any manual edits on the consensus sequences. Editor options
    are: “subl”, “gvim” and “gedit”. Defaults to systems standard
    editor.

-   **report\_checked\_consensus:** Export final consensus sequences
    from the edited pairwise or multiple alignment as FASTAs into a
    separate subdirectory `./checked` in the output directory.

`DR2S` creates bash scripts for the convenient access to important
postprocessing functions:

-   **run\_checkConsensus.sh** Runs the `check_alignment_file` command.
-   **runIGV\_mapInit.sh** Opens an IGV instance of the initial mapping.
-   **runIGV\_mapIter.sh** Opens an IGV instance of the results of the
    mapIter step.
-   **runIGV\_mapFinal.sh** Opens an IGV instance of the final mapping
-   **run\_remap\[X\].sh** Remap the reads of a haplotype to the manuall
    curated sequence to look if it is finally correct. This command is
    available for all found haplotypes
-   **run\_reportCheckedConsensus.sh** report the manually checked
    consensus and state that its finished and can be used.

#### Pipeline control options

A more fine-grained control of the DR2S pipeline is available via a
`json`-based config file. This config file can be created externally and
read in using the `readDR2Sconf()` function. Alternatively the config
file is created by the `createDR2Sconf()` function. Pipeline control
options are set using the `opts` argument in `createDR2Sconf()`, e.g.:

``` r
conf <- createDR2SConf(
  sample = "ID12912701",
  locus = "A*01:01:01:01",
  longreads = list(dir = "pacbio", type = "pacbio", mapper = "minimap"),
  shortreads = list(dir = "illumina", type = "illumina", mapper = "bwamem"),
  datadir = "~/dr2s_data",
  outdir = "~/dr2s_data/output",
  opts = list(
    mapInit = list(topx = "auto",
                   createIgv = FALSE),
    partitionLongreads = list(threshold = 1/5,
                              noGapPartitioning = TRUE,
                              selectCorrelatedPositions = TRUE,
                              selectAllelesBy = "distance"),
    mapIter  = list(iterations = 2),
    mapFinal = list(createIgv = FALSE),
    report   = list(createIgv = FALSE)
  ))
```

The complete set of available options and their defaults are:

``` r
  ##
  ## mapInit() defaults ####
  ##
  mapInit = list(
    ## <includeDeletions>: include deletions in pileup.
    includeDeletions = TRUE,
    ## <includeInsertions>: include insertions in pileup.
    includeInsertions = TRUE,
    ## <callInsertionThreshold>: if <includeInsertions == TRUE>, an insertion
    ## needs to be at frequency <callInsertionThreshold> for it to be included
    ## in the pileup.
    callInsertionThreshold = 1/5,
    ## <microsatellite>: if <pipeline == "SR">, perform a second mapping of
    ## shortreads to the inferred reference. Set to TRUE if you suspect
    ## microsatellites or repetitive regions in your sequence. This extends
    ## the reference to a maximum length and enables a better mapping.
    microsatellite = FALSE,
    ## <forceMapping>: set to TRUE if you want to force processing of "bad"
    ## shortreads when the distribution of coverage is heavily unequal.
    ## Aborts the program if maximum coverage > 75 % quantile * 5.
    forceMapping = FALSE,
    ## <minMapq>: don't filter longreads for mapping quality unless specified.
    ## NOTE: for shortreads <minMapq = 50> is hardcoded.
    minMapq = 0,
    ## <topx>: pick the x top-scoring reads. Set to an integer value to pick
    ## a fixed number of reads. Set to "auto" to use a dynamically determined
    ## number of reads to be selected.
    topx = FALSE,
    ## <pickiness>: if <topx == "auto">: <pickiness < 1>: bias towards higher
    ## scores/less reads; <pickiness > 1>: bias towards lower scores/more reads
    pickiness = 1,
    ## <increasePickiness>: if <topx == "auto">: increase pickiness for the
    ## second iteration of LR mapping
    increasePickiness = 1,
    ## <lowerLimit>: if <topx == "auto"> or <topx > 0>: the  minimum number
    ## of reads to pick if available.
    lowerLimit = 200,
    ## <updateBackgroundModel>: estimate the indel noise in a pileup and use
    ## this information to update the background model for PWM scoring
    updateBackgroundModel = FALSE,
    ## <createIgv>: subsample bam files for visualisation with IgvJs in the
    ## DR2S shiny app.
    createIgv = TRUE,
    ## <plot>: generate diagnostic plots.
    plot = TRUE
  )
  ##
  ## partitionLongreads() defaults ####
  ##
  partitionLongreads = list(
    ## Threshold to call a polymorphic position. A minority nucleotide frequency
    ## below this threshold is considered noise rather than a valid polymorphism.
    threshold = 1/5,
    ## The expected number of distinct alleles in the sample. This should be 2
    ## for heterozygous samples, 1 for homozygous samples may be >2 for some
    ## KIR loci.
    distAlleles = 2,
    ## The minumum frequency of the gap character required to call a gap position.
    skipGapFreq = 2/3,
    ## Don't partition based on gaps. Useful for samples with only few SNPs but
    ## with homopolymers. The falsely called gaps could mask the real variation.
    ## Set to override global default.
    noGapPartitioning = TRUE,
    ## Correlate polymorphic positions and cluster based on the absolute
    ## correlation coefficient. Extract positions from the cluster with the
    ## higher absolute mean correlation coefficient. This gets rid of positions
    ## that are not well distributed across the two alleles.
    selectCorrelatedPositions = FALSE,
    ## if <selectCorrelatedPositions> == TRUE, use <measureOfAssociation>
    ## ("cramer.V" or "spearman") to determine linkage between all polymorphic
    ## positions.
    measureOfAssociation = "cramer.V",
    ## We perform an equivalence test on clusters of polymorhic positions:
    ## Calculate the lower 1-sigma bound of the high-association cluster i.
    ## Calculate the upper 1-sigma bound of the low-association cluster j.
    ## Reject the clusters, if this bounds overlap by more than <proportionOfOverlap>
    ## of the average distance (dij) between clusters.
    proportionOfOverlap = 1/3,
    ## By how much do we expect 2 clusters to minimally differ in mean Cramér's V.
    ## BIC-informed model-based clustering tends to split rather than lump
    ## and this is a heuristical attempt to forestall this.
    minimumExpectedDifference = 0.06,
    ## If more than <distAlleles> clusters are found select clusters based on:
    ## (1) "distance": The hamming distance of the resulting variant consensus
    ## sequences or (2) "count": Take the clusters with the most reads as the
    ## true alleles.
    selectAllelesBy = "distance",
    ## Minimum size of an allele cluster
    minClusterSize = 20,
    ## When selecting reads from allele clusters using a dynamic threshold:
    ## pickiness < 1: bias towards higher scores/less reads
    ## pickiness > 1: bias towards lower scores/more reads
    pickiness = 1,
    ## When selecting reads from allele clusters the minimum number of
    ## reads to pick if available.
    lowerLimit = 40,
    ## Generate diagnostic plots.
    plot = TRUE
  )
  ##
  ## mapIter() defaults ####
  ##
  mapIter = list(
    ## Number of <mapIter> iterations. How often are the
    ## clustered reads remapped to updated reference sequences.
    iterations = 1,
    ## Minimum occupancy (1 - fraction of gap) below which
    ## bases at insertion position are excluded from from consensus calling.
    columnOccupancy = 2/5,
   ## an insertion needs to be at frequency <callInsertionThreshold> for it
    ## to be included in the pileup.
    callInsertionThreshold = 1/5,
    ## Generate diagnostic plots.
    plot = TRUE
  )
  ##
  ## mapFinal() defaults ####
  ##
  mapFinal = list(
    ## include deletions in pileup.
    includeDeletions = TRUE,
    ## include insertions in pileup.
    includeInsertions = TRUE,
    ## an insertion needs to be at frequency <callInsertionThreshold> for it
    ## to be included in the pileup.
    callInsertionThreshold = 1/5,
    ## (for shortreads only) trim softclips and polymorphic ends of reads before
    ## the final mapping
    trimPolymorphicEnds = FALSE,
    ## Subsample bam files for visualisation with IgvJs in the
    ## DR2S shiny app.
    createIgv = TRUE,
    ## Generate diagnostic plots.
    plot = TRUE
  )
  ##
  ## report() defaults ####
  ##
  report = list(
    ## Maximum number of sequence letters per line in pairwise alignment.
    blockWidth = 80,
    ## Suppress remapping of reads against final consensus.
    remap = TRUE,
    ## Subsample bam files for visualisation with IgvJs in the
    ## DR2S shiny app.
    createIgv = TRUE
  )
```

An example config file:

``` json
{
  "sampleId": "ID123",
  "locus": "A",
  "datadir": "/home//dr2s_data",
  "outdir": "/home/user/dr2s_data/A/ID12912701",
  "reference": "HLA-A*01:01:01:01",
  "longreads": {
    "dir": "pacbio",
    "type": "pacbio",
    "mapper": "minimap"
  },
  "shortreads": {
    "dir": "illumina",
    "type": "illumina",
    "mapper": "bwamem"
  },
  "pipeline": "SR",
  "opts": {
    "mapInit": {
      "includeDeletions": true,
      "includeInsertions": true,
      "callInsertionThreshold": 0.2,
      "microsatellite": true,
      "forceMapping": false,
      "minMapq": 0,
      "topx": false,
      "pickiness": 1,
      "increasePickiness": 1,
      "lowerLimit": 200,
      "updateBackgroundModel": false,
      "createIgv": false,
      "plot": true
    },
    "partitionLongreads": {
      "threshold": 0.2,
      "distAlleles": 2,
      "skipGapFreq": 0.6667,
      "noGapPartitioning": true,
      "selectCorrelatedPositions": false,
      "measureOfAssociation": "cramer.V",
      "proportionOfOverlap": 0.3333,
      "minimumExpectedDifference": 0.06,
      "selectAllelesBy": "distance",
      "minClusterSize": 20,
      "pickiness": 1,
      "lowerLimit": 40,
      "plot": true
    },
    "mapIter": {
      "iterations": 2,
      "columnOccupancy": 0.4,
      "callInsertionThreshold": 0.2,
      "plot": true
    },
    "mapFinal": {
      "includeDeletions": true,
      "includeInsertions": true,
      "callInsertionThreshold": 0.2,
      "trimPolymorphicEnds": false,
      "createIgv": false,
      "plot": true
    },
    "report": {
      "blockWidth": 80,
      "remap": true,
      "createIgv": false
    }
  },
  "format": "json"
}
```

### Longread-only workflow *DR2S-LR*

If you want to find allele sequences only based on longreads you just
need to set `shortreads = NULL` in `createDR2SConf()` and skip the the
`partitionShortreads()` step in the DR2S pipeline.

``` r
x <- InitDR2S(createDR2SConf(
  sample = "ID12912701",
  locus = "DPB1",
  longreads = list(dir = "pacbio", type = "pacbio", mapper = "minimap"),
  datadir = "~/dr2s_data",
  outdir = "~/dr2s_data/output"
)) %>% 
  mapInit() %>%
  partitionLongreads() %>%
  mapIter() %>%
  mapFinal() %>%
  report() %>% 
  cache()
```
