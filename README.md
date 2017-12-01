
<!-- README.md is generated from README.Rmd. Please edit that file -->
gDR2S - generic Dual redundant reference sequencing
===================================================

An R package designed to facilitate generating reliable, full-length phase-defined reference sequences for novel HLA and KIR alleles.

**Note, that this package is still maturing. There's no guarantee yet that the API or the under-the-hood workings of the package won't change substantially**

Workflow
--------

### Input and output

`gDR2S` is designed to integrate long-read HLA data (e.g., PacBio or Oxford Nanopore sequences) and short-read shotgun data (Illumina).

As input, we expect long-read and short-read FASTQ files to be placed in separate subdirectories within a working directory and to follow the naming convention `SAMPLEID_LOCUS_.*.fastq(.gz)?`.

`SAMPLEID` can be any arbitrary unique identification code and `LOCUS` should be one of `A`, `B`, `C`, `DQB1`, `DRB1`, or `DPB1`

All output is placed in a directory tree `output\SAMPLEID\READTYPE\LOCUS`.

An example:

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
             +-- ID12912701
                   |
                   +-- pacbio
                         |
                         +-- HLA-DPB1.pacbio.ref.multialign
                     

### Usage

Once all input files are put in place a gDR2S analysis is started with a call to the function `DR2Smap`:

``` r
x <- DR2Smap(
  sample = "ID12912701",
  locus = "DPB1",
  longreads = list(type = "pacbio", dir = "pacbio"),
  shortreads = list(type = "illumina", dir = "illumina"),
  datadir = "~/dr2s_data",
  outdir = "~/dr2s_data/output",
  reference = "04:02:01:01",
  consensus = "mapping",
  iterations = 1,
  partSR = TRUE,
  filterScores = TRUE,
  microsatellite = TRUE,
  threshold = 0.2
)
```

This call generates an `R6` object of class `DR2S` that encapsulates all data and methods for all subsequent analysis steps.

The argument `reference` is the name of the allele against which an initial mapping of the long reads will be performed. A suitable allele can typically be chosen from existing typing information.

The argument `consensus` should be set to `mapping`.

`threshold` is the frequency below which SNP variants are considered noise. num \[0.2\]

`iterations` is the number of iterations to perform in the `mapIter` step. int \[1\]

`microsatellite` decides whether to perform an expansion of the initial reference by an extra mapping step of short reads to the reference. Necessary/useful when you know you have repeats like microsatellites. logical \[TRUE\]

`partSR` decides whether the partitioning is done using SNPs found in shortreads. Usually better to use if you have shortreads. logical \[TRUE\]

`filterScores` parameter for a filtering step in mapInit. Filters scores below a threshold of 0.3\*sequence length. Necessary for KIR, but slightly worse results for e.g. DPB1 genes. logical \[TRUE\]

An analysis proceeds in a number of steps that can be chained together using the pipe `%>%`:

``` r
x %>% 
  mapInit()
  partition_haplotypes() %>%
  split_reads_by_haplotype() %>%
  extract_fastq()
  mapIter()
  partitionShortReads %>%
  mapFinal()
  polish() %>%
  report
```

The individual steps perform the following analyses:

-   **mapInit:** Map the long and short reads against a reference sequence constructed from `reference` using `bwa`.

1.  Map the shortreads against the provided reference. Construct the consensus.
2.  A second mapping of shortreads to the consensus of first Mapping if microsatellites is true. This expands the reference and is necessary of you know that there are repeats like in microsatellites. Construct the consensus.
3.  A third mapping of shortreads to infer polymorphic positions and and read IDs at the polymorphic positions.
4.  Map the longreads to the consensus latest consensus.

-   **partition\_haplotypes:**

1.  Construct a SNP matrix from the longreads at polymorphic positions in shortreads.
2.  Perform hierarchical clustering,
3.  chimera detection and
4.  assign a haplotype and score to each read.

-   **split\_reads\_by\_haplotype:** Pick longreads that best represent the allele haplotypes based on a *haplotype score* from the previous step.
-   **extract\_fastq:** Extract the chosen longreads from the alignment file and write the data as FASTQs into subdirectories for each haplotype.
-   **mapIter:**

1.  Construct a consensus sequence for each haplotype from the initial mapping while using only reads that are assigned to a haplotype.
2.  Iteratively refine the long read consensus sequences by remapping long reads per haplotype to the latest haplotype consensus.

-   **partitionShortReads** Partition the shortreads into haplotypes using the inferred longread haplotypes and initial mapping.
-   **mapFinal:** Map the short reads against the refined long read consensus sequences.
-   **polish:** Disambiguate polymorphic positions in the short read mapping. Check and report inconsistencies, insertions, deletions along the way. This should not be that much with okayish seuence data
-   **report:** Report the finalised short-read-based consensus sequences as FASTA files. Provide a tsv file with suspicious positions that may warrant manual inspection. Report the alignment of all haplotypes.

Throughout this process a number of diagnostic plots are produced and placed in the output directory for later inspection.

While this process works remarkably well, there are situations where alignment artefacts or plain bad luck may introduce errors in the final consensus sequences. You should **never** accept the result as ground truth without some manual and visual consistency checks!

`gDR2S` provides some facilities to aid checking and signing-off of finalised consensus sequences.

A typical post-processing workflow may look as follows:

``` r
plot_diagnostic_alignment(x)
run_igv(x, 3000)
check_alignment_file(x)
report_checked_consensus(x)
```

-   **plot\_diagnostic\_alignment:** Displays an alignment of three preliminary long-read-based consensus sequences and the final short-read-based consensus sequence for all alleles in your browser. This can be used to spot inconsistencies between the long-read and the short-read evidence.

-   **run\_igv:** Opens an instance of the IGV Genome Browser for each haplotype at a specified position (one for each allele) displaying both the long read and short read data for manual inspection.

-   **check\_alignment\_file :** Opens a pairwise alignment of the final consensus sequences in your text editor. Use this to perform any manual edits on the consensus sequences. Editor options are: "subl", "gvim" and "gedit". Defaults to systems standard editor.

-   **report\_checked\_consensus:** Export final consensus sequences from the edited pairwise or multiple alignment as FASTAs into a separate subdirectory `./checked` in the output directory.

### Only longreads

If you want to find allele sequences only based on longreads you need to skip the partitioning of shortreads. The longread partitioning in other steps is automatically done if `partSR` is set to `FALSE` and no shortreads are found in the datadir.

``` r
x %>% 
  mapInit()
  partition_haplotypes() %>%
  split_reads_by_haplotype() %>%
  extract_fastq()
  mapIter()
  mapFinal()
  polish() %>%
  report
```

Installation
------------

This package is only available via gitlab now. It depends on a local installation of `samtools`, `bwa` (&gt;= 0.7.11), python and the pysam library, and a C++11 compliant compiler. An alternative mapper with better results is `minimap2`, a successor of `bwa` For installation over git via SSH, you need to have installed libssh2-1 and libssh2-1-dev prior to install git2r/devtools or reinstall the package afterwards. You do also need to have set up your ssh key in gitlab. Otherwise you need to give your gitlab password in plain.

``` r
install.packages("devtools")  # if not already installed
library("devtools")
install_git("git@srvddgit01.labor.local:rlib/IPDdata.git")
install_git("git@srvddgit01.labor.local:rlib/hlatools.git")
install_git("git@srvddgit01.labor.local:sklas/gDR2S.git")
```
