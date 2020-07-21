# DR2S 1.1.16

* [bugfix] fix isues with modified paired-end fastqs by renaming reads after the first mapping and trimming

# DR2S 1.1.15

* [bugfix] fix issue with bpparallel in the remapping step that only occurs 
  on the slurm cluster.

# DR2S 1.1.14

* [bugfix] fixes multiple issues with parallel access to bam files that only
  arise on the SLURM cluster.
* [bugfix] depend on dtplyr 1.0.1 for correct handling of .data pronouns.

# DR2S 1.1.13

* [bugfix] allow for homozygous samples when downsampling for IGV.

# DR2S 0.0.6

* improved identification of highly associated polymorphic positions.
* addition of diagnistic plots `plot.mclust.png` and `plot-mclust.ovl.png`.
* minor bugfixes.

# DR2S 0.0.5

* implement a JSON-based configuration model and expand configuration options
  for the individual DR2S pipeline steps.
* add optional read filtering to `mapInit()` and add "readpicking" plots.
* add optional restriction of partitioning to correlated (i.e., in-phase) SNPs
  for the lonread-only workflow.
* add correlogram and association plots.
* add a dedicated longread-only pipeline.
* include standard reference data for HLA.
* move from pdf to png for plots.
* cleanup of dependencies.

# DR2S 0.0.4
 
* adapt package for use with a DR2S shiny app.
* create svg plots for use within the shiny app.
* make windows executables for running igv.
* add `refineAlignment()` to usual reporting workflow.
* refine invocation from config files in yaml format.
* move to `cowplot` for saving figures.
* Add seqlogo plot.
* Add option to choose between read count and distance for alleles from LR partitioning.
* Add option to prevent clustering by gaps. Useful if only few other SNPs are present.
* Add homopolymer count by mode value.
* Adapt `refineAlignment()` to cn 1 samples.

# DR2S 0.0.3
 
* Add complete workflow for samples with >2 alleles and 1 allele.
* Add functions to rightshift gaps at homopolymer sites and remove a "backgound
  gap noise" in Sequel and Nanopore reads.
* Choose clusters by read number.
* Move from `hlatools` and `IPDdata` to the new `ipd.Hsapiens.db` package for
  creating references.
* Restructure output folder.
* Add executable scripts for IGV and the most useful postprocessing steps.
* Use relative paths so the dr2s object can be read from everywhere.
* Allow fasta files as custom reference.
* Internal code refactoring.

# DR2S 0.0.2

* Added a `NEWS.md` file to track changes to the package.
* Adapt for HLA-A, -B and -C.
* Fill differing short- and long-read length with gaps. Mostly caused by
  references longer than the sequenced gene.
* More sensitive insert calling.
* Don't allow inserts in the first and last five positions.
* Tested with 3 allele samples.
* Add writing and reading MSA in a phylip-like format.


