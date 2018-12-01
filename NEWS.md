# DR2S 0.0.5

* move to a JSON-based configuration model and expand configuration for the
  individual DR2S pipeline steps.
* add optional read filtering to mapInit().
* add optional restriction of partitioning to correlated (i.e., in-phase) SNPs.
* add a longread-only pipeline.
* cleanup of NAMESPACE

# DR2S 0.0.4
 
* adapt everything to use with the shiny app
* create svg files for all plots to use within the shiny app
* Make windows executables for running igv
* Add refineAlignment to usual reporting workflow
* Refine invocation from config files in yaml format
* move to cowplot for saving figures
* Print a seqlogo
* Add option to choose between read count and distance for alleles from LR partitioning.
* Add option to prevent clustering by gaps. Useful if only few other SNPs are present.
* Add homopolymer count by mode value
* Adapt refineAlignment to cn 1 samples

# DR2S 0.0.3
 
* Add complete workflow for samples with > 2 alleles and 1 allele.
* Add functions to distribute gaps at homopolymer sites and remove a "backgound
gap noise" in Sequel reads.
* Choose clusters by read number.
* Move from `hlatools` and `IPDdata` to the new `ipd.Hsapiens.db` package.
* Restructure output folder.
* Add executable scripts for IGV and most useful postprocessing steps.
* Change everything to relative paths so the dr2s object can be read from 
everywhere.
* Allow fasta files as reference input.
* Refactor.

# DR2S 0.0.2

* Added a `NEWS.md` file to track changes to the package.
* Adapt for HLA-A, -B and -C.
* Fill differing short- and long-read length with gaps. Mostly caused by 
references longer than the sequenced gene.
* More sensitive insert calling.
* Don't allow inserts in the first and last five positions.
* Test with 3 allele samples.
* Add writing and reading MSA in a phylip like format.
* TODO: Test more with KIR.
* TODO: Remove necessety of hlatools.


