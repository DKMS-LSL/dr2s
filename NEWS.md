# DR2S 0.0.3
 
* Add complete workflow for samples with > 2 alleles and 1 allele.
* Add functions to distribute gaps at homopolymer sites and remove a "backgound
gap noise" in Sequel reads.
* Choose clusters by read number .
* Move from `hlatools` and `IPDdata` to the new `ipd.Hsapiens.db` package.
* Restructure output folder.
* Add executable scripts for IGV and most useful postprocessing steps.
* Change everything to relative paths so the dr2s object can be read from 
everywhere
* Allow fasta files as input
* Refactor.


# DR2S 0.0.2

* Added a `NEWS.md` file to track changes to the package.
* Adopt for HLA-A, -B and -C
* Fill differing short- and long-reads length with gaps. Mostly caused by 
references longer than the sequenced gene.
* More sensitive insert calling.
* Don't allow inserts in the first and last five positions.
* Test with 3 allele samples.
* Add writing and reading MSA in a phylip like format.
* TODO: Test more with KIR.
* TODO: Remove necessety of hlatools.


