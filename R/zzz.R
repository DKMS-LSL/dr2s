## Define global variables for vars in "foreach" and "data.table" to avoid 
## R CMD check warnings
globalVariables(c("lrd", "dst", "ref", "h", "hp", "i", "pos", "sampleId",
                  "nucleotide", "freq", "count", "npoly", "clade", "prob",
                  "haplotype", "read"))

