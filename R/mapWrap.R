## mapping wrapper functions
## 
mapReads <- function(maptag, reffile, readfile, allele, readtype, opts, refname,
                  optsname, force, outdir, minMapq, clean, threshold, maxDepth,
                  minBaseQuality, minNucleotideDepth, includeDeletions, 
                  includeInsertions, mapFun) {
  ## Run mapper
  flog.info("  Mapping ...", name = "info")
  if (missing(refname))
    refname <- ""
  samfile <- mapFun(reffile = reffile, readfile = readfile, allele = allele,
    readtype = readtype, opts = opts, refname  = refname, optsname = optsname,
    force = force, outdir = outdir)
    
  ## Run bam - sort - index pipeline
  flog.info("  Indexing ...", name = "info")
  bamfile <- .bamSortIndex(samfile = samfile, reffile = reffile,
                           minMapq = minMapq, force = force)

  ## Calculate pileup from graphmap produced SAM file
  flog.info("  Piling up ...", name = "info")
  pileup <- Pileup( bamfile, threshold, max_depth = maxDepth,
    min_base_quality = minBaseQuality, min_mapq = minMapq,
    min_nucleotide_depth = minNucleotideDepth,
    include_deletions = includeDeletions,
    include_insertions = includeInsertions)

  if (includeInsertions && is.null(ins(pileup$consmat))) {
    ## TODO check threshold
    pileup <- .pileupIncludeInsertions(x = pileup, threshold = 0.1)
  }
  pileup
}
