.bamSortIndex <- function(samfile,
                           reffile,
                           sample = -1,
                           minMapq = 0,
                           threads = 12,
                           threadmem = "4G",
                           force = FALSE,
                           clean = FALSE) {
  samfile <- normalizePath(samfile, mustWork = TRUE)
  reffile <- normalizePath(reffile, mustWork = TRUE)
  ext <- sprintf("%s%ssorted.bam",
                 if (sample > 0 && sample < 1)
                   (sample * 100) %+% "pct."
                 else
                   "",
                 if (minMapq > 0)
                   minMapq %+% "MAPQ."
                 else
                   "")
  sorted <- sub("sam(.gz)?$", ext, samfile)
  ## -F260 exclude 'read unmapped', 'not primary alignment'
  ## Better use -F2308 to also exclude chimeric reads!!!!
  fmt <- paste("samtools view -@%s -F2308 %s -q%s -bT '%s' '%s'",
               "| samtools sort -T /tmp/sorted -m%s -@%s -o '%s' -", 
               "&& samtools index '%s'")
  cmd <- sprintf(fmt, threads, if (sample > 0) "-s" %+% sample else "", minMapq,
            reffile, samfile, threadmem, threads, sorted, sorted)
  ## Don't execute if file exists and force is false
  if (force || !file.exists(sorted)) {
    system(cmd)
  }
  ## Clean up only if the bamfile now exists
  if (clean && file.exists(sorted)) {
    unlink(samfile)
  }

  sorted
}
