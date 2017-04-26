x <- DR2Smap(
  sample = "BOLETH",
  locus = "KIR3DL2",
  longreads = list(type = "pacbio", dir = "pacbio"),
  shortreads = list(type = "illumina", dir = "illumina"),
  datadir = "~/dr2s_data",
  outdir = "~/dr2s_data/output",
  reference = "KIR3DL2*00501",
  consensus = "multialign",
  threshold = 0.2
)

y1 <- x %>%
  map0() %>%
  partition_haplotypes() %>%
  split_reads_by_haplotype() %>%
  extract_fastq()

y2 <- y1 %>% map1()

y3 <- y2 %>%
  map2() %>%
  map3() %>%
  polish() %>%
  report()

