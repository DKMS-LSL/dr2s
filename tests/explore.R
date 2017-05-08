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


x <- DR2S::read_dr2s("~/dr2s_data/output/BOLETH/pacbio/KIR3DL2.pacbio.ref.multialign/")
x

##run_igv(x, 100)

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

x2 <- DR2Smap(
  sample = "HO104",
  locus = "KIR2DL4",
  longreads = list(type = "pacbio", dir = "pacbio"),
  shortreads = list(type = "illumina", dir = "illumina"),
  datadir = "~/dr2s_data",
  outdir = "~/dr2s_data/output",
  reference = "KIR2DL4*00501",
  consensus = "multialign",
  threshold = 0.2
)

x2_r <- x2 %>%
  map0(min_base_quality = 0)# %>%
  #partition_haplotypes() %>%
  #split_reads_by_haplotype() %>%
  #extract_fastq()

#y2 <- y1 %>% map1()

#y3 <- y2 %>%
#  map2() %>%
#  map3() %>%
#  polish() %>%
#  report()


