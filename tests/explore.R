x <- DR2Smap(
  sample = "KGU",
  locus = "KIR3DL2",
  longreads = list(type = "pacbio", dir = "pacbio"),
  shortreads = list(type = "illumina", dir = "illumina"),
  datadir = "~/dr2s_data",
  outdir = "~/dr2s_data/output",
  reference = "KIR3DL2*0010101",
  consensus = "multialign",
  threshold = 0.2
)


##run_igv(x, 100)

y1 <- x %>%
  map0() %>%
  partition_haplotypes() %>%
  split_reads_by_haplotype() %>%
  extract_fastq()

y2 <- y1 %>%
  map1()

y3 <- y2 %>%
  map2()

y4 <- y3 %>%
  map3()

y5 <- y4 %>%
  polish() %>%
  report()

x <- DR2S::read_dr2s("~/dr2s_data/output/JVM/pacbio/KIR3DL2.pacbio.ref.multialign/")
x
run_igv(x, 100)


xx <- DR2Smap(
  sample = "KGU",
  locus = "KIR3DL2",
  longreads = list(type = "pacbio", dir = "pacbio"),
  shortreads = list(type = "illumina", dir = "illumina"),
  datadir = "~/dr2s_data",
  outdir = "~/dr2s_data/output",
  reference = "KIR3DL2*008",
  consensus = "multialign",
  threshold = 0.2
)

yy1 <- xx %>%
  map0() %>%
  partition_haplotypes() %>%
  split_reads_by_haplotype() %>%
  extract_fastq()

yy2  <- yy1 %>% map1()
yy3  <- yy2 %>% map2()
yy4  <- yy3 %>% map3()
yy5  <- yy4 %>%
  polish() %>%
  report()

curr_r <- DR2S::read_dr2s("~/dr2s_data/output/KGU/pacbio/KIR3DL2.pacbio.ref.multialign/")
run_igv(curr_r, 100)
###############################################################################################

se1 <- DR2Smap(
  sample = "JVM",
  locus = "KIR3DL2",
  longreads = list(type = "pacbio", dir = "sequel"),
  shortreads = list(type = "illumina", dir = "illumina"),
  datadir = "~/dr2s_data",
  outdir = "~/dr2s_data/output_sequel",
  reference = "KIR3DL2*0010101",
  consensus = "multialign",
  threshold = 0.2
)

y_se1 <- se1 %>%
  map0(min_base_quality = 0) %>%
  partition_haplotypes() %>%
  split_reads_by_haplotype() %>%
  extract_fastq()

yy2_se1  <- y_se1 %>% map1(min_base_quality = 0)
yy3_se1  <- yy2_se1 %>% map2(min_base_quality = 0)
yy4_se1  <- yy3_se1 %>% map3(min_base_quality = 0)
yy5_se1  <- yy4_se1 %>%
  polish() %>%
  report()



se2 <- DR2Smap(
  sample = "KGU",
  locus = "KIR3DL2",
  longreads = list(type = "pacbio", dir = "sequel"),
  shortreads = list(type = "illumina", dir = "illumina"),
  datadir = "~/dr2s_data",
  outdir = "~/dr2s_data/output_sequel",
  reference = "KIR3DL2*008",
  consensus = "multialign",
  threshold = 0.2
)

y_se2 <- se2 %>%
  map0(min_base_quality = 0) %>%
  partition_haplotypes() %>%
  split_reads_by_haplotype() %>%
  extract_fastq()

yy2_se2  <- y_se2 %>% map1(min_base_quality = 0)
yy3_se2  <- yy2_se2 %>% map2(min_base_quality = 0)
yy4_se2  <- yy3_se2 %>% map3(min_base_quality = 0)
yy5_se2  <- yy4_se2 %>%
  polish() %>%
  report()


se3 <- DR2Smap(
  sample = "DEDKM9562016",
  locus = "KIR3DL2",
  longreads = list(type = "pacbio", dir = "sequel"),
  shortreads = list(type = "illumina", dir = "illumina"),
  datadir = "~/dr2s_data",
  outdir = "~/dr2s_data/output_sequel",
  reference = "KIR3DL2*0010101",
  consensus = "multialign",
  threshold = 0.2
)

y_se3 <- se3 %>%
  map0(min_base_quality = 0) %>%
  partition_haplotypes() %>%
  split_reads_by_haplotype() %>%
  extract_fastq()

yy2_se3  <- y_se3 %>% map1(min_base_quality = 0)
yy3_se3  <- yy2_se3 %>% map2(min_base_quality = 0)
yy4_se3  <- yy3_se3 %>% map3(min_base_quality = 0)
yy5_se3  <- yy4_se3 %>% polish() %>% report()



debug(trim_softclipped_ends)
debug(trim_polymorphic_ends)


samfile <- "~/dr2s_data/output/JVM/pacbio/KIR3DL2.pacbio.ref.multialign/A/"

mapfun <- y2$getMapFun()


maptag   <- "xxx"
refpath  <- y3$map2[["B"]]$seqpath
sreadpath <- y3$getShortreads()
## Run mapper
message("  Mapping short reads against Map2 consensus ...")
samfile <- mapfun(
  reffile  = refpath,
  readfile = sreadpath,
  # if we run shortreads against both pacbio and nanopore data
  # this hack makes sure that we can distinguish the bam files ->
  # we get pacbio.illumina.bwamem.A...bam and nanopore.illumina.bwamem.A...bam
  allele   = "xxx",
  readtype = y2$getSrdType(),
  opts     = list(),
  refname  = "xxx",
  optsname = "xxx",
  force    = FALSE,
  outdir   = "~/tmp/"
)
## Run bam - sort - index pipeline
bamfile <- bam_sort_index(samfile, refpath, clean = TRUE)
## Trim softclips
fq <- trim_softclipped_ends(bam = Rsamtools::scanBam(bamfile)[[1]], preserve_ref_ends = TRUE)
## Trim polymorphic ends
fq <- trim_polymorphic_ends(fq)



fqdir  <- dir_create_if_not_exists(file.path(self$getOutdir(), "merged"))
fqfile <- paste("sread", group, self$getMapper(), "trimmed", "fastq", "gz", sep = ".")
fqout  <- file_delete_if_exists(file.path(fqdir, fqfile))
fqout <- "~/tmp/out.fq"

ShortRead::writeFastq(fq, fqout, compress = TRUE)
file_delete_if_exists(bamfile)


x2 <- DR2Smap(
  sample = "EK",
  locus = "KIR2DL4",
  longreads = list(type = "pacbio", dir = "pacbio"),
  shortreads = list(type = "illumina", dir = "illumina"),
  datadir = "~/dr2s_data",
  outdir = "~/dr2s_data/output",
  reference = "KIR2DL4*0010305",
  consensus = "multialign",
  threshold = 0.2
)

x2_r <- x2 %>%
  map0(min_base_quality = 0)## %>%
  ##partition_haplotypes()## %>%
  ##split_reads_by_haplotype() %>%
  ##extract_fastq()

#y2 <- y1 %>% map1()

#y3 <- y2 %>%
#  map2() %>%
#  map3() %>%
#  polish() %>%
#  report()

x3 <- DR2Smap(
  sample = "JVM",
  locus = "KIR3DL2",
  longreads = list(type = "pacbio", dir = "pacbio"),
  shortreads = list(type = "illumina", dir = "illumina"),
  datadir = "~/dr2s_data",
  outdir = "~/dr2s_data/output",
  reference = "KIR3DL2*0010101",
  consensus = "multialign",
  threshold = 0.2
)

y11 <- x3 %>%
  map0() %>%
  partition_haplotypes() %>%
  split_reads_by_haplotype() %>%
  extract_fastq()

y21 <- y11 %>% map1()

y31 <- y21 %>%
  map2() %>%
  map3() %>%
  polish() %>%
  report()

