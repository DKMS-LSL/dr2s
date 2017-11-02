# Method: MapInit ####

#' @export
mapInit.DR2S <- function(x,
                      opts = list(),
                      optsname = "",
                      partSR = TRUE,
                      pct = 100,
                      threshold = 0.20,
                      min_base_quality = 3,
                      min_mapq = 0,
                      max_depth = 1e4,
                      min_nucleotide_depth = 3,
                      include_deletions = TRUE,
                      include_insertions = TRUE,
                      include_read_ids = TRUE,
                      microsatellite = FALSE,
                      force = FALSE,
                      fullname = TRUE,
                      plot = TRUE) {
  flog.info("Step 0: Initial mapping of long reads ... ", name = "info")
  Sys.sleep(1)
  x$runMapInit(opts = opts, optsname = optsname, partSR = partSR, pct = pct, threshold = threshold,
            min_base_quality = min_base_quality, min_mapq = min_mapq,
            max_depth = max_depth, min_nucleotide_depth = min_nucleotide_depth,
            include_deletions = include_deletions, include_insertions = include_insertions,
            include_read_ids = include_read_ids, microsatellite = microsatellite,
            force = force, fullname = fullname, plot = plot)
  message("  Done!\n")
  invisible(x)
}
DR2S_$set("public", "runMapInit",
         function(opts = list(),
                  optsname = "",
                  partSR = TRUE,
                  pct = 100,
                  threshold = 0.20,
                  min_base_quality = 3,
                  min_mapq = 0,
                  max_depth = 1e4,
                  min_nucleotide_depth = 3,
                  include_deletions = TRUE,
                  include_insertions = TRUE,
                  include_read_ids = TRUE,
                  force = FALSE,
                  fullname = TRUE,
                  microsatellite = FALSE,
                  plot = TRUE) {

           # debug
           # opts = list()
           # partSR = TRUE
           # optsname = ""
           # pct = 100
           # threshold = 0.20
           # min_base_quality = 3
           # min_mapq = 0
           # max_depth = 1e4
           # min_nucleotide_depth = 3
           # include_insertions = TRUE
           # include_deletions = TRUE
           # include_read_ids = TRUE
           # force = FALSE
           # fullname = TRUE
           # plot = TRUE
           # microsatellite = TRUE
           # self <- dpb1_3


           stopifnot(pct > 0 && pct <= 100)
           recode_header <- FALSE
           ## Overide default arguments
           args <- self$getOpts("mapInit")
           if (!is.null(args)) {
             env  <- environment()
             list2env(args, envir = env)
           }

           microsatellite <- self$getMicrosatellite()
           partSR <- self$getPartSR()

           if (recode_header) {
             stopifnot(
               recode_fastq_header(self$getLongreads()),
               recode_fastq_header(self$getShortreads()[1]),
               recode_fastq_header(self$getShortreads()[2])
             )
           }
           if (partSR){
             mapfmt  <- "mapInit1 <%s> <%s> <%s> <%s>"
             maptag  <- sprintf(mapfmt, self$getReference(), self$getSrdType(), self$getMapper(), optstring(opts, optsname))
             if (threshold != self$getThreshold()) {
               self$setThreshold(threshold)
             }
             ## Fetch mapper
             map_fun <- self$getMapFun()
             ## Run mapper
             flog.info(" First mapping of shortreads to provided reference", name = "info")
             samfile <- map_fun(
               reffile  = self$getRefPath(),
               readfile = self$getShortreads(),
               allele   = self$getReference(),
               readtype = self$getSrdType(),
               opts     = opts,
               refname  = "",
               optsname = optsname,
               force    = force,
               outdir   = self$getOutdir()
             )

             ## Run bam - sort - index pipeline
             flog.info("  Indexing ...", name = "info")
             bamfile <- bam_sort_index(
               samfile = samfile,
               reffile = self$getRefPath(),
               sample = pct / 100,
               min_mapq = min_mapq,
               force = force,
               # clean = TRUE
               # debug
               clean = FALSE
             )

             ## Calculate pileup from graphmap produced SAM file
             flog.info("  Piling up ...", name = "info")
             pileup <- Pileup(
               bamfile,
               self$getThreshold(),
               max_depth,
               min_base_quality = min_base_quality,
               min_mapq = min_mapq,
               min_nucleotide_depth = min_nucleotide_depth,
               include_deletions = include_deletions,
               include_insertions = include_insertions
             )
             if (include_insertions && is.null(ins(pileup$consmat))) {
               pileup <- pileup_include_insertions(pileup)
             }

             # calc initial consensus
             flog.info("  Construct initial consensus from shortreads", name = "info")
             conseq <- conseq(pileup$consmat, "mapInit", "prob", exclude_gaps = TRUE, threshold = self$getThreshold())
             conseq_name <- paste0("Init.consensus.", sub(".sam.gz", "", basename(samfile)))
             conseqpath     <- file.path(self$getOutdir(), paste0(conseq_name, ".fa"))
             Biostrings::writeXStringSet(
               Biostrings::DNAStringSet(gsub("[-+]", "N", conseq)),
               conseqpath)

             if (microsatellite){
               mapfmt  <- "mapInit1.2 <%s> <%s> <%s> <%s>"
               maptag  <- sprintf(mapfmt, conseq_name, self$getSrdType(), self$getMapper(), optstring(opts, optsname))
               if (threshold != self$getThreshold()) {
                 self$setThreshold(threshold)
               }
               ## Fetch mapper
               map_fun <- self$getMapFun()
               ## Run mapper
               flog.info(" Refine microsatellites or repeats by extending reference ", name = "info")
               flog.info("  Second mapping of shortreads to reference from first mapping ", name = "info")
               samfile <- map_fun(
                 reffile  = conseqpath,
                 readfile = self$getShortreads(),
                 allele   = conseq_name,
                 readtype = self$getSrdType(),
                 opts     = opts,
                 refname  = "",
                 optsname = optsname,
                 force    = force,
                 outdir   = self$getOutdir()
               )

               ## Run bam - sort - index pipeline
               flog.info("   Indexing ...", name = "info")
               bamfile <- bam_sort_index(
                 samfile = samfile,
                 reffile = conseqpath,
                 sample = pct / 100,
                 min_mapq = min_mapq,
                 force = force,
                 # clean = TRUE
                 # debug
                 clean = FALSE
               )

               ## Calculate pileup from graphmap produced SAM file
               flog.info("   Piling up ...", name = "info")
               pileup <- Pileup(
                 bamfile,
                 self$getThreshold(),
                 max_depth,
                 min_base_quality = min_base_quality,
                 min_mapq = min_mapq,
                 min_nucleotide_depth = min_nucleotide_depth,
                 include_deletions = include_deletions,
                 include_insertions = include_insertions
               )

               if (include_insertions && is.null(ins(pileup$consmat))) {
                 pileup <- pileup_include_insertions(pileup)
               }

               # Infer initial consensus
               flog.info("   Construct second consensus from shortreads with refined repeats ...", name = "info")
               conseq <- conseq(pileup$consmat, "mapInit1.2", "prob", force_exclude_gaps = TRUE, threshold = 0.2)
               conseq_name <- paste0("Init.consensus.2", sub(".sam.gz", "", basename(samfile)))
               conseqpath     <- file.path(self$getOutdir(), paste0(conseq_name, ".fa"))
               Biostrings::writeXStringSet(
                 Biostrings::DNAStringSet(gsub("[-+]", "N", conseq)),
                 conseqpath)
             }
             # Debug
             #browse_seqs(DECIPHER::AlignSeqs(Biostrings::DNAStringSet(c(cs, conseq))))

             mapInitSR1 = structure(
               list(
                 reads   = self$getShortreads(),
                 bamfile = bamfile,
                 pileup  = pileup,
                 tag     = maptag,
                 conseq  = conseq,
                 seqpath = conseqpath,
                 ref     = conseq_name
               ),
               class  = c("mapInit", "list")
             )

             ## Second mapping to infer polymorphic positions from same reference as longreads
             mapfmt  <- "mapInit2 <%s> <%s> <%s> <%s>"
             maptag  <- sprintf(mapfmt, mapInitSR1$ref,self$getSrdType(), self$getMapper(), optstring(opts, optsname))

             ## Fetch mapper
             map_fun <- self$getMapFun()
             ## Run mapper
             flog.info(" Mapping shortreads against mapInit reference for calling SNPs ...", name = "info")
             samfile <- map_fun(
               reffile  = mapInitSR1$seqpath,
               readfile = self$getShortreads(),
               allele   = mapInitSR1$ref,
               readtype = self$getSrdType(),
               opts     = opts,
               refname  = "",
               optsname = optsname,
               force    = force,
               outdir   = self$getOutdir()
             )

             ## Run bam - sort - index pipeline
             flog.info("  Indexing ...", name = "info")
             bamfile <- bam_sort_index(
               samfile = samfile,
               reffile = mapInitSR1$seqpath,
               sample = pct / 100,
               min_mapq = min_mapq,
               force = force,
               # clean = TRUE
               # debug
               clean = FALSE
             )

             ## Calculate pileup from graphmap produced SAM file
             flog.info("  Piling up ...", name = "info")
             pileup <- Pileup(
               bamfile,
               self$getThreshold(),
               max_depth,
               min_base_quality = min_base_quality,
               min_mapq = min_mapq,
               min_nucleotide_depth = min_nucleotide_depth,
               include_deletions = include_deletions,
               include_insertions = include_insertions
             )

             mapInitSR2 = structure(
               list(
                 reads   = self$getShortreads(),
                 bamfile = bamfile,
                 pileup  = pileup,
                 tag     = maptag,
                 conseq  = conseq,
                 seqpath = conseqpath,
                 ref     = conseq_name
               ),
               class  = c("mapInit", "list")
             )
           }



           flog.info(" Map longreads against MapInit reference for clustering ...", name = "info")
           ref <- ifelse(partSR, mapInitSR1$seqpath, self$getRefPath())
           allele <- ifelse(partSR, mapInitSR1$ref, self$getReference())

           mapfmt  <- "mapInit <%s> <%s> <%s> <%s>"
           maptag  <- sprintf(mapfmt, allele, self$getLrdType(), self$getMapper(), optstring(opts, optsname))
           if (threshold != self$getThreshold()) {
             self$setThreshold(threshold)
           }
           ## Fetch mapper
           map_fun <- self$getMapFun()
           ## Run mapper
           flog.info("  Mapping ...", name = "info")
           samfile <- map_fun(
             reffile  = ref,
             readfile = self$getLongreads(),
             allele   = allele,
             readtype = self$getLrdType(),
             opts     = opts,
             refname  = "",
             optsname = optsname,
             force    = force,
             outdir   = self$getOutdir()
           )

           ## Run bam - sort - index pipeline
           flog.info("  Indexing ...", name = "info")
           bamfile <- bam_sort_index(
             samfile = samfile,
             reffile = ref,
             sample = pct / 100,
             min_mapq = min_mapq,
             force = force,
             # clean = TRUE
             # debug
             clean = FALSE
           )

           ## Calculate pileup from graphmap produced SAM file
           flog.info("  Piling up ...", name = "info")
           pileup <- Pileup(
             bamfile,
             self$getThreshold(),
             max_depth,
             min_base_quality = min_base_quality,
             min_mapq = min_mapq,
             min_nucleotide_depth = min_nucleotide_depth,
             include_deletions = TRUE,
             include_insertions = FALSE
           )

           self$mapInit = structure(
             list(
               reads   = self$getLongreads(),
               bamfile = bamfile,
               pileup  = pileup,
               tag     = maptag,
               SR1     = NULL,
               SR2     = NULL
             ),
             class  = c("mapInit", "list")
           )

           if (partSR){
             self$mapInit$SR1 <- mapInitSR1
             self$mapInit$SR2 <- mapInitSR2
           }
           if (plot) {
             flog.info(" Plotting MapInit Summary ...", name = "info")
             ## Coverage and frequency of minor alleles
             gfile <- file.path(self$getOutdir(), paste0("plot.MapInit.LR.", sub("bam$", "pdf", usc(basename(self$mapInit$bamfile)))))
             pdf(file = gfile, width = 12, height = 10, onefile = TRUE, title = paste(self$getLocus(), self$getSampleId(), sep = "." ))
             self$plotmapInitSummary(
               thin = 0.25,
               width = 2,
               label = self$getMapTag("init", "LR"),
               readtype = "LR"
             )
             dev.off()

             if (partSR){
               gfile <- file.path(self$getOutdir(), paste0("plot.MapInit.SR.", sub("bam$", "pdf", usc(basename(self$mapInit$SR2$bamfile)))))
               pdf(file = gfile, width = 12, height = 10, onefile = TRUE, title = paste(self$getLocus(), self$getSampleId(), sep = "." ))
               self$plotmapInitSummary(
                 thin = 0.25,
                 width = 2,
                 label = self$getMapTag("init", "SR"),
                 readtype = "SR"
               )
               dev.off()

             }
           }

           return(invisible(self))
         })

#' @export
print.mapInit <- function(x, ...) {
  msg <- sprintf("An object of class '%s'\n", class(x)[1])
  msg <- sprintf(
    "%s [Tag]     %s\n [Reads]   %s\n [Bamfile] %s\n [Pileup]\n",
    msg,
    x$tag,
    basename(x$reads),
    basename(x$bamfile)
  )
  cat(msg)
}

## Method: partition_haplotypes ####

#' @export
partition_haplotypes.DR2S <- function(x,
                                      max_depth = 1e4,
                                      skip_gap_freq = 2/3,
                                      plot = TRUE,
                                      ...) {
  flog.info("Step 1: Partition reads into haplotypes ...", name = "info")
  Sys.sleep(1)
  x$runHaplotypePartitioning(max_depth = max_depth,
                             threshold = NULL,
                             skip_gap_freq = skip_gap_freq,
                             plot = plot)
  message("  Done!\n")
  invisible(x)

}
#self <- dpb1_3
DR2S_$set("public", "runHaplotypePartitioning", function(max_depth = 1e4,
                  skip_gap_freq = 2/3,
                  threshold = NULL,
                  plot = TRUE) {

                  # debug
                  # skip_gap_freq = 2/3
                  # threshold = NULL
                  # plot = TRUE
           stopifnot(self$hasPileup())

           ## Overide default arguments
           args <- self$getOpts("partition")
           if (!is.null(args)) {
             env  <- environment()
             list2env(args, envir = env)
           }
           if (is.null(threshold)){
             threshold <- self$getThreshold()
           }

           tag  <- self$getMapTag("init", "LR")
           th   <- self$getThreshold()

           if (!is.null(self$mapInit$SR1)){
             useSR <- TRUE
             refseq <- self$mapInit$SR1$conseq
             flog.info(" Generating SNP matrix from shortreads ...", name = "info")
           } else {
             useSR <- FALSE
             refseq <- self$getRefSeq()
             flog.info(" Generating SNP matrix from longreads ...", name = "info")
           }
           ppos <- self$polymorphicPositions(useSR = useSR)

           if (NROW(ppos) == 0){

             stop("No polymorphic position!")
           }

           mat <- if (tryCatch(
             !is(self$partition, "PartList"),
             error = function(e)
               TRUE
           ) ||
           !(all(ppos$position %in% colnames(self$partition$mat)) &&
             all(colnames(self$partition$mat) %in% ppos$position))) {
             SNPmatrix(bamfile = self$mapInit$bamfile, refseq = refseq, polymorphic_positions = ppos)
           } else {
             self$partition$mat
           }
           ## expdev: try seqlogo
           # mat <- hla$partition$mat
           # seqs <- lapply(rownames(mat), function(x) paste(mat[x, ], collapse = ""))#paste0(self$partition$mat[x,])))
           # names(seqs) <- rownames(mat)
           #
           # devtools::install_github("omarwagih/ggseqlogo")
           # library(ggseqlogo)
           # library(ggplot2)
           # p <- ggplot() + geom_logo(unlist(seqs), seq_type = "dna") + theme_logo()
           # DECIPHER::AlignProfiles()
           # ggsave(filename = "tset.pdf", plot= p, width = 15)
           # ggseqlogo(unlist(seqs), seqtype = "DNA")

           flog.info(" Partitioning %s reads over %s SNPs ...", NROW(mat), NCOL(mat), name = "info")
           prt <- partition_reads(x = mat, skip_gap_freq = skip_gap_freq, deepSplit = 1, threshold = threshold)
           browse_seqs(SQS(prt), file = file.path(self$getOutdir(), "partition.fa.html"), openURL = FALSE)

           self$setHapTypes(levels(as.factor(PRT(prt))))

           self$partition = structure(list(
             mat = mat,
             prt = prt,
             hpl = NULL,
             lmt = NULL
           ),
           class = c("PartList", "list"))
           # if (plot) {
           if (FALSE) {
             message("  Plotting ...")
             bnm   <- basename(self$mapInit$bamfile)
             outf  <- file.path(self$getOutdir(), paste0("plot0.partition.", sub("bam$", "pdf", usc(bnm))))
             pdf(file = outf, width = 10, height = 12, onefile = TRUE, title = paste(self$getLocus(), self$getSampleId(), sep = "." ))
             self$plotPartitionSummary(label = tag)
             dev.off()
           }



           return(invisible(self))
         })

#' @export
print.PartList <- function(x, ...) {
  msg <- sprintf("An object of class '%s'.\n", class(x)[1])
  msg <- sprintf("%s [Matrix]    %s reads; %s polymorphic positions\n",
                 msg, NROW(x$mat), NCOL(x$mat))
  cat(msg)
  cat(" [Partition] ")
  print(x$prt, ncols = 4, nrows = 6)
  cat(" [Haplotypes] ")
  print(x$hpl)
}

## Method: split_reads_by_haplotype ####

#' @export
split_reads_by_haplotype.DR2S <- function(x,
                                          limits,
                                          ...) {
  flog.info(" Split partitioned reads by score ...", name = "info")
  x$splitReadsByHaplotype( limits)
  invisible(x)
}

#self <- dpb1_3
DR2S_$set("public", "splitReadsByHaplotype", function(limits,
                                                      plot = TRUE){
           stopifnot(self$hasPartition())

           ## Overide default arguments
           args <- self$getOpts("split")
           if (!is.null(args)) {
             env  <- environment()
             list2env(args, envir = env)
           }

           prt  <- partition(self$getPartition())
           haplotypes <- levels(prt$haplotype)

           # Set all limits to NULL
           self.setLimits <- sapply(haplotypes, function(x) NULL)

           prts <-  lapply(haplotypes, function(x) prt[prt$haplotype == x,])
           names(prts) <- haplotypes
           scores <- lapply(prts, function(x) x$mcoef)

           lmts <- optimal_partition_limits(scores)
           plainlmts <- as.list(lmts$limits$c)
           names(plainlmts) <- lmts$limits$haplotype

           self$setLimits(plainlmts)
           self$partition$lmt <- lmts$plt


           # Get only reads within the limit
           reads <- lapply(names(self$getLimits()), function(x) dplyr::filter(prt, haplotype == x, mcoef >= self$getLimits()[x]))
           names(reads) <- names(self$getLimits())
           for (hp in haplotypes) flog.info("  %s: Using %s longreads with score > %s ...", hp, nrow(reads[[hp]]), plainlmts[hp], name = "info")

           # results Structure
           resStruct <- lapply(names(reads), function(x) {
             tr <- reads[[x]]
             structure(
               tr$read,
               q = tr$mcoef,
               freq = NROW(tr)/NROW(dplyr::bind_rows(reads)),
               limit = self$getLimits()[[x]]
             )})
           names(resStruct) <- haplotypes
           self$partition$hpl = structure(resStruct, class = c("HapList", "list"))
           if (plot) {
             tag   <- self$getMapTag("init", "LR")
             bnm   <- basename(self$mapInit$bamfile)
             outf  <- file.path(self$getOutdir(), paste0("plot.partition.", sub("bam$", "pdf", usc(bnm))))
             pdf(file = outf, width = 10, height = 12, onefile = TRUE, title = paste(self$getLocus(), self$getSampleId(), sep = "." ))
             self$plotPartitionSummary(label = tag, limits = unlist(self$getLimits()))
             # debug
             # label = tag
             # limits = unlist(self$getLimits())
             dev.off()
           }

           return(invisible(self))
         })

#' @export
print.HapList <- function(x, ...) {
  msg <- sprintf("An object of class '%s'.\n", class(x)[1])
  for (haplotype in names(x)) {
    msg <- paste0(msg,
             sprintf("%s: n %s; frequency %s; limit %s\n",
               haplotype,
               length(x[haplotype]),
               round(attr(x[[haplotype]], "freq"), 3),
               attr(x[[haplotype]], "limit")))
  }
  cat(msg)
}

#' @export

## ToDo: change this
summary.HapList <- function(object, ....) {
  data.frame(
    n = c(length(object$A), length(object$B)),
    freq = round(c(
      attr(object$A, "freq"), attr(object$B, "freq")
    ), 2),
    limits = c(attr(object$A, "limit"), attr(object$B, "limit")),
    row.names = c("A", "B")
  )
}

## Method: extract_fastq ####

#' @export
extract_fastq.DR2S <- function(x, nreads = NULL, replace = FALSE, nalign = 40, ...) {
  flog.info("Step2: Extracting haplotyped reads ...", name = "info")
  Sys.sleep(1)
  x$extractFastq(nreads = nreads, replace = replace, nalign = nalign)
  message("  Done!\n")
  invisible(x)
}
# debug
#self <- dpb1_3
DR2S_$set("public", "extractFastq",
          function(nreads = NULL, replace = FALSE, nalign = 40) {

            # debug
            # nreads = NULL
            # replace = FALSE
            # nalign = 40

            stopifnot(self$hasHapList())

            ## Overide default arguments
            args <- self$getOpts("extract")
            if (!is.null(args)) {
              env  <- environment()
              list2env(args, envir = env)
            }

            ## override nreads in config if not NULL
            if (!is.null(nreads))
              self$setNreads(nreads)

            ## do this for each haptype
            hptypes <- self$getHapTypes()
            for (hptype in hptypes){
              dir <- dir_create_if_not_exists(file.path(self$getOutdir(), hptype))

              qnames <- self$getHapList(hptype)

                # x = self$mapInit$bamfile
                # qnames = qnames
                # n = self$getNreads()
                # replace = replace

              fq  <- extract_fastq(
                x = self$mapInit$bamfile,
                qnames = qnames,
                n = self$getNreads(),
                replace = replace
              )
              if (!is.null(nreads)) {
                attr(self$partition$hpl[hptype], "index") <- which(qnames %in% as.character(ShortRead::id(fq)))
                attr(self$partition$hpl[hptype], "n") <- self$getNreads()
              }
              file <- paste("hap", hptype, self$getLrdType(), self$getMapper(),
                             paste0("lim", 100 * abs(attr(self$getHapList(hptype), "limit"))),
                             paste0("n", length(fq)),
                             "fastq", "gz", sep = ".")
              out <- file_delete_if_exists(file.path(dir, file))
              ShortRead::writeFastq(fq, out, compress = TRUE)
              self$mapIter[["0"]][[hptype]] = structure(
                list(
                  dir     = dir,
                  reads   = out,
                  ref     = NULL,
                  bamfile = NULL,
                  pileup  = NULL,
                  conseq  = NULL,
                  seqpath = NULL,
                  params  = NULL,
                  tag     = NULL
                ),
                class = c("mapIter", "list")
              )

              # if (self$hasMultialign()) {
              # if (FALSE) {
              #   message("   Constructing initial ", hptype, " consensus from ", nalign, " reads ...")
              #   # msa <- multialign(self, hptype, nalign)
              #   # DECIPHER::BrowseSeqs(msa)
              #   cons <- multialign_consensus(multialign(self, hptype, nalign))
              #   consfile <- paste("consensus", self$getLrdType(), "multialign", paste0(hptype, nalign), "fa", sep = ".")
              #   names(cons) <- consfile
              #   consout <- file_delete_if_exists(file.path(dir, consfile))
              #   Biostrings::writeXStringSet(cons, consout)
              #   self$mapIter[["0"]][[hptype]]$ref <- "multialign"
              #   self$mapIter[["0"]][[hptype]]$conseq <- cons
              #   self$mapIter[["0"]][[hptype]]$seqpath <- consout
              #   self$mapIter[["0"]][[hptype]]$tag <- "multialign"
              # } else {
              #}
            }

            return(invisible(self))
          })

## Method: mapIter ####

#' @export
mapIter.DR2S <- function(x,
                      opts = list(),
                      iterations = 1,
                      pct = 100,
                      min_base_quality = 3,
                      min_mapq = 0,
                      max_depth = 1e4,
                      min_nucleotide_depth = 3,
                      include_deletions = TRUE,
                      include_insertions = TRUE,
                      pruning_cutoff = 0.75,
                      force = FALSE,
                      fullname = TRUE,
                      plot = TRUE) {
  flog.info("Step 3: Iteratively mapping partitioned long reads against consensus sequences from initial mapping using only partitioned reads... ", name = "info")
  Sys.sleep(1)
  x$runMapIter(opts = opts, iterations = 1, pct = pct, min_base_quality = min_base_quality,
            min_mapq = min_mapq, max_depth = max_depth,
            min_nucleotide_depth = min_nucleotide_depth,
            include_deletions = include_deletions,
            include_insertions = include_insertions,
            pruning_cutoff = pruning_cutoff, force = force,
            fullname = fullname, plot = plot)
  message("  Done!\n")
  invisible(x)
}

# self <- hla.fq
#self <- dpb1_3
DR2S_$set("public", "runMapIter",
         function(opts = list(),
                  iterations = 1,
                  pct = 100,
                  min_base_quality = 3,
                  min_mapq = 0,
                  max_depth = 1e4,
                  min_nucleotide_depth = 3,
                  include_deletions = TRUE,
                  include_insertions = TRUE,
                  pruning_cutoff = 0.75,
                  force = FALSE,
                  fullname = TRUE,
                  plot = TRUE) {
         # # debug

         # self <- hla.fq

         # opts = list()
         # pct = 100
         # min_base_quality = 3
         # min_mapq = 0
         # max_depth = 1e4
         # min_nucleotide_depth = 3
         # include_deletions = TRUE
         # include_insertions = TRUE
         # pruning_cutoff = 0.75
         # force = FALSE
         # fullname = TRUE
         # plot = TRUE
         # iterations = 1
         # ##

           ## Overide default arguments
           args <- self$getOpts("mapIter")
           if (!is.null(args)) {
             env  <- environment()
             list2env(args, envir = env)
           }

           ## Mapper
           map_fun <- self$getMapFun()

           hptypes <- self$getHapTypes()
           iterations <- self$getIterations()

           # Construct consensus from initial mapping with the clustered reads
           flog.info(" Constructing a consensus from initial mapping ...", name = "info")
           bamfile <- self$mapInit$bamfile
           if (self$getPartSR()){
             ref <- self$mapInit$SR1$conseq
           } else {
             ref <- self$getRefSeq()
           }
           mat <- msa_from_bam(bamfile, ref, paddingLetter = "-")

           foreach(hptype = hptypes) %do% {
             flog.info("  Constructing a consensus from initial mapping for haplotype %s ...", hptype, name = "info")
             readIds <- self$getHapList(hptype)
             cmat <- Biostrings::consensusMatrix(mat[readIds], as.prob = TRUE)[VALID_DNA(include = "indel"),]
             conseq_name <- paste0("consensus.mapIter.0.", hptype)
             conseq <- conseq(x = t(cmat), name = conseq_name, type = "prob",
                    force_exclude_gaps = TRUE, prune_matrix = TRUE,
                    cutoff = pruning_cutoff)

             seqpath     <- file.path(self$mapIter[["0"]][[hptype]]$dir, paste0(conseq_name, ".fa"))
             self$mapIter[["0"]][[hptype]]$ref <- "mapIter0"
             self$mapIter[["0"]][[hptype]]$conseq <- conseq
             self$mapIter[["0"]][[hptype]]$seqpath <- seqpath
             self$mapIter[["0"]][[hptype]]$tag <- "mapIter0"
             Biostrings::writeXStringSet(
               Biostrings::DNAStringSet(gsub("[-+]", "", conseq)),
               seqpath
             )
           }

           #iteration = 1
           for (iteration in 1:iterations) {
             flog.info(" Iteration %s ", iteration, " of ", iterations, name = "info")

             iterationC = as.character(iteration)
             prevIteration <- self$mapIter[[as.character(as.numeric(iteration-1))]]

             # hptype = "A"
             foreach (hptype = hptypes) %do% {
               reftag <- prevIteration[[hptype]]$ref
               outdir   <- prevIteration[[hptype]]$dir
               readpath <- prevIteration[[hptype]]$reads
               refpath  <- prevIteration[[hptype]]$seqpath

               nreads <- if (!is.null(attr(self$getHapList(hptype), "n"))) {
                 attr(self$getHapList(hptype), "n")
               } else {
                 length(self$getHapList(hptype))
               }
               optsname <- sprintf("%s [%s]", hptype, nreads)
               mapfmt  <- "mapIter <%s> <%s> <%s> <%s>"
               maptag  <- sprintf(mapfmt, iteration, hptype, self$getLrdType(), self$getMapper(), optstring(opts, optsname))

               self$mapIter[[as.character(iteration)]][[hptype]] = structure(
                 list(
                   dir     = outdir,
                   reads   = readpath,
                   ref     = refpath,
                   bamfile = list(),
                   pileup  = list(),
                   conseq  = list(),
                   seqpath = list(),
                   params  = list(pruning_cutoff = pruning_cutoff),
                   tag     = list()
                 ),
                 class = c("mapIter", "list")
               )

                # self$mapIter[["0"]][[hptype]]$ref <- "multialign"
                # self$mapIter[["0"]][[hptype]]$conseq <- cons
                # self$mapIter[["0"]][[hptype]]$seqpath <- consout
                # self$mapIter[["0"]][[hptype]]$tag <- "multialign"

               ## Run mapper
               flog.info("  Mapping of partitioned longreads of haplotype %s ...", hptype, name = "info")
               samfile <- map_fun(
                 reffile  = refpath,
                 readfile = readpath,
                 allele   = paste0("mapIter", iteration),
                 readtype = self$getLrdType(),
                 opts     = opts,
                 refname  = reftag,
                 optsname = optsname,
                 force    = force,
                 outdir   = outdir
               )
#
#
#                  reffile  = refpath
#                  readfile = readpath
#                  allele   = paste0("mapIter", iteration)
#                  readtype = self$getLrdType()
#                  opts     = opts
#                  refname  = reftag
#                  optsname = optsname
#                  force    = force
#                  outdir   = outdir

               # samfile
               ## Run bam - sort - index pipeline
               flog.info("   Indexing ...", name = "info")
               bamfile <- bam_sort_index(
                 samfile,
                 refpath,
                 pct / 100,
                 min_mapq,
                 force = force,
                 clean = FALSE
                 # debug
                 # clean = TRUE
               )
               self$mapIter[[iterationC]][[hptype]]$bamfile = bamfile

               ## Calculate pileup from graphmap produced SAM file
               flog.info("   Piling up ...", name = "info")
               pileup <- Pileup(
                 bamfile,
                 self$getThreshold(),
                 max_depth = max_depth,
                 min_base_quality = min_base_quality,
                 min_mapq = min_mapq,
                 min_nucleotide_depth = min_nucleotide_depth,
                 include_deletions = include_deletions,
                 include_insertions = include_insertions
               )
               if (include_insertions && is.null(ins(pileup$consmat))) {
                 pileup <- pileup_include_insertions(x = pileup, threshold = 0.15)
               }
               self$mapIter[[iterationC]][[hptype]]$pileup = pileup

               # ## Construct consensus sequence
               flog.info("   Constructing a consensus ...", name = "info")
               conseq_name <- paste0("consensus.", sub(".sam.gz", "", basename(samfile)))
               conseq      <- conseq(x = pileup, name = conseq_name, type = "prob",
                                     exclude_gaps = TRUE, prune_matrix = TRUE,
                                     cutoff = pruning_cutoff)
               seqpath     <- file.path(outdir, paste0(conseq_name, ".fa"))
               self$mapIter[[iterationC]][[hptype]]$seqpath = seqpath
               self$mapIter[[iterationC]][[hptype]]$conseq = conseq
               self$mapIter[[iterationC]][[hptype]]$ref = conseq_name
               Biostrings::writeXStringSet(
                 Biostrings::DNAStringSet(gsub("[-+]", "", conseq)),
                 seqpath
               )

               ## Set maptag
               self$mapIter[[iterationC]][[hptype]]$tag = maptag
             }

             if (plot) {
               flog.info(" Plotting MapIter Summary ...", name = "info")
               ## Coverage and base frequency
               gfile <- file.path(
                 self$getOutdir(),
                 paste("plot.MapIter", iteration, self$getLrdType(), self$getMapper(), "pdf", sep = ".")
               )
               pdf(file = gfile, width = 8 * length(hptypes), height = 8, onefile = TRUE, title = paste(self$getLocus(), self$getSampleId(), sep = "." ))
               self$plotmapIterSummary(thin = 0.1, width = 4, iteration = iteration, drop.indels = TRUE)
               ## map1: 0.1; 4
               dev.off()
               ## Consensus sequence probability
               gfile <- file.path(
                 self$getOutdir(),
                 paste("plot.MapIter.conseq", iteration, self$getLrdType(), self$getMapper(), "pdf", sep = ".")
               )
               pdf(file = gfile, width = 16, height = 4*length(hptypes), onefile = TRUE, title = paste(self$getLocus(), self$getSampleId(), sep = "." ))
               self$plotmapIterSummaryConseqProb(text_size = 1.75,
                                              iteration = iteration,
                                              point_size = 0.75,
                                              threshold = "auto")
               # text_size = 1.75
               #                                iteration = iteration
               #                                point_size = 0.75
               #                                threshold = "auto"
               dev.off()
             }
           }

           invisible(self)
         })

#' @export

print.mapIter <- function(x, ...) {
  msg <- sprintf("An object of class '%s'.\\n", class(x)[1])
  bamf <- ifelse(is.null(x$bamfile), " no bamfile", basename(x$bamfile %||% ""))
  msg <- sprintf(
    "%s [Dir] %s\\n [Reads] %s\\n [Reference] %s\\n [Bamfile] %s\\n",
    msg,
    x$dir,
    basename(x$reads),
    basename(x$ref),
    bamf
  )
  cat(msg)
}



## Method: partition Shortreads ####

#' @export
partitionShortReads.DR2S <- function(x,
                      opts = list(),
                      force = FALSE) {
  flog.info("Step 4: Partition shortreads based on initial mapping and longread clustering ... ", name = "info")
  Sys.sleep(1)
  x$runPartitionShortReads(opts = opts,
            force = force)
  message("  Done!\n")
  invisible(x)
}
#self <- dpb1_3.map
DR2S_$set("public", "runPartitionShortReads",
         function(opts = list(),
                  force = FALSE,
                  optsname = "",
                  pct = 100,
                  threshold = 0.20,
                  min_mapq = 0) {


                  ## debug
                    # opts = list()
                    # force = FALSE
                    # optsname = ""
                    # pct = 100
                    # threshold = 0.20
                    # min_mapq = 0

           ## Overide default arguments
           args <- self$getOpts("partitionSR")
           if (!is.null(args)) {
             env  <- environment()
             list2env(args, envir = env)
           }

           ## stop if no shortreads provided
           if (is.null(self$getConfig("shortreads"))) {
             flog.warn(" Cannot partition shortreads. No shortreads provided", name = "info")
             return(invisible(self))
           }

           if (self$getPartSR()){
             flog.info(" Found shortread mapping from MapInit ...", name = "info")
             bamfile <- self$mapInit$SR2$bamfile
           } else {
             mapfmt  <- "mapPartSR <%s> <%s> <%s> <%s>"
             maptag  <- sprintf(mapfmt, self$mapInit$SR1$ref, self$getSrdType(), self$getMapper(), optstring(opts, optsname))
             if (threshold != self$getThreshold()) {
               self$setThreshold(threshold)
             }
             ## Fetch mapper
             map_fun <- self$getMapFun()
             ## Run mapper
             flog.warn(" Found no shortread mapping from MapInit ...", name = "info")
             flog.info(" Mapping short reads against provided reference ...", name = "info")
             samfile <- map_fun(
               reffile  = self$mapInit$SR1$seqpath,
               readfile = self$getShortreads(),
               allele   = self$mapInit$SR1$ref,
               readtype = self$getSrdType(),
               opts     = opts,
               refname  = "",
               optsname = optsname,
               force    = force,
               outdir   = self$getOutdir()
             )
             ## Run bam - sort - index pipeline
             flog.info("  Indexing ...", name = "info")
             bamfile <- bam_sort_index(
               samfile = samfile,
               reffile = self$mapInit$SR1$seqpath,
               sample = pct / 100,
               min_mapq = min_mapq,
               force = force,
               # clean = TRUE
               # debug
               clean = FALSE
             )
           }
           ref <- self$mapInit$SR1$conseq
           refname <- names(ref)
           hptypes <- self$getHapTypes()

           prt_mat <- self$partition$mat
           seqs <- lapply(self$partition$hpl, function(x) get_seqs_from_mat(as.matrix(prt_mat[x,])))
           names(seqs) <- hptypes

           mats <- lapply(seqs, function(x) create_PWM(x))
           mats <- foreach(m = mats) %do% {
             colnames(m) <- colnames(prt_mat)
             m
           }
           names(mats) <- names(seqs)

           ppos <- colnames(mats[[1]])


           # Run partitioning
           srpartition <- get_SR_partition_scores(ppos, refname, bamfile, mats, cores = "auto")

           ##' Assign read to haplotype with highest probability,
           ##' i.e. prod over probabilities of each haplotype and choose max
           flog.info(" Get highest scoring haplotype for each read", name = "info")
           srpartition$haplotypes <- score_highest_SR(srpartition$srpartition, diffThreshold = 0.001)

           # Write fastqs
           foreach(hptype = hptypes ) %do% {
             srfilenames <- c()
             flog.info(" Write shortread fastq for haplotype %s ...", hptype, name = "info")
             fqs <- self$getShortreads()
             dontUseReads <- srpartition$haplotypes$read[!srpartition$haplotypes$read %in% dplyr::filter(srpartition$haplotypes, haplotype == hptype)$read]
             #fq <- fqs[1]
             # write fastq's
             foreach(fq = fqs) %do% {
               srFastqHap = file.path(self$getOutdir(), hptype, paste0(c(strsplit(basename(fq), "\\.")[[1]][1], hptype, "fastq.gz"), collapse = "."))
               write_part_fq(fq = fq, srFastqHap = srFastqHap, dontUseReads = dontUseReads)
               srfilenames <- c(srfilenames, srFastqHap)
               self$srpartition[[hptype]]$srpartition <- srpartition
            }
            self$srpartition[[hptype]]$SR[[fq]] <- srfilenames
          }
           return(invisible(self))
         })

#' @export
print.partitionShortReads <- function(x, ...) {
  msg  <- sprintf("An object of class '%s'.\n", class(x)[1])
  bamf <- paste0(basename(unlist(x$bamfile) %||% ""), collapse = ", ")
  seqp <- paste0(basename(unlist(x$seqpath) %||% ""), collapse = ", ")

  msg <- sprintf(
    "%s [Dir] %s\n [Longreads] %s\n [Shortreads] %s\n [References] %s\n [Bamfile] %s\n [Seqpath] %s\n",
    msg, x$dir,
    paste0(basename(unlist(x$reads)), collapse = ", "),
    paste0(basename(x$sreads), collapse = ", "),
    paste0(basename(unlist(x$ref)), collapse = ", "),
    bamf, seqp
  )
  cat(msg)
}


## Method: mapFinal ####

#' @export
mapFinal.DR2S <- function(x,
                      opts = list(),
                      pct = 100,
                      min_base_quality = 3,
                      min_mapq = 0,
                      max_depth = 1e5,
                      min_nucleotide_depth = 3,
                      include_deletions = TRUE,
                      include_insertions = TRUE,
                      include_read_ids = TRUE,
                      force = FALSE,
                      fullname = TRUE,
                      plot = TRUE,
                      clip = TRUE) {
  flog.info("Step 5: Mapping short and long reads against refined consensus sequences ...", name = "info")
  Sys.sleep(1)
  x$runMapFinal(opts = opts, pct = pct, min_base_quality = min_base_quality,
            min_mapq = min_mapq, max_depth = max_depth,
            min_nucleotide_depth = min_nucleotide_depth,
            include_deletions = include_deletions,
            include_insertions = include_insertions,
            include_read_ids = include_read_ids, force = force,
            fullname = fullname, plot = plot, clip = clip)
  message("  Done!\n")
  invisible(x)
}

# debug
# self <- dpb1_3.map
# opts = list()
# pct = 100
# min_base_quality = 3
# min_mapq = 0
# max_depth = 1e5
# min_nucleotide_depth = 3
# include_deletions = TRUE
# include_insertions = TRUE
# include_read_ids = TRUE
# force = FALSE
# fullname = TRUE
# plot = TRUE
# clip = TRUE

DR2S_$set("public", "runMapFinal",
         function(opts = list(),
                  pct = 100,
                  min_base_quality = 3,
                  min_mapq = 0,
                  max_depth = 1e5,
                  min_nucleotide_depth = 3,
                  include_deletions = TRUE,
                  include_insertions = TRUE,
                  include_read_ids = TRUE,
                  force = FALSE,
                  fullname = TRUE,
                  plot = TRUE,
                  clip = TRUE) {

           ## Overide default arguments
           args <- self$getOpts("mapFinal")
           if (!is.null(args)) {
             env  <- environment()
             list2env(args, envir = env)
           }

           ## stop if no shortreads provided
           if (is.null(self$getConfig("shortreads"))) {
             flog.warn(" Cannot run mapFinal. No shortreads provided", name = "info")
             return(invisible(self))
           }

           reftag    <- "final"
           outdir    <- dir_create_if_not_exists(file.path(self$getOutdir(), reftag))

           lastIter  <- self$mapIter[[max(names(self$mapIter))]]
           hptypes   <- self$getHapTypes()
           readpathsLR  <- sapply(hptypes, function(x) lastIter[[x]]$reads)
           names(readpathsLR) <- hptypes
           refpaths   <- sapply(hptypes, function(x) lastIter[[x]]$seqpath)
           names(refpaths) <- hptypes
           readpathsSR <- sapply(hptypes, function(x) self$srpartition[[x]]$SR)
           names(readpathsSR) <- hptypes
           ## Mapper
           map_fun <- self$getMapFun()

           self$mapFinal = structure(
             list(
               dir     = outdir,
               sreads  = readpathsSR,
               lreads   = readpathsLR,
               ref     = refpaths,
               bamfile = list(),
               pileup  = list(),
               tag     = list(),
               seqpath = list()
             ), class = c("mapFinal", "list")
           )

           ## Remap long reads to the same reference sequences as short reads
           # debug
           # hptype = "A"
           for (hptype in hptypes) {
             flog.info(" Run mapFinal for haplotype %s ...", hptype, name = "info" )
             refpath  <- refpaths[[hptype]]

             mapgroupLR <- paste0("LR", hptype)
             maptagLR   <- paste("mapFinal", mapgroupLR, self$getLrdType(), self$getMapper(),
                               optstring(opts), sep = ".")
             readpathLR <- readpathsLR[[hptype]]
             ## Run mapper
             flog.info("   Mapping long reads against latest consensus ...", name = "info")
             samfile <- map_fun(
               reffile  = refpath,
               readfile = readpathLR,
               allele   = mapgroupLR,
               readtype = self$getLrdType(),
               opts     = opts,
               refname  = hptype,
               optsname = optstring(opts),
               force    = force,
               outdir   = outdir
             )

             ## Run bam - sort - index pipeline
             flog.info("  Indexing ...", name = "info")
             bamfile <- bam_sort_index(
               samfile,
               refpath,
               pct / 100,
               min_mapq,
               force = force,
               clean = FALSE
               # debug
               # clean = TRUE
             )
             self$mapFinal$bamfile[[mapgroupLR]] = bamfile

             ## Calculate pileup from graphmap produced SAM file
             flog.info("  Piling up ...", name = "info")
             pileup <- Pileup(
               bamfile,
               self$getThreshold(),
               max_depth = max_depth,
               min_base_quality = min_base_quality,
               min_mapq = min_mapq,
               min_nucleotide_depth = min_nucleotide_depth,
               include_deletions = include_deletions,
               include_insertions = FALSE
             )
             self$mapFinal$pileup[[mapgroupLR]] = pileup

             ## Set maptag
             self$mapFinal$tag[[mapgroupLR]] = maptagLR


           ## Map short reads
             if (!is.null(readpathsSR[[hptype]])){
               mapgroupSR <- paste0("SR", hptype)
               maptagSR   <- paste("mapFinal", mapgroupSR, self$getLrdType(), self$getMapper(),
                                 optstring(opts), sep = ".")

               readfiles <- readpathsSR[[hptype]]
               # debug:
               #readfiles <- self$getShortreads()

               ## Run mapper
               flog.info("  Mapping short reads against latest consensus ...", name = "info")

               samfile <- map_fun(
                 reffile  = refpath,
                 readfile = readfiles,
                 # if we run shortreads against both pacbio and nanopore data
                 # this hack makes sure that we can distinguish the bam files ->
                 # we get pacbio.illumina.bwamem.A...bam and nanopore.illumina.bwamem.A...bam
                 allele   = paste0(mapgroupSR, ".", self$getLrdType()),
                 readtype = self$getSrdType(),
                 opts     = opts,
                 refname  = hptype,
                 optsname = optstring(opts),
                 force    = force,
                 outdir   = outdir
               )

               if (clip) {
                 flog.info(" Trimming softclips and polymorphic ends ...", name = "info")
                 ## Run bam - sort - index pipeline
                 bamfile <- bam_sort_index(samfile, refpath, pct / 100, min_mapq, force = force, clean = TRUE)
                 ## Trim softclips
                 fq <- trim_softclipped_ends(bam = Rsamtools::scanBam(bamfile)[[1]], preserve_ref_ends = TRUE)
                 ## Trim polymorphic ends
                 fq <- trim_polymorphic_ends(fq)
                 ## Write new shortread file to disc
                 fqdir  <- dir_create_if_not_exists(file.path(self$getOutdir(), "final"))
                 fqfile <- paste("sread", hptype, self$getMapper(), "trimmed", "fastq", "gz", sep = ".")
                 fqout  <- file_delete_if_exists(file.path(fqdir, fqfile))
                 ShortRead::writeFastq(fq, fqout, compress = TRUE)
                 file_delete_if_exists(bamfile)
                 ## Rerun mapper
                 flog.info("  Mapping trimmed short reads against latest consensus ... ", name = "info")
                 samfile <- map_fun(
                   reffile  = refpath,
                   readfile = fqout,
                   allele   = paste0(mapgroupSR, ".", self$getLrdType()),
                   readtype = self$getSrdType(),
                   opts     = list(A = 1, B = 4, O = 2),
                   refname  = hptype,
                   optsname = optstring(opts),
                   force    = force,
                   outdir   = outdir
                 )
                 # cleanup
                 file_delete_if_exists(fqout)
               }

               ## Run bam - sort - index pipeline
               flog.info("  Indexing ...", name = "info")
               bamfile <- bam_sort_index(
                 samfile,
                 refpath,
                 pct / 100,
                 min_mapq,
                 force = force,
                 clean = FALSE

                 # debug
                 # clean = TRUE
               )
               self$mapFinal$bamfile[[mapgroupSR]] = bamfile

               ## Calculate pileup from graphmap produced SAM file
               flog.info("  Piling up ...", name = "info")
               pileup <- Pileup(
                 bamfile,
                 self$getThreshold(),
                 max_depth = max_depth,
                 min_base_quality = min_base_quality + 10,
                 min_mapq = min_mapq,
                 min_nucleotide_depth = min_nucleotide_depth,
                 include_deletions = include_deletions,
                 include_insertions = include_insertions
               )
               if (include_read_ids && is.null(ids(pileup$consmat))) {
                 pileup <- pileup_include_read_ids(pileup)
               }

               if (include_insertions && is.null(ins(pileup$consmat))) {
                 pileup <- pileup_include_insertions(pileup)
               }
               self$mapFinal$pileup[[mapgroupSR]] = pileup

               ## Set maptag
               self$mapFinal$tag[[mapgroupSR]] = maptagSR
             }

             # calc new consensus
             cseq <- conseq(pileup$consmat, paste0("mapFinal", hptype), "ambig", exclude_gaps = TRUE, threshold = self$getThreshold())
             self$mapFinal$seq[[hptype]] <- cseq
           }


           if (plot) {
             flog.info(" Plotting final summary figures ...", name = "info")
             ## Coverage and base frequency

             if (!is.null(self$mapFinal$sreads$A)){
               for (readtype in c("LR", "SR")){
                 gfile <- file.path(self$getOutdir(), paste("plotFinal", readtype, self$getLrdType(), self$getMapper(), "pdf", sep = "."))
                 pdf(file = gfile, width = 8 * length(hptypes), height = 8, onefile = TRUE, title = paste(self$getLocus(), self$getSampleId(), sep = "." ))
                 self$plotmapFinalSummary(iteration = "final", readtype = readtype, thin = 0.25, width = 20)
                 dev.off()
               }
             } else {
               gfile <- file.path(self$getOutdir(), paste("plotFinal.LR", self$getLrdType(), self$getMapper(), "pdf", sep = "."))
               pdf(file = gfile, width = 8 * length(hptypes), height = 8, onefile = TRUE, title = paste(self$getLocus(), self$getSampleId(), sep = "." ))
               self$plotmapFinalSummary(iteration = "final", readtype = "LR", thin = 0.25, width = 20)
               dev.off()
             }
           }

           return(invisible(self))
         })

#' @export
print.mapFinal <- function(x, ...) {
  msg  <- sprintf("An object of class '%s'.\n", class(x)[1])
  bamf <- paste0(basename(unlist(x$bamfile) %||% ""), collapse = ", ")
  msg <- sprintf(
    "%s [Dir] %s\n [Longreads] %s\n [Shortreads] %s\n [References] %s\n [Bamfile] %s",
    msg, x$dir,
    paste0(basename(unlist(x$lreads)), collapse = ", "),
    paste0(basename(unlist(x$sreads)), collapse = ", "),
    paste0(basename(unlist(x$ref)), collapse = ", "),
    bamf
  )
  cat(msg)
}

## Method: runPipeline ####

DR2S_$set("public", "runPipeline", function() {
  steps_ <- self$getPipeline()
  while (length(steps_) > 0) {
    step <- steps_[1]
    private$run_(step)
    steps_ <- steps_[-1]
  }
  self
})


#' ## Method:  polish ####
#'
#' #' @export
#' polish.DR2S <- function(x) {
#'   flog.info("Step 6: Infer problematic positions and consensus sequence ...", name = "info")
#'   Sys.sleep(1)
#'   x$polish()
#'   message("  Done!\n")
#'   invisible(x)
#' }
#'
#' DR2S_$set("public", "polish", function() {
#'   self <- polish(self)
#'   return(invisible(self))
#' }
#' ## Method: report  ####
#'
#' #' @export
#' mapFinal.DR2S <- function(x) {
#'   flog.info("Step 7: report consensus sequences and problematic positions ...", name = "info")
#'   Sys.sleep(1)
#'   x$report()
#'   message("  Done!\n")
#'   invisible(x)
#' }
#'
#' DR2S_$set("public", "report", function() {
#'   self <- report(self)
#'   return(invisible(self))
#' }
