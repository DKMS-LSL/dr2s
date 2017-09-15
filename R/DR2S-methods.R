# Method: map0 ####

#' @export
map0.DR2S <- function(x,
                      opts = list(),
                      optsname = "",
                      pct = 100,
                      threshold = 0.20,
                      min_base_quality = 3,
                      min_mapq = 0,
                      max_depth = 1e4,
                      min_nucleotide_depth = 3,
                      include_deletions = FALSE,
                      include_insertions = FALSE,
                      force = FALSE,
                      fullname = TRUE,
                      # for now set plot to FALSE; cant plot with the clustering
                      plot = FALSE) {
                      #plot = TRUE) {
  message("\nStep 0: Initial mapping of long reads ... \n")
  Sys.sleep(1)
  x$runMap0(opts = opts, optsname = optsname, pct = pct, threshold = threshold,
            min_base_quality = min_base_quality, min_mapq = min_mapq,
            max_depth = max_depth, min_nucleotide_depth = min_nucleotide_depth,
            include_deletions = include_deletions, include_insertions = include_insertions,
            force = force, fullname = fullname, plot = plot)
  message("  Done!\n")
  invisible(x)
}

DR2S_$set("public", "runMap0",
         function(opts = list(),
                  optsname = "",
                  pct = 100,
                  threshold = 0.20,
                  min_base_quality = 3,
                  min_mapq = 0,
                  max_depth = 1e4,
                  min_nucleotide_depth = 3,
                  include_deletions = FALSE,
                  include_insertions = FALSE,
                  force = FALSE,
                  fullname = TRUE,
                  plot = TRUE) {
           stopifnot(pct > 0 && pct <= 100)
           recode_header <- FALSE

           ## Overide default arguments
           args <- self$getOpts("map0")
           if (!is.null(args)) {
             env  <- environment()
             list2env(args, envir = env)
           }

           if (recode_header) {
             stopifnot(
               recode_fastq_header(self$getLongreads()),
               recode_fastq_header(self$getShortreads()[1]),
               recode_fastq_header(self$getShortreads()[2])
             )
           }

           mapfmt  <- "map0 <%s> <%s> <%s> <%s>"
           maptag  <- sprintf(mapfmt, self$getReference(), self$getLrdType(), self$getMapper(), optstring(opts, optsname))
           if (threshold != self$getThreshold()) {
             self$setThreshold(threshold)
           }
           ## Fetch mapper
           map_fun <- self$getMapFun()
           ## Run mapper
           message("  Mapping ...\n")
           samfile <- map_fun(
             reffile  = self$getRefPath(),
             readfile = self$getLongreads(),
             allele   = self$getReference(),
             readtype = self$getLrdType(),
             opts     = opts,
             refname  = "",
             optsname = optsname,
             force    = force,
             outdir   = self$getOutdir()
           )

           ## Run bam - sort - index pipeline
           message("  Indexing ...")
           bamfile <- bam_sort_index(
             samfile = samfile,
             reffile = self$getRefPath(),
             sample = pct / 100,
             min_mapq = min_mapq,
             force = force,
             clean = TRUE
           )

           ## Calculate pileup from graphmap produced SAM file
           message("  Piling up ...")
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

           self$map0 = structure(
             list(
               reads   = self$getLongreads(),
               bamfile = bamfile,
               pileup  = pileup,
               tag     = maptag
             ),
             class  = c("map0", "list")
           )

           if (plot) {
             message("  Plotting ...")
             ## Coverage and frequency of minor alleles
             gfile <- file.path(self$getOutdir(), paste0("plot0.", sub("bam$", "pdf", usc(basename(bamfile)))))
             pdf(file = gfile, width = 12, height = 10, onefile = TRUE)
             self$plotMap0Summary(
               thin = 0.25,
               width = 2,
               label = maptag
             )
             dev.off()
           }

           return(invisible(self))
         })

#' @export
print.map0 <- function(x, ...) {
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
                                      # shuffle = TRUE,
                                      skip_gap_freq = 2/3,
                                      plot = FALSE,
                                      # plot = TRUE,
                                      ...) {
  message("\nStep 1: Partition reads into haplotypes ...")
  Sys.sleep(1)
  x$runHaplotypePartitioning(max_depth = max_depth,
                             # shuffle = shuffle,
                             skip_gap_freq = skip_gap_freq,
                             plot = plot)
  message("  Done!\n")
  invisible(x)

}

#self <- dedk
DR2S_$set("public", "runHaplotypePartitioning", function(max_depth = 1e4,
                  skip_gap_freq = 2/3,
                  plot = TRUE) {

                  # debug
                  # max_depth = 1e4
                  # shuffle = TRUE,
                  # skip_gap_freq = 2/3
                  # plot = TRUE
           stopifnot(self$hasPileup())

           ## Overide default arguments
           args <- self$getOpts("partition")
           if (!is.null(args)) {
             env  <- environment()
             list2env(args, envir = env)
           }

           tag  <- self$getMapTag(0)
           th   <- self$getThreshold()
           ppos <- self$polymorphicPositions()
           message("  Generating SNP matrix ...")
           mat <- if (tryCatch(
             !is(self$partition, "PartList"),
             error = function(e)
               TRUE
           ) ||
           !(all(ppos$position %in% colnames(self$partition$mat)) &&
             all(colnames(self$partition$mat) %in% ppos$position))) {
             DR2S:::SNPmatrix(bamfile = self$map0$bamfile, max_depth, polymorphic_positions = ppos)
           } else {
             self$partition$mat
           }

           message(sprintf("  Partitioning %s reads over %s SNPs ...", NROW(mat), NCOL(mat)))
           prt <- partition_reads(x = mat, skip_gap_freq = skip_gap_freq)
           self$setHapTypes(levels(as.factor(PRT(self$partition$prt))))

           self$partition = structure(list(
             mat = mat,
             prt = prt,
             hpl = NULL,
             lmt = NULL
           ),
           class = c("PartList", "list"))

           if (FALSE){
           # if (plot) {
             message("  Plotting ...")
             bnm   <- basename(self$map0$bamfile)
             outf  <- file.path(self$getOutdir(), paste0("plot0.partition.", sub("bam$", "pdf", usc(bnm))))
             pdf(file = outf, width = 10, height = 12, onefile = TRUE)
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
  x$splitReadsByHaplotype( limits)
  invisible(x)
}

#self <- dl3.part
# x <- dl3.part

DR2S_$set("public", "splitReadsByHaplotype", function(limits,
                                                      plot = FALSE){
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
           reads <- lapply(names(self$getLimits()), function(x) dplyr::filter(prt, haplotype == x, mcoef <= self$getLimits()[x]))
           names(reads) <- names(self$getLimits())

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

           if (FALSE) {
           # if (plot) {
             tag   <- self$getMapTag(0)
             bnm   <- basename(self$map0$bamfile)
             outf  <- file.path(self$getOutdir(), paste0("plot0.partition.", sub("bam$", "pdf", usc(bnm))))
             pdf(file = outf, width = 10, height = 12, onefile = TRUE)
             self$plotPartitionSummary(label = tag, limits = c(self$getLimitA(), self$getLimitB()))
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
  message("\nStep2: Extracting haplotyped reads ...")
  Sys.sleep(1)
  x$extractFastq(nreads = nreads, replace = replace, nalign = nalign)
  message("  Done!\n")
  invisible(x)
}
#self <- hla.split
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
              ## extract haplotype A reads from bamfile
              fq  <- extract_fastq(
                x = self$map0$bamfile,
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
              self$map1[[hptype]] = structure(
                list(
                  dir     = dir,
                  reads   = out,
                  ref     = NULL,
                  bamfile = list(reference = NULL, alternate = NULL),
                  pileup  = list(reference = NULL, alternate = NULL),
                  conseq  = list(
                    reference  = NULL,
                    alternate  = NULL,
                    merged     = NULL,
                    multialign = NULL
                  ),
                  seqpath = list(
                    reference  = NULL,
                    alternate  = NULL,
                    merged     = NULL,
                    multialign = NULL
                  ),
                  tag     = list(
                    reference  = NULL,
                    alternate  = NULL,
                    merged     = NULL,
                    multialign = NULL
                  )
                ),
                class = c("map1", "list")
              )
              if (self$hasMultialign()) {
                message("   Constructing initial ", hptype, " consensus from ", nalign, " reads ...")
                cons <- multialign_consensus(multialign(self, hptype, nalign))
                hptype
                consfile <- paste("consensus", self$getLrdType(), "multialign", paste0(hptype, nalign), "fa", sep = ".")
                names(cons) <- consfile
                consout <- file_delete_if_exists(file.path(dir, consfile))
                Biostrings::writeXStringSet(cons, consout)
                self$map1[[hptype]]$ref <- list(refpath = consout, refseq = cons)
              }
            }

            return(invisible(self))
          })

## Method: map1 ####

#' @export
map1.DR2S <- function(x,
                      opts = list(),
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
  message("\nStep 3: Mapping partitioned long reads against the reference ... \n")
  Sys.sleep(1)
  x$runMap1(opts = opts, pct = pct, min_base_quality = min_base_quality,
            min_mapq = min_mapq, max_depth = max_depth,
            min_nucleotide_depth = min_nucleotide_depth,
            include_deletions = include_deletions,
            include_insertions = include_insertions,
            pruning_cutoff = pruning_cutoff, force = force,
            fullname = fullname, plot = plot)
  message("  Done!\n")
  invisible(x)
}

# debug
#self <- dedk.fq
#self$map1$
DR2S_$set("public", "runMap1",
         function(opts = list(),
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

                  # debug
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
           stopifnot(pct > 0 && pct <= 100)

           ## Overide default arguments
           args <- self$getOpts("map1")
           if (!is.null(args)) {
             env  <- environment()
             list2env(args, envir = env)
           }

           ## Mapper
           map_fun <- self$getMapFun()

           ## Probably no need to use anymore as we can use multialign
           if(FALSE){
           refs <- if (!is.null(self$getAltPath())) {
             c("reference", "alternate")
           } else if (self$hasMultialign()) {
             "multialign"
           } else {
             "reference"
           }
           # group = "A"
           # reftag = "multialign"
           }
           hptypes <- levels(as.factor(PRT(self$partition$prt)))
           # foreach(group = hptypes) %:%
           #   foreach(reftag = refs) %do% {
           reftag <- "multialign"

           foreach(group = hptypes) %do% {
               if (is.null(outdir <- self$map1[[group]]$dir)) {
                 stop("Readfiles missing for readgroup ", dQuote(group))
               }
               nreads <- if (!is.null(attr(self$getHapList(group), "n"))) {
                 attr(self$getHapList(group), "n")
               } else {
                 length(self$getHapList(group))
               }
               optsname <- sprintf("%s [%s]", group, nreads)
               allele  <- switch(reftag,
                                 reference  = self$getReference(),
                                 alternate  = self$getAlternate(),
                                 multialign = "multialign")
               mapfmt  <- "map1 <%s> <%s> <%s> <%s>"
               maptag  <- sprintf(mapfmt, allele, self$getLrdType(), self$getMapper(), optstring(opts, optsname))

               ## Run mapper
               message("  Mapping allele group ", dQuote(group), " against ", dQuote(allele),  " ...")
               samfile <- map_fun(
                 reffile  = switch(reftag, reference = self$getRefPath(), alternate = self$getAltPath(), multialign = self$map1[[group]]$ref$refpath),
                 readfile = self$map1[[group]]$reads,
                 allele   = "map1",
                 readtype = self$getLrdType(),
                 opts     = opts,
                 refname  = reftag,
                 optsname = optsname,
                 force    = force,
                 outdir   = outdir
               )

               ## Run bam - sort - index pipeline
               message("  Indexing ...")
               bamfile <- bam_sort_index(
                 samfile,
                 switch(reftag, reference = self$getRefPath(), alternate = self$getAltPath(), multialign = self$map1[[group]]$ref$refpath),
                 pct / 100,
                 min_mapq,
                 force = force,
                 clean = TRUE
               )
               self$map1[[group]]$bamfile[[reftag]] = bamfile

               ## Calculate pileup from graphmap produced SAM file
               message("  Piling up ...")
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
               self$map1[[group]]$pileup[[reftag]] = pileup

               ## Construct consensus sequence
               message("  Constructing a consensus ...")
               conseq_name <- paste0("consensus.", sub(".sam.gz", "", basename(samfile)))
               conseq      <- conseq(x = pileup, name = conseq_name, type = "prob",
                                     exclude_gaps = TRUE, prune_matrix = TRUE,
                                     cutoff = pruning_cutoff)
               seqpath     <- file.path(outdir, paste0(conseq_name, ".fa"))
               self$map1[[group]]$conseq[[reftag]]  = conseq
               self$map1[[group]]$seqpath[[reftag]] = seqpath
               Biostrings::writeXStringSet(
                 Biostrings::DNAStringSet(gsub("[-+]", "", conseq)),
                 seqpath
               )

               ## Set maptag
               self$map1[[group]]$tag[[reftag]] = maptag
             }


           ## NO need to use as we use multialign
           if (FALSE){
           ## infer which haplogroup is closer to the reference allele
           if (length(refs) == 1) {
             nA <- n_polymorphic_positions(self$map1$A$pileup[[reftag]]$consmat)
             nB <- n_polymorphic_positions(self$map1$B$pileup[[reftag]]$consmat)
             if (nA < nB) {
               self$setARefType("ref")
               self$setBRefType("alt")
             } else {
               self$setARefType("alt")
               self$setBRefType("ref")
             }
           } else if (length(refs) == 2) {
             nAref <- n_polymorphic_positions(self$map1$A$pileup$reference$consmat)
             nAalt <- n_polymorphic_positions(self$map1$A$pileup$alternate$consmat)
             nBref <- n_polymorphic_positions(self$map1$B$pileup$reference$consmat)
             nBalt <- n_polymorphic_positions(self$map1$B$pileup$alternate$consmat)
             if (nAref < nAalt && nBref > nBalt) {
               self$setARefType("ref")
               self$setBRefType("alt")
             } else {
               self$setARefType("alt")
               self$setBRefType("ref")
             }
           }

           ## merge consensus sequences
           if (length(refs) == 2) {
             message("  Merging the consensus sequences ...")

             mseqpathA <- sub("\\.reference\\.", ".merged.", self$map1$A$seqpath$reference)
             mconseqA <- self$mergeMap1Conseqs(group = "A")
             self$map1$A$conseq$merged = mconseqA
             self$map1$A$seqpath$merged = mseqpathA
             Biostrings::writeXStringSet(mconseqA, mseqpathA)
             rtag <- strsplit(self$map1$A$tag$reference, " ")[[1]]
             atag <- strsplit(self$map1$A$tag$alternate, " ")[[1]]
             mtag <- intersect(rtag, atag)
             self$map1$A$tag$merged <-
               paste0(c(mtag[1],
                        c(setdiff(rtag, atag), setdiff(atag, rtag)),
                        mtag[2:length(mtag)]), collapse = " ")

             mseqpathB <- sub("\\.reference\\.", ".merged.", self$map1$B$seqpath$reference)
             mconseqB <- self$mergeMap1Conseqs("B")
             self$map1$B$conseq$merged = mconseqB
             self$map1$B$seqpath$merged = mseqpathB
             Biostrings::writeXStringSet(mconseqB, mseqpathB)
             rtag <- strsplit(self$map1$B$tag$reference, " ")[[1]]
             atag <- strsplit(self$map1$B$tag$alternate, " ")[[1]]
             mtag <- intersect(rtag, atag)
             self$map1$B$tag$merged <-
               paste0(c(mtag[1],
                        c(setdiff(rtag, atag), setdiff(atag, rtag)),
                        mtag[2:length(mtag)]), collapse = " ")
           }
           }


           ## Work here:
           # First clean error
           ##self <- dedk.map1
           ##hptypes <- c("A", "B")
           #self$map1$A <- self$map1$C
           #self$map1$B <- self$map1$D

           if (plot) {
             message("  Plotting ...")
             if (self$hasMap1Alternate()) {
               foreach(group = hptypes) %do% {
                 ## Coverage and base frequency
                 gfile <- file.path(
                   self$getOutdir(),
                   paste("plot1", group, self$getLrdType(), self$getMapper(), "pdf", sep = ".")
                 )
                 pdf(file = gfile, width = 16, height = 8, onefile = TRUE)
                 self$plotMap1Summary(group = group, thin = 0.1, width = 4)
                 dev.off()
               }
             } else {
               gfile <- file.path(
                 self$getOutdir(),
                 paste("plot1", self$getLrdType(), self$getMapper(), "pdf", sep = ".")
               )
               pdf(file = gfile, width = 16, height = 8, onefile = TRUE)
               self$plotMap1Summary(group = NULL, thin = 0.1, width = 4)
               dev.off()
             }


             ## Consensus sequence probability
             gfile <-
               file.path(self$getOutdir(),
                         paste("plot1.conseq", self$getLrdType(), self$getMapper(), "pdf", sep = "."))
             if (self$hasMap1Alternate()) {
               pdf(file = gfile, width = 20, height = 8, onefile = TRUE)
             } else {
               pdf(file = gfile, width = 12, height = 8, onefile = TRUE)
             }

             self$plotMap1SummaryConseqProb(text_size = 1.5,
                                            point_size = 0.5,
                                            threshold = 0.75)
             dev.off()
           }

           return(invisible(self))
         })

#' @export
print.map1 <- function(x, ...) {
  msg <- sprintf("An object of class '%s'\n", class(x)[1])
  bamf <-
    paste0(sprintf("%s: %s", names(compact(x$bamfile)), sapply(compact(x$bamfile), basename)), collapse = "\n           ")
  seqp <-
    paste0(sprintf("%s: %s", names(compact(x$seqpath)), sapply(compact(x$seqpath), basename)), collapse = "\n          ")
  msg <- sprintf(
    "%s [Dir] %s\n [Reads] %s\n [Bamfile] %s\n [Pileup]\n [Conseq]\n [Seqpath] %s\n",
    msg, x$dir, basename(x$reads), bamf, seqp
  )
  cat(msg)
}

## Method: map2 ####

#' @export
map2.DR2S <- function(x,
                      opts = list(),
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
  message("\nStep 4: Mapping partitioned long reads against inferred consensus sequences ... \n")
  Sys.sleep(1)
  x$runMap2(opts = opts, pct = pct, min_base_quality = min_base_quality,
            min_mapq = min_mapq, max_depth = max_depth,
            min_nucleotide_depth = min_nucleotide_depth,
            include_deletions = include_deletions,
            include_insertions = include_insertions,
            pruning_cutoff = pruning_cutoff, force = force,
            fullname = fullname, plot = plot)
  message("  Done!\n")
  invisible(x)
}

DR2S_$set("public", "runMap2",
         function(opts = list(),
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

           ## Overide default arguments
           args <- self$getOpts("map2")
           if (!is.null(args)) {
             env  <- environment()
             list2env(args, envir = env)
           }

           ## Mapper
           map_fun <- self$getMapFun()

           reftag <- if (self$hasMap1Alternate()) {
             "merged"
           } else if (self$hasMultialign()) {
             "multialign"
           } else {
             "reference"
           }
           # group = "A"
           # group = "B"
           for (group in c("A", "B")) {
             outdir   <- self$map1[[group]]$dir
             readpath <- self$map1[[group]]$reads
             refpath  <- self$map1[[group]]$seqpath[[reftag]]

             nreads <- if (!is.null(attr(self$getHapList(group), "n"))) {
               attr(self$getHapList(group), "n")
             } else {
               length(self$getHapList(group))
             }
             optsname <- sprintf("%s [%s]", group, nreads)
             mapfmt  <- "map2 <%s> <%s> <%s> <%s>"
             maptag  <- sprintf(mapfmt, group, self$getLrdType(), self$getMapper(), optstring(opts, optsname))

             self$map2[[group]] = structure(
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
               class = c("map2", "list")
             )

             ## Run mapper
             message("  Mapping allele group ", dQuote(group), " against merged consensus ...")
             samfile <- map_fun(
               reffile  = refpath,
               readfile = readpath,
               allele   = "map2",
               readtype = self$getLrdType(),
               opts     = opts,
               refname  = reftag,
               optsname = optsname,
               force    = force,
               outdir   = outdir
             )

             ## Run bam - sort - index pipeline
             message("  Indexing ...")
             bamfile <- bam_sort_index(
               samfile,
               refpath,
               pct / 100,
               min_mapq,
               force = force,
               clean = TRUE
             )
             self$map2[[group]]$bamfile = bamfile

             ## Calculate pileup from graphmap produced SAM file
             message("  Piling up ...")
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
             self$map2[[group]]$pileup = pileup

             ## Construct consensus sequence
             message("  Constructing a consensus ...")
             conseq_name <- paste0("consensus.", sub(".sam.gz", "", basename(samfile)))
             conseq      <- conseq(x = pileup, name = conseq_name, type = "prob",
                                   exclude_gaps = TRUE, prune_matrix = TRUE,
                                   cutoff = pruning_cutoff)
             seqpath     <- file.path(outdir, paste0(conseq_name, ".fa"))
             self$map2[[group]]$seqpath = seqpath
             self$map2[[group]]$conseq = conseq
             Biostrings::writeXStringSet(
               Biostrings::DNAStringSet(gsub("[-+]", "", conseq)),
               seqpath
             )

             ## Set maptag
             self$map2[[group]]$tag[[reftag]] = maptag
           }

           if (plot) {
             message("  Plotting ...")
             ## Coverage and base frequency
             gfile <- file.path(
               self$getOutdir(),
               paste("plot2", self$getLrdType(), self$getMapper(), "pdf", sep = ".")
             )
             pdf(file = gfile, width = 16, height = 8, onefile = TRUE)
             self$plotMap2Summary(thin = 0.1, width = 10)
             dev.off()
             ## Consensus sequence probability
             gfile <- file.path(
               self$getOutdir(),
               paste("plot2.conseq", self$getLrdType(), self$getMapper(), "pdf", sep = ".")
             )
             pdf(file = gfile, width = 16, height = 8, onefile = TRUE)
             self$plotMap2SummaryConseqProb(text_size = 1.75,
                                            point_size = 0.75,
                                            threshold = 0.75)
             dev.off()
           }

           invisible(self)
         })

#' @export
print.map2 <- function(x, ...) {
  msg <- sprintf("An object of class '%s'.\n", class(x)[1])
  bamf <- basename(x$bamfile %||% "")
  msg <- sprintf(
    "%s [Dir] %s\n [Reads] %s\n [Reference] %s\n [Bamfile] %s\n",
    msg,
    x$dir,
    basename(x$reads),
    basename(x$ref),
    bamf
  )
  cat(msg)
}

## Method: map3 ####

#' @export
map3.DR2S <- function(x,
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
  message("\nStep 5: Mapping short and long reads against refined consensus sequences ... \n")
  Sys.sleep(1)
  x$runMap3(opts = opts, pct = pct, min_base_quality = min_base_quality,
            min_mapq = min_mapq, max_depth = max_depth,
            min_nucleotide_depth = min_nucleotide_depth,
            include_deletions = include_deletions,
            include_insertions = include_insertions,
            include_read_ids = include_read_ids, force = force,
            fullname = fullname, plot = plot, clip = clip)
  message("  Done!\n")
  invisible(x)
}

DR2S_$set("public", "runMap3",
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
           args <- self$getOpts("map3")
           if (!is.null(args)) {
             env  <- environment()
             list2env(args, envir = env)
           }

           ## stop if no shortreads provided
           if (is.null(self$getConfig("shortreads"))) {
             message("Cannot run map3. No shortreads provided")
             return(invisible(self))
           }

           reftag    <- "merged"
           outdir    <- dir_create_if_not_exists(file.path(self$getOutdir(), reftag))
           sreadpath <- self$getShortreads()
           areadpath <- self$map1$A$reads
           breadpath <- self$map1$B$reads
           arefpath  <- self$map2$A$seqpath
           brefpath  <- self$map2$B$seqpath

           ## Mapper
           map_fun <- self$getMapFun()

           self$map3 = structure(
             list(
               dir     = outdir,
               sreads  = sreadpath,
               areads  = areadpath,
               breads  = breadpath,
               aref    = arefpath,
               bref    = brefpath,
               bamfile = list(),
               pileup  = list(),
               tag     = list()
             ),
             class = c("map3", "list")
           )

           ## Remap long reads to the same reference sequences as short reads
           ## group = "B"
           for (group in c("A", "B")) {
             mapgroup <- paste0("LR", group)
             maptag   <- paste("map3", mapgroup, self$getLrdType(), self$getMapper(),
                               optstring(opts), sep = ".")
             refpath  <- self$map2[[group]]$seqpath
             lreadpath <- self$map1[[group]]$read

             ## Run mapper
             message("  Mapping long reads against Map2 consensus ...")
             samfile <- map_fun(
               reffile  = refpath,
               readfile = lreadpath,
               allele   = mapgroup,
               readtype = self$getLrdType(),
               opts     = opts,
               refname  = group,
               optsname = optstring(opts),
               force    = force,
               outdir   = outdir
             )

             ## Run bam - sort - index pipeline
             message("  Indexing ...")
             bamfile <- bam_sort_index(
               samfile,
               refpath,
               pct / 100,
               min_mapq,
               force = force,
               clean = TRUE
             )
             self$map3$bamfile[[mapgroup]] = bamfile

             ## Calculate pileup from graphmap produced SAM file
             message("  Piling up ...")
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
             self$map3$pileup[[mapgroup]] = pileup

             ## Set maptag
             self$map3$tag[[mapgroup]] = maptag
           }

           ## Map short reads
           ## group = "A"
           for (group in c("A", "B")) {
             mapgroup <- paste0("SR", group)
             maptag   <- paste("map3", mapgroup, self$getLrdType(), self$getMapper(),
                               optstring(opts), sep = ".")
             refpath  <- self$map2[[group]]$seqpath
             ## Run mapper
             message("  Mapping short reads against Map2 consensus ...")
             samfile <- map_fun(
               reffile  = refpath,
               readfile = sreadpath,
               # if we run shortreads against both pacbio and nanopore data
               # this hack makes sure that we can distinguish the bam files ->
               # we get pacbio.illumina.bwamem.A...bam and nanopore.illumina.bwamem.A...bam
               allele   = paste0(mapgroup, ".", self$getLrdType()),
               readtype = self$getSrdType(),
               opts     = opts,
               refname  = group,
               optsname = optstring(opts),
               force    = force,
               outdir   = outdir
             )

             if (clip) {
               message("  Trimming softclips and polymorphic ends ...")
               ## Run bam - sort - index pipeline
               bamfile <- bam_sort_index(samfile, refpath, pct / 100, min_mapq, force = force, clean = TRUE)
               ## Trim softclips
               fq <- trim_softclipped_ends(bam = Rsamtools::scanBam(bamfile)[[1]], preserve_ref_ends = TRUE)
               ## Trim polymorphic ends
               fq <- trim_polymorphic_ends(fq)
               ## Write new shortread file to disc
               fqdir  <- dir_create_if_not_exists(file.path(self$getOutdir(), "merged"))
               fqfile <- paste("sread", group, self$getMapper(), "trimmed", "fastq", "gz", sep = ".")
               fqout  <- file_delete_if_exists(file.path(fqdir, fqfile))
               ShortRead::writeFastq(fq, fqout, compress = TRUE)
               file_delete_if_exists(bamfile)
               ## Rerun mapper
               message("  Mapping trimmed short reads against Map2 consensus ... ")
               samfile <- map_fun(
                 reffile  = refpath,
                 readfile = fqout,
                 allele   = paste0(mapgroup, ".", self$getLrdType()),
                 readtype = self$getSrdType(),
                 opts     = list(A = 1, B = 4, O = 2),
                 refname  = group,
                 optsname = optstring(opts),
                 force    = force,
                 outdir   = outdir
               )
               # cleanup
               file_delete_if_exists(fqout)
             }

             ## Run bam - sort - index pipeline
             message("  Indexing ...")
             bamfile <- bam_sort_index(
               samfile,
               refpath,
               pct / 100,
               min_mapq,
               force = force,
               clean = TRUE
             )
             self$map3$bamfile[[mapgroup]] = bamfile

             ## Calculate pileup from graphmap produced SAM file
             message("  Piling up ...")
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
             self$map3$pileup[[mapgroup]] = pileup

             ## Set maptag
             self$map3$tag[[mapgroup]] = maptag
           }

           if (plot) {
             message("  Plotting ...")
             ## Coverage and base frequency
             gfile <- file.path(self$getOutdir(), paste("plot3", "LR", self$getLrdType(), self$getMapper(), "pdf", sep = "."))
             pdf(file = gfile, width = 16, height = 8, onefile = TRUE)
             self$plotMap3SummaryLR(thin = 0.25, width = 20)
             dev.off()

             gfile <- file.path(self$getOutdir(), paste("plot3", "SR", self$getLrdType(), self$getMapper(), "pdf", sep = "."))
             pdf(file = gfile, width = 16, height = 8, onefile = TRUE)
             self$plotMap3SummarySR(thin = 0.25, width = 20)
             dev.off()
           }

           return(invisible(self))
         })

#' @export
print.map3 <- function(x, ...) {
  msg  <- sprintf("An object of class '%s'.\n", class(x)[1])
  bamf <- paste0(basename(unlist(x$bamfile) %||% ""), collapse = ", ")
  seqp <- paste0(basename(unlist(x$seqpath) %||% ""), collapse = ", ")
  msg <- sprintf(
    "%s [Dir] %s\n [Longreads] %s\n [Shortreads] %s\n [References] %s, %s\n [Bamfile] %s\n [Seqpath] %s\n",
    msg, x$dir,
    paste(basename(x$areads), basename(x$breads), sep = ", "),
    paste0(basename(x$sreads), collapse = ", "),
    basename(x$aref), basename(x$bref),
    bamf, seqp
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
