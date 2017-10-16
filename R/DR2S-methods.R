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
                      include_deletions = FALSE,
                      include_insertions = FALSE,
                      force = FALSE,
                      fullname = TRUE,
                      # for now set plot to FALSE; cant plot with the clustering
                      plot = TRUE) {
  message("\nStep 0: Initial mapping of long reads ... \n")
  Sys.sleep(1)
  x$runMapInit(opts = opts, optsname = optsname, partSR = partSR, pct = pct, threshold = threshold,
            min_base_quality = min_base_quality, min_mapq = min_mapq,
            max_depth = max_depth, min_nucleotide_depth = min_nucleotide_depth,
            include_deletions = include_deletions, include_insertions = include_insertions,
            force = force, fullname = fullname, plot = plot)
  message("  Done!\n")
  invisible(x)
}
# self <- hla
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
                  include_deletions = FALSE,
                  include_insertions = FALSE,
                  force = FALSE,
                  fullname = TRUE,
                  plot = TRUE) {
           stopifnot(pct > 0 && pct <= 100)
           recode_header <- FALSE

           # debug
           # opts = list()
           # partSR = TRUE
           # optsname = ""
           # pct = 100
           # threshold = 0.20
           # min_base_quality = 3
           # min_mapq = 0
           # max_depth = 1e4
           # min_nucleotide_depth = 3,
           # include_insertions = FALSE
           # force = FALSE
           # fullname = TRUE
           # plot = TRUE

           ## Overide default arguments
           args <- self$getOpts("mapInit")
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
           mapfmt  <- "mapInit <%s> <%s> <%s> <%s>"
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
             # clean = TRUE
             # debug
             clean = FALSE
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

           self$mapInit = structure(
             list(
               reads   = self$getLongreads(),
               bamfile = bamfile,
               pileup  = pileup,
               tag     = maptag
             ),
             class  = c("mapInit", "list")
           )
           if (partSR){
             mapfmt  <- "mapInit <%s> <%s> <%s> <%s>"
             maptag  <- sprintf(mapfmt, self$getReference(), self$getSrdType(), self$getMapper(), optstring(opts, optsname))
             if (threshold != self$getThreshold()) {
               self$setThreshold(threshold)
             }
             ## Fetch mapper
             map_fun <- self$getMapFun()
             ## Run mapper
             message("  Mapping short reads against provided reference...\n")
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
             message("  Indexing ...")
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
             message("  Piling up ...")
             pileup <- Pileup(
               bamfile,
               self$getThreshold(),
               max_depth,
               min_base_quality = min_base_quality,
               min_mapq = min_mapq,
               min_nucleotide_depth = min_nucleotide_depth,
               include_deletions = include_deletions,
               include_insertions = FALSE#include_insertions
             )

             self$mapInit$SR = structure(
               list(
                 reads   = self$getShortreads(),
                 bamfile = bamfile,
                 pileup  = pileup,
                 tag     = maptag
               ),
               class  = c("mapInit", "list")
             )
           }
           if (plot) {
             message("  Plotting ...")
             ## Coverage and frequency of minor alleles
             gfile <- file.path(self$getOutdir(), paste0("plot.MapInit.LR.", sub("bam$", "pdf", usc(basename(self$mapInit$bamfile)))))
             pdf(file = gfile, width = 12, height = 10, onefile = TRUE)
             self$plotmapInitSummary(
               thin = 0.25,
               width = 2,
               label = self$getMapTag("init", "LR"),
               readtype = "LR"
             )
             dev.off()

             if (partSR){
               gfile <- file.path(self$getOutdir(), paste0("plot.MapInit.SR.", sub("bam$", "pdf", usc(basename(self$mapInit$SR$bamfile)))))
               pdf(file = gfile, width = 12, height = 10, onefile = TRUE)
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
                                      # shuffle = TRUE,
                                      skip_gap_freq = 2/3,
                                      plot = TRUE,
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
#self <- hla
DR2S_$set("public", "runHaplotypePartitioning", function(max_depth = 1e4,
                  skip_gap_freq = 2/3,
                  plot = TRUE) {

                  # debug
                  max_depth = 1e4
                  skip_gap_freq = 2/3
                  plot = TRUE
           stopifnot(self$hasPileup())

           ## Overide default arguments
           args <- self$getOpts("partition")
           if (!is.null(args)) {
             env  <- environment()
             list2env(args, envir = env)
           }

           tag  <- self$getMapTag("init", "LR")
           th   <- self$getThreshold()

          if (!is.null(self$mapInit$SR)){
            useSR <- TRUE
          } else {
            useSR <- FALSE
          }
           ppos <- self$polymorphicPositions(useSR = useSR)

           if (NROW(ppos) < 2){
             stop("Less than two polymorphic positions")
           }
           message("  Generating SNP matrix ...")
           mat <- if (tryCatch(
             !is(self$partition, "PartList"),
             error = function(e)
               TRUE
           ) ||
           !(all(ppos$position %in% colnames(self$partition$mat)) &&
             all(colnames(self$partition$mat) %in% ppos$position))) {
             SNPmatrix(bamfile = self$mapInit$bamfile, refseq = self$getRefSeq(), polymorphic_positions = ppos)
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

           message(sprintf("  Partitioning %s reads over %s SNPs ...", NROW(mat), NCOL(mat)))
           prt <- partition_reads(x = mat, skip_gap_freq = skip_gap_freq, deepSplit = 1, outdir = self$getOutdir())

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

# self <- dedk.part
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
             pdf(file = outf, width = 10, height = 12, onefile = TRUE)
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
  message("\nStep2: Extracting haplotyped reads ...")
  Sys.sleep(1)
  x$extractFastq(nreads = nreads, replace = replace, nalign = nalign)
  message("  Done!\n")
  invisible(x)
}
# debug
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
              if (FALSE) {
                message("   Constructing initial ", hptype, " consensus from ", nalign, " reads ...")
                # msa <- multialign(self, hptype, nalign)
                # DECIPHER::BrowseSeqs(msa)
                cons <- multialign_consensus(multialign(self, hptype, nalign))
                consfile <- paste("consensus", self$getLrdType(), "multialign", paste0(hptype, nalign), "fa", sep = ".")
                names(cons) <- consfile
                consout <- file_delete_if_exists(file.path(dir, consfile))
                Biostrings::writeXStringSet(cons, consout)
                self$mapIter[["0"]][[hptype]]$ref <- "multialign"
                self$mapIter[["0"]][[hptype]]$conseq <- cons
                self$mapIter[["0"]][[hptype]]$seqpath <- consout
                self$mapIter[["0"]][[hptype]]$tag <- "multialign"
              } else {
                self$mapIter[["0"]][[hptype]]$ref <- "reference"
                self$mapIter[["0"]][[hptype]]$conseq <- self$getRefSeq()
                self$mapIter[["0"]][[hptype]]$seqpath <- self$getRefPath()
                self$mapIter[["0"]][[hptype]]$tag <- "reference"
              }
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
  message("\nStep 4: Mapping partitioned long reads against inferred consensus sequences ... \n")
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
#self <- dl3
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

           # iterations = 2
           for (iteration in 1:iterations) {
             message("Iteration ", iteration, " of ", iterations)

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
               message("  Mapping Haplotype ", dQuote(hptype), " against merged consensus ...")
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
               message("  Indexing ...")
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
               self$mapIter[[iterationC]][[hptype]]$pileup = pileup

               # ## Construct consensus sequence
               message("  Constructing a consensus ...")
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
               message("  Plotting ...")
               ## Coverage and base frequency
               gfile <- file.path(
                 self$getOutdir(),
                 paste("plot.MapIter", iteration, self$getLrdType(), self$getMapper(), "pdf", sep = ".")
               )
               pdf(file = gfile, width = 8 * length(hptypes), height = 8, onefile = TRUE)
               self$plotmapIterSummary(thin = 0.1, width = 4, iteration = iteration, drop.indels = TRUE)
               ## map1: 0.1; 4
               dev.off()
               ## Consensus sequence probability
               gfile <- file.path(
                 self$getOutdir(),
                 paste("plot.MapIter.conseq", iteration, self$getLrdType(), self$getMapper(), "pdf", sep = ".")
               )
               pdf(file = gfile, width = 16, height = 4*length(hptypes), onefile = TRUE)
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
                      force = FALSE,
                      plot = TRUE) {
  message("\nStep 5: Partition shortreads ... \n")
  Sys.sleep(1)
  x$runPartitionShortReads(opts = opts,
            force = force, plot = plot)
  message("  Done!\n")
  invisible(x)
}
#self <- dl3.map
DR2S_$set("public", "runPartitionShortReads",
         function(opts = list(),
                  force = FALSE,
                  plot = TRUE) {

           ## Overide default arguments
           args <- self$getOpts("partitionSR")
           if (!is.null(args)) {
             env  <- environment()
             list2env(args, envir = env)
           }

           ## stop if no shortreads provided
           if (is.null(self$getConfig("shortreads"))) {
             message("Cannot partition shortreads. No shortreads provided")
             return(invisible(self))
           }

           bamfile <- Rsamtools::BamFile(self$mapInit$SR$bamfile)
           ref <- self$getRefSeq()
           refname <- names(ref)
           hptypes <- self$getHapTypes()
           #' dont use the first matrix, but use only the high scoring
           #' reads form split_reads_by_haplotype
           # mats <- mats(self$partition$prt)
           prt_mat <- self$partition$mat
           seqs <- lapply(self$partition$hpl, function(x) get_seqs_from_mat(prt_mat[x,]))
           names(seqs) <- hptypes
           mats <-  lapply(seqs,  function(x) Biostrings::consensusMatrix(x,  as.prob = TRUE)[VALID_DNA(),])
           mats <- foreach(m = mats) %do% {
             m[m == 0] <- 0.001
             colnames(m) <- colnames(prt_mat)
             m
           }
           names(mats) <- names(seqs)
           ppos <- colnames(mats[[1]])

           # Run partitioning
           srpartition <- get_SR_partition_scores(ppos, refname, bamfile, mats, cores = "auto")

           ##' Assign read to haplotype with highest probability,
           ##' i.e. prod over probabilities of each haplotype and choose max
           message("get highest scoring haplotype for each read")
           # debug
           #srpartition <- self$srpartition
           # srpartition <- srpartition$A$srpartition
           # srpartition
           srpartition$haplotypes <- score_highest_SR(srpartition$srpartition)

           # Write fastqs
           #hptype = "B"
           foreach(hptype = hptypes ) %do% {
             srfilenames <- c()
             message("Write shortread fastq for haplotype ", hptype, " ...")
             fqs <- self$getShortreads()
             dontUseReads <- srpartition$haplotypes$read[!srpartition$haplotypes$read %in% dplyr::filter(srpartition$haplotypes, haplotype == hptype)$read]
             #fq <- fqs[1]
             # write fastq's
             foreach(fq = fqs) %do% {
               srFastqHap = file.path(self$getOutdir(), hptype, paste0(c(strsplit(basename(fq), "\\.")[[1]][1], hptype, ".fastq.gz"), collapse = "."))
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
  message("\nStep 5: Mapping short and long reads against refined consensus sequences ... \n")
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
# self <- dedk.partSR
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
             message("Cannot run mapFinal. No shortreads provided")
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
             refpath  <- refpaths[[hptype]]

             mapgroupLR <- paste0("LR", hptype)
             maptagLR   <- paste("mapFinal", mapgroupLR, self$getLrdType(), self$getMapper(),
                               optstring(opts), sep = ".")
             readpathLR <- readpathsLR[[hptype]]
             ## Run mapper
             message("  Mapping long reads against latest consensus ...")
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
             message("  Indexing ...")
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
             self$mapFinal$pileup[[mapgroupLR]] = pileup

             ## Set maptag
             self$mapFinal$tag[[mapgroupLR]] = maptagLR


           ## Map short reads
             mapgroupSR <- paste0("SR", hptype)
             maptagSR   <- paste("mapFinal", mapgroupSR, self$getLrdType(), self$getMapper(),
                               optstring(opts), sep = ".")

             readfiles <- readpathsSR[[hptype]]
             # debug:
             #readfiles <- self$getShortreads()

             ## Run mapper
             message("  Mapping short reads against latest consensus ...")

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
               message("  Trimming softclips and polymorphic ends ...")
               ## Run bam - sort - index pipeline
               bamfile <- bam_sort_index(samfile, refpath, pct / 100, min_mapq, force = force, clean = TRUE)
               ## Trim softclips
               fq <- trim_softclipped_ends(bam = Rsamtools::scanBam(bamfile)[[1]], preserve_ref_ends = TRUE)
               ## Trim polymorphic ends
               fq <- trim_polymorphic_ends(fq)
               ## Write new shortread file to disc
               fqdir  <- dir_create_if_not_exists(file.path(self$getOutdir(), "merged"))
               fqfile <- paste("sread", hptype, self$getMapper(), "trimmed", "fastq", "gz", sep = ".")
               fqout  <- file_delete_if_exists(file.path(fqdir, fqfile))
               ShortRead::writeFastq(fq, fqout, compress = TRUE)
               file_delete_if_exists(bamfile)
               ## Rerun mapper
               message("  Mapping trimmed short reads against latest consensus ... ")
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
             message("  Indexing ...")
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
             self$mapFinal$pileup[[mapgroupSR]] = pileup

             ## Set maptag
             self$mapFinal$tag[[mapgroupSR]] = maptagSR

             # calc new consensus
             cseq <- conseq(pileup$consmat, paste0("mapFinal", hptype), "ambig", exclude_gaps = TRUE, threshold = self$getThreshold())
             self$mapFinal$seq[[hptype]] <- cseq
           }


           if (plot) {
             message("  Plotting ...")
             ## Coverage and base frequency

             for (readtype in c("LR", "SR")){
               gfile <- file.path(self$getOutdir(), paste("plotFinal", readtype, self$getLrdType(), self$getMapper(), "pdf", sep = "."))
               pdf(file = gfile, width = 8 * length(hptypes), height = 8, onefile = TRUE)
               self$plotmapFinalSummary(iteration = "final", readtype = readtype, thin = 0.25, width = 20)
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
