

# Consensus sequences -----------------------------------------------------


#' Construct a consensus matrix from a pile-up
#'
#' @param x \code{pileup} object.
#' @param name Name for consensus sequence.
#' @param type One of "prob", "freq", or "ambig".
#' @param threshold If \code{type == "ambig"}, threshold to call an ambiguous
#' consensus call.
#' @param exclude_gaps Exclude gaps at insertion position from consensus
#' calling.
#' @param prune_matrix Exclude positions with a sudden drop in coverage from
#' consensus calling.
#' @param ... Additional arguments to \code{prune_consensus_matrix}.
#' Defaults \code{n_look_behind = 12} and \code{cutoff = 0.5}.
#'
#' @return A \code{\linkS4class{BStringSet}} object with a metadata list
#' containing the slots:
#' \describe{
#'   \item{zscore}{}
#'   \item{freq}{}
#'   \item{ambigs}{}
#'   \item{insertions}{}
#'   \item{deletions}{}
#'   \item{consmat}{}
#' }
#' @export
#'
#' @examples
#' ###
conseq <- function(x, ...) UseMethod("conseq")

#' @export
conseq.pileup <- function(x, name = "conseq", type = c("prob", "freq", "ambig"),
                          threshold = NULL, exclude_gaps = FALSE,
                          force_exclude_gaps = FALSE,
                          prune_matrix = FALSE, ...) {
  if (is.null(threshold)) {
    threshold <- x$threshold
  }
  x <- consmat(x, freq = FALSE)
  conseq(x, name = name, type = type, threshold = threshold,
         force_exclude_gaps = force_exclude_gaps,
         exclude_gaps = exclude_gaps, prune_matrix = prune_matrix, ... )
}

# ToDo: RM
# conseq.hmm <- function(x, name = NULL, outdir = getwd()) {
#   # debug
#   # x <- model
#   cons <- run_hmmer("hmmemit", x$outfile, outdir, force = TRUE)
#   seq <- Biostrings::readDNAStringSet(cons$outfile)
#   logo <- run_hmmer("hmmlogo", x$outfile, outdir, force = TRUE)
#   scores <- readr::read_table2(logo$outfile , skip = 2, col_names = FALSE)
#   scores <- scores %>%
#     dplyr::select(-X6) %>%
#     dplyr::mutate(X7 = as.numeric(gsub(")", "", X7)))
#   names(scores) <- c("Position", "A", "C", "G", "T", "height")
#   scores <- scores %>%
#     dplyr::select(A, C, G, T, height) %>%
#     dplyr::rowwise() %>%
#     dplyr::mutate(score = max(c(A, C, G, T))/ height)
#
#   metadata(seq) <- list(
#     zscore     = NULL,
#     freq       = scores$score,
#     ambigs     = NULL,
#     insertions = NULL,
#     deletions  = NULL,
#     consmat    = NULL
#   )
#   seq
# }

#' @export
conseq.matrix <- function(x, name = NULL, type = c("prob", "freq", "ambig"),
                          threshold = NULL, exclude_gaps = FALSE,
                          force_exclude_gaps = FALSE,
                          prune_matrix = FALSE, ...) {
  type <- match.arg(type, c("prob", "freq", "ambig"))
  if (type == "ambig" && is.null(threshold)) {
    stop("Must set threshold for ambiguity consensus calling!")
  }
  x <- consmat(x, freq = FALSE)
  if (prune_matrix) {
    x <- prune_consensus_matrix(cm = x, ...)
  }
  conseq <- switch(type,
    prob  = make_prob_consensus_(x, exclude_gaps, force_exclude_gaps),
    freq  = make_freq_consensus_(x, exclude_gaps, force_exclude_gaps),
    ambig = make_ambig_consensus_(x, threshold, exclude_gaps)
  )
  names(conseq) <- name
  conseq
}

#' @keywords internal
#' @export
simple_consensus <- function(x) {
  cmf <- consmat(x, freq = TRUE)
  if (any(na_ <- is.na(rowSums(cmf)))) {
    cmf[na_, ] <- 0
    cmf <- cbind(cmf, N = ifelse(na_, 1, 0))
  }
  paste0(colnames(cmf)[apply(cmf, 1, which.max)], collapse = "")
}

## strict consensus based on z-scores
make_prob_consensus_ <- function(x, exclude_gaps, as_string = FALSE) {
  x_ori <- x
  if (exclude_gaps && length(ins_ <- as.character(ins(x))) > 0) {
    x[dimnames(x)$pos %in% ins_, "-"] <- 0
  }

  # don't allow gaps at beginning and start
  maxbases <- names(unlist(apply(x, 1, function(a) list(which(a == max(a))[1]))))
  maxbase <- which(maxbases != "-")
  maxgap <- which(maxbases == "-")
  if (min(maxgap) < min(maxbase)){
    exclude_from_start <- min(which(maxbases == "-")):(min(which(maxbases != "-"))-1)
    x[exclude_from_start,"-"] <- 0
  }
  if (max(maxgap) > max(maxbase)) {
    exclude_from_end <- (max(which(maxbases != "-"))+1):max(which(maxbases == "-"))
    x[exclude_from_end,"-"] <- 0
  }

  ## remove rows which have been set to zero
  if (length(i <- which(n(x) == 0L)) > 0) {
    x <- x[-i, ]
  }
  rowsd <- mean(n(x)) / 2
  z <- sweep(sweep(x, 1, .rowMeans(x, NROW(x), NCOL(x)), `-`), 1, rowsd, `/`)
  bases <- colnames(x)[apply(z, 1, which.max)]
  if (as_string) {
    return(paste0(bases, collapse = ""))
  }
  dels <- bases == "-"
  seq  <- Biostrings::BStringSet(paste0(bases[!dels], collapse = ""))

  # fix zscore; rm del positions
  z <- z[!dels,]

  metadata(seq) <- list(
    zscore     = unname(apply(z, 1, max)),
    freq       = NULL,
    ambigs     = NULL,
    insertions = ins(x_ori),
    deletions  = unname(which(dels)),
    consmat    = x_ori
  )
  seq
}

## strict majority rule consensus
make_freq_consensus_ <- function(x, exclude_gaps, force_exclude_gaps, as_string = FALSE) {
  x_ori <- x
  if (exclude_gaps && length(ins_ <- ins(x)) > 0) {
    x[ins_, "-"] <- 0
  }
  if (force_exclude_gaps){
    x[,"-"] <- 0
  }
  cmf <- consmat(x, freq = TRUE)
  ## remove rows which have been set to zero
  if (length(i <- which(n(cmf) == 0L)) > 0) {
    cmf <- cmf[-i, ]
  }
  bases <- colnames(x)[apply(cmf, 1, which.max)]
  if (as_string) {
    return(paste0(bases, collapse = ""))
  }
  dels <- bases == "-"
  seq  <- Biostrings::BStringSet(paste0(bases[!dels], collapse = ""))
  metadata(seq) <- list(
    zscore     = NULL,
    freq       = unname(apply(cmf, 1, max)),
    ambigs     = NULL,
    insertions = ins(x_ori),
    deletions  = unname(which(dels)),
    consmat    = x_ori
  )
  seq
}

## consensus with ambiguities
## <exclude_gaps> affects behaviour at polymorphic positions:
## if <!exclude_gaps> a gap ambiguity will be called (small letter)
## if <exclude_gaps> the alternate base(s) will be called irrespective of the frequency of the gap
make_ambig_consensus_ <- function(x, threshold, exclude_gaps = FALSE, as_string = FALSE) {
  ## Filter all bases with a frequency > threshold
  cmf <- consmat(x, freq = TRUE)
  ## remove rows which have been set to zero
  if (length(i <- which(n(cmf) == 0L)) > 0) {
    cmf <- cmf[-i, ]
  }
  baselist <- apply(cmf, 1, function(m) {
    rs <- m[i <- m > threshold]
    names(rs) <- names(m)[i]
    list(rs)
  })
  #   baselist <- foreach(m = iter(consmat(x, freq = TRUE), by = "row")) %do% {
  #     i <- m > threshold
  #     rs <- m[i]
  #     names(rs) <- nm[i]
  #     rs
  #   }
  ## for each position return a BASE and a FREQ
  ## b <- baselist[[4020]]
  ## b <- baselist[ins(cmf)]
  ## n(cmf)[ins(cmf)]
  ## ins(cmf)
  ## baselist2 <- baselist[sapply(baselist, length) > 1]
  ## b <- baselist2[[1]]
  ## b <- baselist[[68]]
  s <- lapply(baselist, function(b) {
    b <- unlist(b)
    if (length(b) == 1) {
      list(base = names(b), freq = unname(b))
    }
    else if (length(b) > 1) {
      ## sort by name
      b <- b[order(names(b))]
      ## if we have a gap and exclude_gaps == TRUE
      if (any(gap <- names(b) == "-") && exclude_gaps) {
        b <- b[!gap]
      }
      NUC <- paste0(names(b), collapse = "")
      list(
        base = names(CODE_MAP())[charmatch(NUC, CODE_MAP())] %|na|% "N",
        freq = b
      )
    }
    else {
      stop("No bases?")
    }
  })
  bases <- vapply(s, `[[`, "base", FUN.VALUE = "")
  if (as_string) {
    return(paste0(bases, collapse = ""))
  }
  ambigs <- x[which(!bases %in% DNA_BASES()), ]
  attr(ambigs, "ambiguities") <- unname(bases[which(!bases %in% DNA_BASES())])
  dels <- bases == "-"
  seq  <- Biostrings::BStringSet(paste0(bases[!dels], collapse = ""))
  metadata(seq) <- list(
    zscore     = NULL,
    freq       = unname(vapply(s, function(x) sum(x[["freq"]]), 0)),
    ambigs     = ambigs,
    insertions = ins(x),
    deletions  = unname(which(dels)),
    consmat    = x
  )
  seq
}


# Summarise and plot ------------------------------------------------------


#' Plot probability of consensus bases.
#'
#' @param seqA consensus sequence A.
#' @param seqB consensus sequence B.
#' @param text_size Text size.
#' @param point_size Point size.
#' @param labelA Label A.
#' @param labelB Label B.
#'
#' @return A \code{ggplot} object.
#' @export
#'
#' @examples
#' ##
                                    # threshold = "auto"
                                    # text_size = 3
                                    # point_size = 1
plot_conseq_probability <- function(cseqs,
                                    threshold = "auto",
                                    text_size = 3,
                                    point_size = 1) {
  labels <- sapply(cseqs, function(x) x$label)
  seqs <- sapply(cseqs, function(x) x$cseq)
  # get labels and tags
  if (all(nzchar(labels))){
    tags <- lapply(labels, function(x) gsub("[<>]", "", strsplit(x, " ", fixed = TRUE)[[1]]))
    groups <- lapply(tags, function(x) paste0(setdiff(x, Reduce(intersect, tags)), collapse = "."))
    label <- paste0(Reduce(intersect, tags), collapse = ".")
  } else{
    label <- ""
    groups <- 1:length(cseqs)
  }

  if (all(unlist(lapply(seqs, function(x) !is.null(metadata(x)$zscore))))) {
    ylabel <- "Z-score probability"
  } else {
    ylabel <- "Frequency"
  }
  df <- dplyr::bind_rows(
    foreach(hp = names(cseqs)) %do% {
      seq <- seqs[[hp]]
      dplyr::data_frame(
        group = groups[[hp]],
        pos   = seq_len(Biostrings::width(seq)),
        base  = strsplit(as.character(seq), split = "")[[1]],
        prob  = if (!is.null(metadata(seq)$zscore)) {
          as.numeric(pnorm(metadata(seq)$zscore))
        } else {
          metadata(seq)$freq
        }
      )
    }
  )

  lower <- if (threshold == "auto") {
    dplyr::summarise(dplyr::group_by(df, group), lower = quantile(prob, 0.0025))
  } else {
    dplyr::data_frame(group = unlist(groups), lower = threshold)
  }
  df2 <- dplyr::filter(dplyr::left_join(df, lower, by = "group"), prob <= lower)
  ggplot(df, aes(pos, prob)) +
    facet_grid(group ~ .) +
    geom_point(aes(colour = base), alpha = 0.5, size = point_size) +
    scale_color_manual(values = NUCCOL()) +
    geom_label(aes(x = pos, y = prob - prob*0.05, label = pos), data = df2,
               colour = "black", size = text_size) +
    geom_hline(aes(yintercept = lower), linetype = "dashed", colour = "gray20", data = lower) +
    ylim(c(0.25, 1)) +
    xlab("Position [bp]") +
    ylab(ylabel) +
    ggtitle(label) +
    theme_bw()
}

# plot_logo <- function(file) {
#   file <- file.path(outdir, "t.l")
#   logo$height[1:10]
#   ggplot
#   logo
#   l <- Biostrings::width(conseq)
#   logo[1:l,]
#   View(logo)
#
# }



