##
## Phasing of short reads during the assembly step
##
phaseplot <- function(mat) {
  mat1 <- phasebreaks(mat)
  ggplot(mat, aes(xmin = posx - 15, xmax = posx + 15, ymin = posy_lower, ymax = posy_upper)) +
    geom_rect(aes(fill = linkage), alpha = 0.5) +
    geom_rect(data = mat1, fill = "darkred", alpha = 0.5) +
    geom_text(aes(x = posx, y = posx, label = posx), dplyr::filter(mat1, breakpos),
              angle = -60, hjust = -0.25, vjust = 0.5, size = 3) +
    scale_fill_gradient2() +
    geom_abline(linetype = "dashed", colour = "gray60") +
    xlab("Position (bp)") +
    ylab("Position (bp)")
}

phasebreaks <- function(mat) {
  if (is.null(mat)) {
    return(NULL)
  }
  tmp <- mat %>%
    dplyr::filter(phase == -1) %>%
    dplyr::mutate(dist = y - x - 1) %>%
    dplyr::filter(dist == 0)
  mat1 <- dplyr::filter(mat, x %in% tmp$x | y %in% tmp$y, phase == -1)
  z <- intersect(tmp$x, tmp$y)
  breakpoints <- sort(c(z, dplyr::filter(tmp, !x %in% z & !y %in% z)$x))
  breakpos <- dplyr::filter(tmp, x %in% breakpoints) %>%
    dplyr::select(x, y) %>%
    dplyr::mutate(breakpos = TRUE)
  dplyr::left_join(mat1, breakpos, by = c("x", "y")) %>%
    dplyr::mutate(breakpos = ifelse(is.na(breakpos), FALSE, breakpos))
}

phasematrix <- function(var) {
  if ((n <- length(var)) == 0) {
    df <- data.frame(matrix(ncol = 8, nrow = 0))
    names(df) <- c("x", "y", "posx", "posy", "posy_lower", "posy_upper", "linkage", "phase")
    return(tibble::as_tibble(df))
  }
  mat <- matrix(rep(0, 5*((n*n) - n)/2), ncol = 5)
  for (i in seq_len(n - 1)) {
    for (j in seq.int(i + 1, n)) {
      mat[(i - 1)*n - ((i*i) - i)/2 + j - i, ] <- c(i, j, phasereads_(var[[i]], var[[j]]))
    }
  }
  # Set phasebreaks with linkage < 0.001 to 0. Most cases are caused by the span between positions near readlength.
  mat[,5][abs(mat[,5]) < 0.001] <- 0
  
  df <- data.frame(mat)
  colnames(df) <- c("x", "y", "posx", "posy", "linkage")
  df <- dplyr::tbl_df(df) %>%
    dplyr::filter(linkage != 0) %>%
    dplyr::mutate(phase = ifelse(linkage > 0, 1, -1)) %>%
    dplyr::mutate(upper = posy - posx) %>%
    dplyr::group_by(x) %>%
    dplyr::mutate(lower = lag(upper, 1, default = 0)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      posy_lower = posx + lower + 1,
      posy_upper = posx + upper
    ) %>%
    dplyr::select(x, y, posx, posy, posy_lower, posy_upper, linkage, phase)
  df
}

phasereads_ <- function(p1, p2) {
  p1pos <- attr(p1, "position")
  p2pos <- attr(p2, "position")

  r1r <- unique(attr(p1, "reads_ref"))
  r1a <- unique(attr(p1, "reads_alt"))
  r2r <- unique(attr(p2, "reads_ref"))
  r2a <- unique(attr(p2, "reads_alt"))

  l1r <- length(r1r)
  l1a <- length(r1a)
  l2r <- length(r2r)
  l2a <- length(r2a)

  prr <- length(intersect(r1r, r2r))/min(l1r, l2r)
  pra <- length(intersect(r1r, r2a))/min(l1r, l2a)
  par <- length(intersect(r1a, r2r))/min(l1a, l2r)
  paa <- length(intersect(r1a, r2a))/min(l1a, l2a)

  c(
    pos1 = p1pos,
    pos2 = p2pos,
    linkage = prr*paa - pra*par
  )
}
