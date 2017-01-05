#' Multiple plot function
#'
#' ggplot objects can be passed in \code{...}, or to \code{plotlist} (as a list of ggplot objects)
#'
#' If the layout is something like \code{matrix(c(1, 2, 3, 3), nrow = 2, byrow = TRUE)},
#' then plot 1 will go in the upper left, 2 will go in the upper right, and
#' 3 will go all the way across the bottom.
#'
#' @param ... ggplot objects.
#' @param plotlist A list of ggplot objects.
#' @param cols Number of columns in layout.
#' @param layout  A matrix specifying the layout. If present, 'cols' is ignored.
#' @author Winston Chang
#' @references \url{http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/}
#' @keywords internal
multiplot <- function(..., plotlist = NULL, cols = 1, layout = NULL) {
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  num_plots <- length(plots)
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(num_plots / cols)),
                     ncol = cols, nrow = ceiling(num_plots / cols))
  }
  if (num_plots == 1) {
    suppressWarnings(print(plots[[1]]))
  } else {
    # Set up the page
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow(layout), ncol(layout))))
    # Make each plot, in the correct location
    for (i in seq_len(num_plots)) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      suppressWarnings(
        print(plots[[i]], vp = grid::viewport(layout.pos.row = matchidx$row,
                                              layout.pos.col = matchidx$col))
      )
    }
  }
}
