### QQ and MD plots 

plot_quant <- function(x, 
                       y,
                       method = 'QQ',
                       pts = 1e3L,
                       main = NULL,
                       xlab = NULL,
                       ylab = NULL) {
  
  # Preliminaries
  x <- x[is.finite(x)]
  if (length(x) < 1L) {
    stop('x must have at least one finite, non-missing value.')
  }
  y <- y[is.finite(y)]
  if (length(y) < 1L) {
    stop('y must have at least one finite, non-missing value.')
  }
  if (!method %in% c('QQ', 'MD')) {
    stop('method must be either "QQ" or "MD".')
  }
  if (is.null(main)) {
    if (method == 'QQ') main <- 'QQ Plot'
    else main <- 'MD Plot'
  }
  if (is.null(xlab)) xlab <- 'X'
  if (is.null(ylab)) ylab <- 'Y'
  
  # Tidy dadta
  x <- quantile(x, probs = seq(0L, 1L, length.out = pts))
  y <- quantile(y, probs = seq(0L, 1L, length.out = pts))
  if (method == 'QQ') df <- data_frame(X = x, Y = y)
  else df <- data_frame(X = (x + y) / 2L, Y = x - y)
  
  # Build plot
  p <- ggplot(df, aes(X, Y)) + 
    geom_point(size = 0.5) + 
    labs(title = main, x = xlab, y = ylab) + 
    theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5))
  
  if (method == 'QQ') {
    p <- p + geom_abline(intercept = 0L, slope = 1L, color = 'red', size = 0.2)
  } else {
    p <- p + geom_hline(yintercept = 0L, color= 'red', size = 0.2)
  }
  
  print(p)
}


