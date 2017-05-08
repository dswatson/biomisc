filt_dist <- function(dat, 
                      top, 
                      dist = 'euclidean') {
  
  # Preliminaries
  
  # 
  
  dm <- matrix(0L, nrow = ncol(dat), ncol = ncol(dat))
  for (i in 2L:ncol(dat)) {
    for (j in 1L:(i - 1L)) {
      if (dist == 'euclidean') {
        top_idx <- nrow(dat) - top + 1L
        dm[i, j] <- sqrt(sum(sort.int((dat[, i] - dat[, j])^2L,
                                      partial = top_idx)[top_idx:nrow(dat)]))
      } else {
        tops <- order((dat[, i] - dat[, j])^2, decreasing = TRUE)[seq_len(top)]
        if (dist == 'pearson') {
          dm[i, j] <- 1 - cor(dat[tops, i], dat[tops, j])
        } else if (dist == 'MI') {
          dm[i, j] <- max(as.matrix(MIdist(t(dat[tops, c(i, j)]))))
        } else if (dist == 'KLD') {
          dm[i, j] <- max(as.matrix(KLdist.matrix(t(dat[tops, c(i, j)]))))
        }
      }
    }
  }
  dm <- pmax(dm, t(dm))
  
  # Export
  return(dm)
  
}
    
    
    
    
    
    
    
    