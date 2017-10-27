### GENERALIZED FEATURE WEIGHT ESTIMATION ###

imp <- function(mod, 
                dat, 
                nsim = 100L, 
              robust = TRUE, 
            parallel = TRUE) {

  # Load libraries
  require(matrixStats)

  # Prepare output
  n <- nrow(dat)
  p <- ncol(dat)
  out <- matrix(nrow = n, ncol = p)

  # Windows
  if (robust) {
    intervals <- colMads(as.matrix(dat))
  } else {
    intervals <- colSds(as.matrix(dat))
  }

  # Loop
  y_hat <- predict(mod, as.data.frame(dat))
  for (i in seq_len(n)) {    # Initialize
    x_tilde <- matrix(rep(dat[i, ], nsim), ncol = p, byrow = TRUE,
                      dimnames = list(NULL, colnames(dat)))
    for (j in seq_len(p)) {  # Simulate 
      x_tilde[, j] <- runif(nsim, 
                            min = dat[i, j] - 2L * intervals[j], 
                            max = dat[i, j] + 2L * intervals[j])
      y_tilde <- try(predict(mod, x_tilde), silent = TRUE)
      if (is(y_tilde, 'try-error')) {
        y_tilde <- predict(mod, as.data.frame(x_tilde))
      }
      if (robust) {          # Estimate partial derivative
        out[i, j] <- median((y_tilde - y_hat[i]) / (x_tilde[, j] - dat[i, j]))
      } else {
        out[i, j] <- mean((y_tilde - y_hat[i]) / (x_tilde[, j] - dat[i, j]))
      }
      
    }
  }

  # Export
  out <- colMeans(out)
  names(out) <- colnames(dat)
  return(out)

}


