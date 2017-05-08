#' Simulate Data Using a Limma Model Fit
#'
#' This function evaluates the impact of different expression filters on the 
#' residuals of a fitted DESeqDataSet.
#'
#' @param dat Omic data matrix or matrix-like object with rows corresponding to probes
#'   and columns to samples.
#' @param fit An object of class \code{limma::\link[limma]{MArrayLM}}, as created by 
#'   a call to \code{\link[limma]{eBayes}}.
#' @param nsim Number of simulations to run. Each simulation generates a matrix of 
#'   dimensionality equal to \code{dat}. For \code{nsim > 1}, matrices are \code{
#'   cbind}ed together as if adding new samples to the dataset.   
#'
#' @details
#' 
#'
#' @examples
#' library(limma)
#' mat <- matrix(rnorm(1000 * 6), nrow = 1000, ncol = 6)
#' grp <- rep(c("A", "B"), 3)
#' des <- model.matrix(~ grp)
#' fit <- eBayes(lmFit(mat, des))
#' sims <- sim_lim(mat, fit, nsim = 10)
#'
#' @export
#'

sim_lim <- function(dat, 
                    fit, 
                    nsim = 1) {
  
  p <- nrow(dat)
  n <- ncol(dat)
  signal_mat <- matrix(rep(fitted(fit), nsim), nrow = p, ncol = n * nsim)
  noise_mat <- matrix(rnorm(n = p * n * nsim, mean = 0L, sd = sqrt(fit$s2.post)), 
                      nrow = p, ncol = n * nsim)
  out <- signal_mat + noise_mat
  dimnames(out) <- list(rownames(dat), 
                        paste0(colnames(dat), '_sim', rep(seq_len(nsim), each = n)))
  return(out) 
  
}