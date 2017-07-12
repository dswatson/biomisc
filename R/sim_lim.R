#' Simulate Data Using a Limma Model Fit
#'
#' This function generates simulated data with the same signal structure as a
#' provided dataset.
#'
#' @param dat Omic data matrix or matrix-like object with rows corresponding to
#'   probes and columns to samples.
#' @param fit An object of class \code{\link[limma]{MArrayLM}}, as created by a
#'   call to \code{\link[limma]{eBayes}}.
#' @param nsim Number of simulations to run. Each simulation generates a matrix
#'   of dimensionality equal to \code{dat}. For \code{nsim > 1}, matrices are
#'   \code{cbind}ed together as if adding new samples to the dataset.
#' @param rmse Method for estimating the root mean square error of probewise
#'   regressions. Must be one of either \code{lmFit}, which uses observed model
#'   residuals, or \code{eBayes}, which uses shrunken estimates. 
#' @param fix_coefs Logical indicating whether model coefficients should be 
#'   taken as fixed. If \code{FALSE}, fold changes are drawn from random normal 
#'   distributions with mean equal to the observed coefficients and standard
#'   deviations equal to the standard error of the estimates.
#' @param cv Perform cross validation to estimate model error and fold changes?
#'   Must be either \code{NULL}, in which case no cross validation is performed, 
#'   or an integer less than or equal to the study sample size indicating the 
#'   number of folds.
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
#' @importFrom matrixStats 
#'

sim_lim <- function(dat,
                    fit,
                    nsim = 1,
                    rmse = 'lmFit',
                fix_coef = TRUE,
                      cv = NULL) {

  g <- nrow(dat)
  n <- ncol(dat)
  if (fix_coef) {
    signal_mat <- matrix(rep(fitted(fit), nsim), nrow = g, ncol = n * nsim)
  } else {
    
  }
  if (rmse == 'eBayes') {
    SD <- sqrt(fit$s2.post)
  } else if (rmse == 'lmFit') {
    residual_mat <- residuals(fit, dat)
    SD <- rowSds(residual_mat)
  } else {
    
  }
  noise_mat <- matrix(rnorm(g * n * nsim, mean = 0L, sd = SD),
                      nrow = g, ncol = n * nsim)
  out <- signal_mat + noise_mat
  dimnames(out) <- list(rownames(dat),
                        paste0(colnames(dat), '_sim', rep(seq_len(nsim), each = n)))
  return(out)

}
