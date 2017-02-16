#' Run QuSAGE on any limma model fit
#'
#' This function is a wrapper for the QuSAGE algorithm, which tests for pathway 
#' enrichment, designed for easy integration with limma.
#'
#' @param dat Omic data matrix or matrix-like object with rows corresponding to 
#'   probes and columns to samples. Any object that can be processed by 
#'   \code{\link[limma]{getEAWP}} is acceptable.
#' @param fit An object of class \code{\link[limma]{MArrayLM}}, as created by a call to
#'   \code{\link[limma]{eBayes}}.
#' @param coef Column name or number specifying which coefficient or contrast
#'   of the linear model is of interest. 
#' @param geneSets Either a named list of pathways to be compared, or a vector of
#'   probe names representing a single gene set.
#' @param n.points The number of points at which to sample the convoluted
#'   \emph{t}-distribution. See Details.
#'
#' @details
#' \code{qlim} combines the statistical flexibility and empirical Bayes methods of 
#' \code{limma} with the specificity and sensitivity of the QuSAGE algorithm for 
#' detecting pathway enrichment. Simulations have shown that this pipeline 
#' outperforms each package's independent methods for gene set analysis in 
#' complex experimental designs, namely \code{\link[limma]{camera}} and 
#' \code{link[qusage]{qgen}}. See Watson & John, forthcoming.
#'
#' By default \code{n.points} is set to 2^14, or 16,384 points, which will give very
#' accurate \emph{p}-values in most cases. Sampling at more points will increase the
#' accuracy of the resulting \emph{p}-values, but will also linearly increase the
#' amount of time needed to calculate the result. With larger sample sizes, as few
#' as 1/4 this number of points can be used without seriously affecting the accuracy
#' of the resulting \emph{p}-values. However, when there are a small number of samples
#' (i.e., fewer than 8 samples total), the \emph{t}-distribution must be sampled over
#' a much wider range, and the number of points needed for sampling should be
#' increased accordingly. It is recommended that when running \code{qlim} with fewer
#' than 8 samples, the number of points be increased to at least 2^15 or 2^16. It may
#' also be useful to test higher values of this parameter, as it can often result in
#' a much more significant \emph{p}-value with small sample sizes.
#'
#' @return A data frame with the following columns:
#' \itemize{
#'   \item{'Pathway'} Name of the pathway.
#'   \item{'logFC'} Average log2 fold change of genes in the pathway.
#'   \item{'p.value'} The gene set \emph{p}-value, as calculated using
#'     \code{pdf.pVal}.
#'   \item{'q.value'} The false discovery rate, as estimated using Storey's
#'     \emph{q}-value method.
#' }
#'
#' @references
#' Smyth, G.K. (2004). "Linear models and empirical Bayes methods for assessing
#' differential expression in microarray experiments." \emph{Stat. Appl. Genet. 
#' Molec. Biol.}, \emph{3}(1).
#' \url{http://www.statsci.org/smyth/pubs/ebayes.pdf}
#'
#' Yaari, G. Bolen, C.R., Thakar, J. & Kleinstein, S.H. (2013). "Quantitative set
#' analysis for gene expression: a method to quantify gene set differential expression
#' including gene-gene correlations." \emph{Nucleic Acids Res.}, \emph{41}(18): e170.
#' \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3794608/}
#'
#' Turner, J.A., Bolen, C.R. & Blankenship, D.M. (2015). "Quantitative gene set analysis
#' generalized for repeated measures, confounder adjustment, and continuous covariates."
#' \emph{BMC Bioinformatics}, \emph{16}(1): 272.
#' \url{http://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-015-0707-9}
#'
#' Storey, J. & Tibshirani, R. (2003). "Statistical significance for genomewide
#' studies." \emph{Proc. Natl. Acad. Sci.}, \emph{100}(16): 9440–9445.
#' \url{http://www.pnas.org/content/100/16/9440.full}
#'
#' @examples
#' # Fit limma model
#' library(limma)
#' eset <- matrix(rnorm(5000 * 10), nrow = 5000, ncol = 10)
#' dimnames(eset) <- list(seq_len(nrow(eset)), 
#'                        paste0('S', seq_len(ncol(eset))))
#' clin <- data.frame(treat = rep(c("ctl", "trt"), each = 5),
#'                    batch = rep(c("b1", "b2"), times = 5),
#'                       x1 = rnorm(10),
#'                       x2 = runif(10))
#' des <- model.matrix(~ treat + batch + x1 + x2, data = clin)
#' fit <- eBayes(lmFit(eset, des))
#' 
#' # Create list of differentially expressed pathways
#' geneSets = list()
#' for (i in 0:10) {
#'   genes <- ((30 * i) + 1):(30 * (i + 1))
#'   eset[genes, clin$treat == "trt"] <- eset[genes, clin$treat == "trt"] + rnorm(1)
#'   geneSets[[paste("Set", i)]] <- genes
#' }
#' 
#' # Run qlim
#' qlim(eset, fit, coef = 2, geneSets)
#'
#' @export
#' @import limma
#' @import qusage
#' @import qvalue
#' @import dplyr
#'

qlim <- function(dat,
                 fit,
                 coef,
                 geneSets,
                 n.points = 2^14) {

  # Preliminaries
  if (nrow(dat) < 3) {
    stop('dat must have at least three probes.')
  }
  if (nrow(dat) != nrow(fit) || ncol(dat) != nrow(fit$design)) {
    stop('dat is not conformal with fit.')
  }
  if (!identical(rownames(dat), rownames(fit))) {
    stop('dat and fit must have identical rownames.')
  }
  if (!is(fit, 'MArrayLM')) {
    stop('fit must be an MArrayLM object.')
  }
  if (is.null(fit$t) && is.null(fit$F)) {
    stop('fit must undergo variance moderation through eBayes or treat before ',
         'running qlim.')
  }
  if (is.null(fit$coefficients)) {
    stop('Coefficients not found in fit object.')
  }
  if (is.character(coef) && !coef %in% fit$coefficients) {
    stop(paste0('"', coef, '" not found in fit$coefficients. If passing a string to ',
                'coef, make sure it matches one of colnames(fit$coefficients).'))
  } 
  
  # Prep data
  dat <- getEAWP(dat)
  dat <- dat$exprs
  se <- sqrt(fit$s2.post) * fit$stdev.unscaled[, coef]
  sd.a <- se / (fit$sigma * fit$stdev.unscaled[, coef])
  sd.a[is.infinite(sd.a)] <- 1
  resid_mat <- residuals.MArrayLM(fit, dat)
  overlap <- sapply(geneSets, function(p) sum(p %in% rownames(dat)))
  geneSets <- geneSets[overlap > 1]

  # Run QuSAGE functions
  res <- newQSarray(mean = fit$coefficients[, coef],  # Create QSarray obj
                      SD = se,
                sd.alpha = sd.a,
                     dof = fit$df.total,
                  labels = rep('resid', ncol(dat)))
  res <- aggregateGeneSet(res, geneSets, n.points)    # PDF per gene set
  res <- calcVIF(resid_mat, res, useCAMERA = FALSE)   # VIF on resid_mat

  # Export
  out <- qsTable(res, number = Inf, sort.by = 'p') %>%
    rename(p.value = p.Value,
           Pathway = pathway.name,
             logFC = log.fold.change) %>%
    mutate(q.value = qvalue(p.value)$qvalues) %>%
    select(Pathway:p.value, q.value)
  return(out)

}


