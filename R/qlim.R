#' Run QuSAGE on any limma model fit
#'
#' This function is a wrapper for the QuSAGE algorithm, which tests for enrichment
#' across pathways or gene sets.
#'
#' @param dat Omic data matrix with rows corresponding to probes and columns
#'   to samples.
#' @param fit An object of class \code{MArrayLM}, as created by a call to
#'   \code{limma::lmFit} or \code{limma::eBayes}.
#' @param coef Column name specifying which coefficient or contrast
#'   of the linear model is of interest.
#' @param geneSets Either a named list of pathways to be compared, or a vector of
#'   probe names representing a single gene set.
#' @param n.points The number of points at which to sample the convoluted
#'   \emph{t}-distribution. See Details.
#'
#' @details
#' \code{qlim} extends the methodology of \code{qusage::qgen} to studies of
#' arbitrary design complexity by operating directly on the output of a \code{limma}
#' model fit. Taking advantage of that package's extremely flexible statistical
#' framework, which allows for continuous and/or categorical covariates, probe and/or
#' sample weights, intra-class correlation estimates, and empirical Bayes shrinkage,
#' \code{qlim} can take the results of any omic experiment analyzable in \code{limma}
#' and produce a table of enriched pathways for any given coefficient. See \pkg{limma}
#' and \pkg{qusage} for more details.
#'
#' By default \code{n.points} is set to 2^14, or 16,384 points, which will give very
#' accurate \emph{p}-values in most cases. Sampling at more points will increase the
#' accuracy of the resulting \emph{p}-values, but will also linearly increase the
#' amount of time needed to calculate the result. With larger sample sizes, as few
#' as 1/4 this number of points can be used without seriously affecting the accuracy
#' of the resulting \emph{p}-values, however, when there are a small number of samples
#' (i.e., fewer than 8 samples total), the \emph{t}-distribution must be sampled over
#' a much wider range, and the number of points needed for sampling should be
#' increased accordingly. It is recommended that when running \code{qlim} with fewer
#' than 8 samples the number of points be increased to at least 2^15 or 2^16. It may
#' also be useful to test higher values of this parameter, as it can often result in
#' a much more significant \emph{p}-value with small sample sizes.
#'
#' @return A data frame with the following columns:
#' \describe{
#'   \item{One}{Pathway} Name of the pathway.
#'   \item{Two}{logFC} Average log2 fold change of genes in the pathway.
#'   \item{Three}{p.value} The gene set \emph{p}-value, as calculated using
#'     \code{pdf.pVal}.
#'   \item{Four}{q.value} The false discovery rate, as estimated using Storey's
#'     \emph{q}-value method.
#' }
#'
#' @references
#' Smyth, G.K. (2004). "Linear models and empirical Bayes methods for assessing
#' differential expression in microarray experiments." \emph{Statistical Applications
#' in Genetics and Molecular Biology}, \emph{3}(1).
#' \url{http://www.statsci.org/smyth/pubs/ebayes.pdf}
#'
#' Yaari, G. Bolen, C.R., Thakar, J. & Kleinstein, S.H. (2013). "Quantitative set
#' analysis for gene expression: a method to quantify gene set differential expression
#' including gene-gene correlations." \emph{Nucleic Acids Research}, \emph{41}(18): e170.
#' \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3794608/}
#'
#' Turner, J.A., Bolen, C.R. & Blankenship, D.M. (2015). "Quantitative gene set analysis
#' generalized for repeated measures, confounder adjustment, and continuous covariates."
#' \emph{BMC Bioinformatics}, \emph{16}(1): 272.
#' \url{http://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-015-0707-9}
#'
#' Storey, J. & Tibshirani, R. (2003). "Statistical significance for genomewide
#' studies." \emph{Proc Natl Acad Sci}, \emph{100}(16): 9440–9445.
#' \url{http://www.pnas.org/content/100/16/9440.full}
#'
#' @examples
#' # Fit limma model
#' library(limma)
#' eset <- matrix(rnorm(5000 * 10), nrow = 5000, ncol = 10)
#' dimnames(eset) <- list(seq_len(nrow(eset)), paste0('S', seq_len(ncol(eset))))
#' clin <- data.frame(treat = rep(c("A", "B"), each = 5),
#'                       x1 = rnorm(10),
#'                       x2 = runif(10))
#' des <- model.matrix(~ treat + x1 + x2, data = clin)
#' fit <- eBayes(lmFit(eset, des))
#' 
#' # Create list of differentially expressed pathways
#' geneSets = list()
#' for (i in 0:10) {
#'   genes <- ((30 * i) + 1):(30 * (i + 1))
#'   eset[genes, clin$treat == "B"] <- eset[genes, clin$treat == "B"] + rnorm(1)
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

  # Prep data
  se <- sqrt(fit$s2.post) * fit$stdev.unscaled[, coef]
  sd_a <- se / (fit$sigma * fit$stdev.unscaled[, coef])
  sd_a[is.infinite(sd_a)] <- 1
  resid_mat <- residuals.MArrayLM(fit, dat)
  overlap <- sapply(geneSets, function(x) sum(x %in% rownames(dat)))
  geneSets <- geneSets[overlap > 0]

  # Run QuSAGE functions
  res <- newQSarray(mean = fit$coefficients[, coef],  # Create QSarray obj
                      SD = se,
                sd.alpha = sd_a,
                     dof = fit$df.total,
                  labels = rep('Resid', ncol(dat)))
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


