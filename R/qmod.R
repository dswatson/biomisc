#' Run QuSAGE on limma or DESeq2 model objects
#'
#' This function is a wrapper for the QuSAGE algorithm, which tests for pathway 
#' enrichment, designed for easy integration with limma and DESeq2.
#'
#' @param fit An object of class \code{limma::\link[limma]{MArrayLM}}, as created by 
#'   a call to \code{\link[limma]{eBayes}}, or a \code{\link[DESeq2]{DESeqDataSet}} 
#'   that has been fit with a negative binomial GLM. See Details.
#' @param dat An expression matrix or matrix-like object, with rows corresponding to
#'   probes and columns to samples. Only necessary if \code{fit} is an \code{MArrayLM}
#'   object.
#' @param coef Column name or number specifying which coefficient of the model is of
#'   interest. 
#' @param contrast Character or numeric vector of length two, specifying the column
#'   names or numbers to be contrasted. The first and second elements will be the 
#'   numerator and denominator, respectively, of the fold change calculation. 
#'   Exactly one of \code{coef} or \code{contrast} must be \code{NULL}. 
#' @param geneSets Either a named list of pathways to be compared, or a vector of
#'   probe names representing a single gene set.
#' @param n.points The number of points at which to sample the convoluted
#'   \emph{t}-distribution. See Details.
#'
#' @details
#' ### INTRO TO QMOD? ###
#' 
#' ### SOMETHING ON MArrayLM vs. DESeqDataSet OBJECTS ###
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
#' Smyth, G.K. (2004). "Linear models and empirical Bayes methods for assessing
#' differential expression in microarray experiments." \emph{Stat. Appl. Genet. 
#' Molec. Biol.}, \emph{3}(1).
#' \url{http://www.statsci.org/smyth/pubs/ebayes.pdf}
#' 
#' Love, M., Huber, W., & Anders, S. (2014). "Moderated estimation of fold change and 
#' dispersion for RNA-seq data with DESeq2." \emph{Genome Biology}, \emph{15}:550.
#' \url{https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8}
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
#' qlim(fit, eset, coef = 2, geneSets)
#'
#' @export
#' @importFrom limma eBayes getEAWP 
#' @import DESeq2 
#' @import qusage
#' @import qvalue
#' @import dplyr
#'

qmod <- function(fit,
                 dat = NULL,
                 coef,
                 contrast = NULL,
                 geneSets,
                 n.points = 2^14) {

  # Preliminaries
  if (nrow(fit) < 3) {
    stop('fit must have at least three probes.')
  }
  if (is(fit, 'MArrayLM')) {
    coefs <- colnames(fit$coefficients)
    p <- length(coefs)
    if (is.null(dat)) {
      stop('dat must be provided when fit is an MArrayLM object.')
    }
    if (nrow(fit) != nrow(dat) || nrow(fit$design) != ncol(dat)) {
      stop('dat is not conformal with fit.')
    }
    if (!identical(rownames(fit), rownames(dat))) {
      stop('dat and fit must have identical rownames.')
    }
    if (is.null(rownames(fit)) && is.null(rownames(dat))) {
      rownames(fit) <- rownames(dat) <- seq_len(nrow(fit))
    } 
  } else if (is(fit, 'DESeqDataSet')) {
    coefs <- resultsNames(fit)
    p <- ncol(model.matrix(design(fit), colData(fit)))
    if (is.null(rownames(fit))) {
      rownames(fit) <- seq_len(nrow(fit))
    }
    if (!is.null(dat)) {
      warning('dat is ignored when fit is a DESeqDataSet.')
    }
  } else {
    stop('fit must be an object of class MArrayLM or DESeqDataSet.')
  }
  if (is.null(contrast)) {
    if (is.character(coef) && !coef %in% coefs) {
      stop(paste0("'", coef, "' not found in fit's design matrix."))
    } 
    if (is.numeric(coef) && !coef %in% seq_len(p)) {
      stop(paste("No coef number", coef, "found in fit's design matrix."))
    }
  } else if (is.null(coef)) {
    if (length(contrast) != 2) {
      stop('contrast must be a vector of length 2.')
    }
    if ((is.character(contrast) && any(!contrast %in% coef)) ||
        (is.numeric(contrast) && any(!contrast %in% seq_len(p)))) {
      stop("Both coefficients passed to contrast must be in fit's design matrix.")
    }
  } else {
    stop('Exactly one of coef or contrast must be NULL.')
  }
  if (is(fit, 'MArrayLM')) {
    if (is.null(fit$t) && is.null(fit$F) && is.null(contrast)) {
      fit <- eBayes(fit)
      warning('Standard errors for fit had not been moderated. Running eBayes before ',
              'testing for enrichment. See "?eBayes" for more info.')
    }
    if (!is.null(fit$t) && !is.null(fit$F) && is.null(coef)) {
      stop('Standard errors for fit must not be moderated when passing a contrast ',
           'to qmod. Use an lmFit output instead. The function will internally ',
           'create the appropriate contrast matrix and run eBayes on that. See ',
           '"?contrasts.fit" for more info.')
    }
  }
  
  # Prep data
  if (is(fit, 'MArrayLM')) {
  
    dat <- getEAWP(dat)
    dat <- dat$exprs
    resid_mat <- residuals(fit, dat)
    if (is.null(coef)) {
      coef <- 'Contrast'
      suppressWarnings(
        cm <- makeContrasts(coef = paste(contrast[2], '-', contrast[1]), 
                            levels = coefs)
      )
      fit <- contrasts.fit(fit, cm)
      fit <- eBayes(fit)
    } 
    mean <- fit$coefficients[, coef] 
    SD <- sqrt(fit$s2.post) * fit$stdev.unscaled[, coef]
    sd.alpha <- se / (fit$sigma * fit$stdev.unscaled[, coef])
    sd.alpha[is.infinite(sd.alpha)] <- 1
    dof <- fit$df.total
    
  } else {
    
    if (is.null(contrast)) {
      dds_res <- results(fit, name = coef)
    } else {
      dds_res <- results(fit, contrast = list(contrast))
    }
    dds_res <- as.data.frame(dds_res) %>% na.omit()
    fit <- fit[rownames(dds_res), ]
    if (is.null(sizeFactors(fit))) {
      signal_mat <- assays(fit)[['mu']]
    } else {
      signal_mat <- t(t(assays(fit)[['mu']]) / sizeFactors(fit))
    }
    resid_mat <- counts(fit, normalized = TRUE) - signal_mat
    mean <- dds_res$log2FoldChange
    SD <- dds_res$lfcSE
    sd.alpha <- rep(1, nrow(fit))
    dof <- ncol(fit) - p
    
  }
  overlap <- sapply(geneSets, function(p) sum(p %in% rownames(fit)))
  geneSets <- geneSets[overlap > 0]

  # Run QuSAGE functions
  res <- newQSarray(mean = mean,                     # Create QSarray obj
                      SD = SD,
                sd.alpha = sd.alpha,
                     dof = dof,
                  labels = rep('resid', ncol(fit)))
  res <- aggregateGeneSet(res, geneSets, n.points)   # PDF per gene set
  res <- calcVIF(resid_mat, res, useCAMERA = FALSE)  # VIF on resid_mat

  # Export
  out <- qsTable(res, number = Inf, sort.by = 'p') %>%
    rename(p.value = p.Value,
           Pathway = pathway.name,
             logFC = log.fold.change) %>%
    mutate(q.value = qvalue(p.value)$qvalues) %>%
    select(Pathway:p.value, q.value)
  return(out)

}

# Add arguments to tweak internal calls to limma::eBayes and/or DESeq::results functions?
