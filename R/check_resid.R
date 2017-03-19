#' Check the residuals of a DESeq model
#'
#' This function evaluates the impact of different expression filters on the 
#' residuals of a fitted DESeqDataSet.
#'
#' @param dds A \code{DESeqDataSet} object that has been fit with a negative
#'   binomial GLM.
#' @param filt Numeric vector of length two specifying the filter criterion. Each 
#'   probe must have at least \code{filt[1]} log2-counts per million in at least 
#'   \code{filt[2]} libraries to pass the expression threshold.
#'
#' @details
#' The statistical tests upon which \code{qmod} is based presume that residuals
#' for a fitted model are approximately normally distributed. This is not generally
#' true of negative binomal GLMs, the family of models used by DESeq. The 
#' non-normality of residuals is especially pronounced for low count probes, which are
#' by default not filtered out until after modeling in the DESeq pipeline. (See 
#' \code{\link[DESeq2]{results}} for more details.) To run \code{qmod} on 
#' \code{DESeqDataSet} objects, it is necessary to filter out underexpressed probes 
#' and apply a variance stabilizing transformation. We recommend applying the lightest
#' possible expression filter, although there is no precise algorithm for determining
#' what this should be. 
#' 
#' As a general rule of thumb, the \code{limma} authors advise setting \code{filt[1]} 
#' to 10 / (\emph{L} / 1,000,000), where \emph{L} = the minimum library size for a 
#' given count matrix. They also recommend setting \code{filt[2]} to the number of 
#' replicates in the largest group. They note, however, that these are general 
#' guidelines, not strict rules. 
#'
#' @examples
#' library(DESeq2)
#' dds <- makeExampleDEESeqDataSet()
#' dds <- DESeq(dds)
#' check_resid(dds)
#'
#' @export
#' @importFrom edgeR cpm
#' @importFrom DESeq2 counts normalizationFactors sizeFactors  
#' @importFrom tidyr gather
#' @import dplyr 
#' @import ggplot2
#'

check_resid <- function(dds, 
                        filt = c(1, 1)) {
  
  # Preliminaries
  if (!is(dds, 'DESeqDataSet')) {
    stop('dds must be a DESeqDataSet.')
  }
  if (is.null(assays(dds)[['mu']])) {
    stop('dds must be fit with a negative binomial GLM.')
  }

  # Create resid_mat
  p <- nrow(dds)
  keep <- rowSums(cpm(counts(dds)) >= filt[1]) >= filt[2]
  dds <- dds[keep, , drop = FALSE]
  cnts <- counts(dds)
  if (is.null(sizeFactors(dds))) {
    cnts <- log2((cnts + 0.5) / (normalizationFactors(dds) + 1L) * 1e6L)
    fit <- log2((assays(dds)[['mu']] + 0.5) / 
                 (normalizationFactors(dds) + 1L) * 1e6L)
  } else {
    cnts <- cpm(cnts, lib.size = sizeFactors(dds), log = TRUE, prior.count = 0.5)
    fit <- cpm(assays(dds)[['mu']], lib.size = sizeFactors(dds), 
               log = TRUE, prior.count = 0.5)
  }
  resid_mat <- cnts - fit

  # Output
  bad <- p - sum(keep)
  bad_p <- round((bad / p) * 100L, 2L)
  cat(paste0('Filter criterion removes ', bad, ' (', bad_p, '%) of probes.'))
  
  # Plot
  resid_df <- gather(tbl_df(resid_mat), 'Sample', 'Residual')
  data_frame(Expected = qnorm(ppoints(1e4L)),
             Observed = quantile(resid_df$Residual, 
                                 probs = seq(0L, 1L, length.out = 1e4L))) %>%
    ggplot(aes(Expected, Observed)) + 
    geom_point(size = 0.5) + 
    labs(main = 'Normal Q-Q Plot',
         x = 'Expected Quantiles',
         y = 'Observed Quantiles') + 
    theme_bw()
  
}



