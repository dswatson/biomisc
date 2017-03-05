#' Check the distribution of residuals
#'
#' This function...
#'
#' @param dds Fitted \code{DESeqDataSet} object.
#' @param filt Numeric vector of length 2 specifying the filter criterion. Each 
#'   probe must have at least \code{filt[1]} log2-counts per million in at least 
#'   \code{filt[2]} libraries to pass the expression threshold.
#'
#' @details
#' 
#' @references
#'
#'
#' @examples
#' library(DESeq2)
#' dds <- makeExampleDEESeqDataSet()
#' dds <- DESeq(dds)
#' check_resid(dds)
#'
#' @export
#' @importFrom edgeR cpm
#' @importFrom DESeq2 counts normalizationFactors sizeFactors assays 
#' @importFrom tidyr gather
#' @importFrom dplyr tbl_df
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
  cnts <- counts(dds)
  keep <- rowSums(cpm(cnts) >= filt[1L]) >= filt[2L]
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
  gather(tbl_df(resid_mat), 'Sample', 'Expression') %>%
    ggplot(aes(Expression)) + 
    geom_path(stat = 'density') + 
    labs(main = 'Residual Density', 
         x = expression(log[2]*' CPM'),
         y = 'Density') + 
    theme_bw()
}



