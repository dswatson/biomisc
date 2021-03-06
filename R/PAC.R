#' Calculate PAC Index
#'
#' This function calculates the proportion of ambiguous clustering (PAC) for
#' each value of \emph{k} tested via consensus clustering.
#'
#' @param cc A list created by a call to
#'   \code{\link[ConsensusClusterPlus]{ConsensusClusterPlus}}.
#' @param window Lower and upper bounds for the consensus index sub-interval
#'   over which to calculate the PAC. Must be on (0, 1).
#' @param plot Plot PAC scores as a bar plot? Default is \code{TRUE}.
#'
#' @details
#' Consensus clustering is a method for testing the stability of cluster
#' membership under resampling (Monti et al., 2003). Senbabaoglu et al. (2014)
#' demonstrated that traditional methods for estimating optimal cluster number
#' fail when probes are not independent, which they rarely are in omic data. The
#' authors propose a new statistic, the proportion of ambiguous clustering
#' (PAC), which measures the increase in the empirical CDF curve for each
#' potential cluster number \emph{k} over a user-defined sub-interval of the
#' consensus index generated by the consensus cluster algorithm. The minimal PAC
#' score for a given range of \emph{k} is taken to be the optimal cluster number
#' for that data set. The default settings of \code{window = c(0.1, 0.9)} are
#' taken from the original PAC paper, and generally lead to stable results.
#'
#' @return A data frame with PAC scores for each value of \emph{k} in \code{cc}.
#'
#' @references
#' Monti, S., Tamayo, P., Mesirov, J., & Golub, T. (2003).
#' \href{http://link.springer.com/article/10.1023/A:1023949509487}{Consensus
#' Clustering: A Resampling-Based Method for Class Discovery and Visualization
#' of Gene Expression Microarray Data}. \emph{Machine Learning}, \emph{52}:
#' 91-118.
#'
#' Senbabaoglu, Y., Michailidis, G. & Li, J.Z. (2014).
#' \href{http://www.nature.com/articles/srep06207}{Critical limitations of
#' consensus clustering in class discovery}. \emph{Scientific Reports}, \emph{
#' 4}:6207.
#'
#' @examples
#' # Build consensus cluster object
#' library(ConsensusClusterPlus)
#' mat <- matrix(rnorm(100 * 10), nrow = 100, ncol = 10)
#' cc <- ConsensusClusterPlus(mat, maxK = 5)
#'
#' # Run PAC()
#' PAC(cc)
#'
#' @export
#' @importFrom purrr map_lgl
#' @importFrom ggsci scale_fill_d3
#' @import dplyr
#' @import ggplot2
#'

PAC <- function(cc,
                window = c(0.1, 0.9),
                  plot = TRUE) {

  # Preliminaries
  maxK <- length(cc)
  if (!(cc %>% is.list) ||
      any(2:maxK %>% map_lgl(~ !'consensusMatrix' %in% names(cc[[.x]])))) {
    stop('cc must be a list object created by a call to ConsensClusterPlus.')
  }
  if (min(window) <= 0L || max(window) >= 1L) {
    stop('Upper and lower bounds must be on (0, 1), exclusive.')
  }

  # Tidy data
  suppressWarnings(
    PAC <- expand.grid(k = 2:maxK,
                     Idx = window) %>%
    rowwise(.) %>%
    mutate(CDF = ecdf(cc[[k]]$consensusMatrix %>% keep(lower.tri(.)))(Idx),
             k = as.factor(k)) %>%
    group_by(k) %>%
    mutate(PAC = diff(CDF)) %>%
    distinct(k, PAC) %>%
    as.data.frame(.)
  )
  
  # Export
  if (plot) {
    p <- ggplot(PAC, aes(k, PAC, fill = k)) +
      geom_bar(stat = 'identity') +
      scale_fill_d3() +
      labs(title = 'PAC Statistics') +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5))
    print(p)
  }
  cat('Optimal cluster number is k = ', PAC$k[which.min(PAC$PAC)], ', with a ',
      'PAC index of ', round(min(PAC$PAC), 6), '.\n')
  return(PAC)

}


