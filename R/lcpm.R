#' Transform Counts
#'
#' This convenient wrapper function converts raw counts to the log2-counts per million 
#' scale, following normalization and a minimal count shift.
#' 
#' @param mat Probe by sample matrix of raw counts.
#' @param method Normalization method to be used. See Details.
#' @param prior.count Average count to be added to each observation to avoid taking 
#'   log of zero.
#' 
#' @details 
#' \code{method = "TMM"} is the weighted trimmed mean of M-values (to the reference) 
#' proposed by Robinson & Oshlack (2010), where the weights are from the delta method 
#' on Binomial data. If \code{refColumn} is unspecified, the library whose upper 
#' quartile is closest to the mean upper quartile is used.
#' 
#' \code{method = "RLE"} is the scaling factor method proposed by Anders & Huber 
#' (2010). We call it "relative log expression", as median library is calculated from 
#' the geometric mean of all columns and the median ratio of each sample to the median 
#' library is taken as the scale factor.
#' 
#' \code{method = "upperquartile"} is the upper-quartile normalization method of 
#' Bullard, et al. (2010), in which the scale factors are calculated from the 75% 
#' quantile of the counts for each library, after removing genes which are zero in all 
#' libraries. This idea is generalized here to allow scaling by any quantile of the 
#' distributions.
#' 
#' If \code{method = "none"}, then the normalization factors are set to 1.
#' 
#' For symmetry, normalization factors are adjusted to multiply to 1. The effective 
#' library size is then the original library size multiplied by the scaling factor.
#' 
#' Note that rows that have zero counts for all columns are trimmed before 
#' normalization factors are computed. Therefore rows with all zero counts do not 
#' affect the estimated factors.
#' 
#' @return A numeric matrix of normalized counts on the log2-CPM scale. 
#' 
#' @references 
#' Anders, S. & Huber, W. (2010). "Differential expression analysis for sequence 
#' count data." \emph{Genome Biology}, 11:R106.
#' \url{https://genomebiology.biomedcentral.com/articles/10.1186/gb-2010-11-10-r106}
#' 
#' Bullard, J.H., Purdom, E., Hansen, K.D. & Dudoit, S. (2010). "Evaluation of 
#' statistical methods for normalization and differential expression in mRNA-Seq 
#' experiments." \emph{BMC Bioinformatics}, 11:94.
#' \url{http://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-94}
#' 
#' Robinson, M.D. & Oshlack, A. (2010). "A scaling normalization method for 
#' differential expression analysis of RNA-seq data." \emph{Genome Biology}, 11:R25.
#' \url{https://genomebiology.biomedcentral.com/articles/10.1186/gb-2010-11-3-r25}
#' 
#' @examples
#' # Simulate count data
#' mat <- matrix(rnbinom(5000, mu = 4, size = 1), nrow = 1000, ncol = 5)
#' 
#' # Plot raw counts
#' library(bioplotr)
#' plot_density(mat)
#' 
#' # Plot transformed counts
#' trans_mat <- lcpm(mat)
#' plot_density(trans_mat)
#'
#' @export
#' @importFrom edgeR DGEList calcNormFactors cpm
#'

lcpm <- function(mat, 
                 method = 'RLE', 
            prior.count = 0.5) {
  
  # Transform
  mat <- DGEList(mat)
  mat <- calcNormFactors(mat, method = method)
  mat <- cpm(mat, log = TRUE, prior.count = prior.count)
  
  # Export
  return(mat)
}
