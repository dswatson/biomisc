#' Transform Counts
#'
#' This convenient wrapper function converts raw counts (or count-like data) to the
#' log2-counts per million scale, following normalization.
#'
#' @importFrom edgeR DGEList calcNormFactors cpm
#'

lcpm <- function(mat, method = 'RLE', prior.count = 0.5) {
  mat <- DGEList(mat)
  mat <- calcNormFactors(mat, method = method)
  mat <- cpm(mat, log = TRUE, prior.count = prior.count)
}
