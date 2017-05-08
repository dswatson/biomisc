#' Run voom on kallisto counts using tximport
#'
#' This function 
#'
#' @param txi
#' @param filter
#' @param design 
#' @param normalize.method Normalization method to be applied to the logCPM values. 
#'   Choices are the same as those for the method argument of 
#'   \code{normalizeBetweenArrays} when the data is single-channel.
#' @param span Width of the lowess smoothing window as a proportion.
#' @param plot Plot the mean-variance trend of the normalized counts?
#' @param wts Include sample weights for a heteroskedastic fit?
#' @param var.design
#' @param method
#' @param maxiter
#' @param tol
#'
#' @details
#' 
#' @references
#'
#'
#' @examples
#' 
#'
#' @export
#' @importFrom limma normalizeBetweenArrays lmFit arrayWeights voom
#' @importFrom DESeq2 DESeqDataSetFromTximport estimateSizeFactors normalizationFactors
#' @importFrom edgeR DGEList calcNormFactors cpm
#' @importFrom bioplotr plot_mv
#' @import dplyr
#' 
#' @suggests??? tximport, DESeq2
#'

voom_txi <- function(txi,
                     filter, 
                     design = NULL,
                     normalize.method = 'none',
                     span = 0.5,
                     plot = FALSE,
                     wts = FALSE,
                     var.design = NULL,
                     method = 'genebygene',
                     maxiter = 50,
                     tol = 1e-10) {
  
  out <- list()
  
  # Check counts
  
  
  #	Check design
  if (is.null(design)) {
    design <- matrix(1L, nrow = ncol(txi$counts), ncol = 1L, 
                     dimnames = list(colnames(txi$counts), 'GrandMean'))
  }
  
  # Check filter
  if (is.null(filter)) {
    keep <- rowSums(cpm(txi$counts) >= 1L) >= 1L
  } else {
    keep <- rowSums(cpm(txi$counts) >= filter[1]) >= filter[2]
  }
  
  # Use DESeq to offset counts
  dds <- DESeqDataSetFromTximport(txi, colData = as.data.frame(design), design = ~ 1)
  dds <- estimateSizeFactors(dds)
  cnts <- counts(dds, normalized = TRUE)
  
  # Plug into regular voom pipeline
  y <- cnts[keep, , drop = FALSE]
  lib.size <- colSums(y)
  y <- log2((y + 0.5) / (lib.size + 1L) * 1e6L)
  y <- normalizeBetweenArrays(y, method = normalize.method)
  fit <- lmFit(y, design, ...)
  if (is.null(fit$Amean)) fit$Amean <- rowMeans(y, na.rm = TRUE)
  
  #	Fit lowess trend to sqrt-standard-deviations by log-count-size
  sx <- fit$Amean + rowMeans(log2(lib.size + 1L)) - log2(1e6L)
  sy <- sqrt(fit$sigma)
  allzero <- rowSums(counts) == 0L
  if (any(allzero)) {
    sx <- sx[!allzero]
    sy <- sy[!allzero]
  }
  l <- lowess(sx, sy, f = span)
  
  # Plot
  if (plot) {
    ?
  }
  
  # Interpolate
  f <- approxfun(l, rule = 2)
  
  #	Find individual quarter-root fitted counts
  if (fit$rank < ncol(design)) {
    j <- fit$pivot[1:fit$rank]
    fitted.values <- fit$coef[, j, drop = FALSE] %*% t(fit$design[, j, drop = FALSE])
  } else {
    fitted.values <- fit$coef %*% t(fit$design)
  }
  fitted.cpm <- 2L^fitted.values
  fitted.count <- 1e-6L * fitted.cpm * (nf + 1)
  fitted.logcount <- log2(fitted.count)
  
  #	Apply trend to individual observations
  w <- 1L / f(fitted.logcount)^4L
  dim(w) <- dim(fitted.logcount)
  
  #	Output
  out$E <- y
  out$weights <- w
  out$design <- design
  out$targets <- nf
  
  v <- new('EList', out)
  if (!wts) return(v)
  else {
    aw <- arrayWeights(v, design = design, var.design = var.design, method = method, 
                       maxiter = maxiter, tol = tol)
    v <- voom(cnts, design = design, weights = aw, lib.size = lib.size, 
              normalize.method = normalize.method, plot = plot, span = span, ...)
    aw <- arrayWeights(v, design = design, var.design = var.design, method = method, 
                       maxiter = maxiter, tol = tol, trace = trace)
    v$weights <- t(aw * t(v$weights))
    v$sample.weights <- aw
    return(v)
  }
  
}


