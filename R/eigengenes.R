#' Create a module eigengene matrix
#'
#' This function reduces an expression matrix to the module level.
#'
#' @param dat An expression matrix or matrix-like object, with rows corresponding to
#'   probes and columns to samples. Data is presumed to be filtered and normalized
#'   prior to dimensionality reduction. For count data, this means undergoing some
#'   sort of variance stabilizing transformation, such as \code{\link{lcpm}, 
#'   \link[DESeq2]{vst}, or \link[DESeq2]{rlog}}.
#' @param geneSets A named list of one or several gene sets.
#'
#' @details
#' Eigengenes are a convenient way to summarize data from a given gene set into a 
#' single vector of length equal to the study sample size. The eigengene of a module 
#' represents the first principle component of the sample by probe expression matrix 
#' for all genes of which the module is composed. This is useful for correlating gene 
#' sets with one another, or associating them with clinical variables of interest. 
#' They may even be used to create a meta-network of module eigengenes.
#' 
#' @references
#' Zhang, B. & Horvath, S. (2005). "A General Framework for Weighted Gene 
#' Co-Expression Network Analysis". \emph{Stat. Appl. Genet. Molec. Biol.}, \emph{4}: 
#' 1, 17.
#' \url{https://www.ncbi.nlm.nih.gov/pubmed/16646834}
#' 
#' Langfelder, P. & Horvath, S. (2007). "Eigengene networks for studying the 
#' relationships between co-expression modules." \emph{BMC Bioinformatics}, \emph{1}:
#' 54.
#' \url{http://bmcsystbiol.biomedcentral.com/articles/10.1186/1752-0509-1-54}
#'
#' Horvath S. & Dong, J. (2008). "Geometric Interpretation of Gene Coexpression Network 
#' Analysis." \emph{PLoS Comput. Biol.}, \emph{4}(8): e1000117.
#' \url{http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000117}
#'
#' @examples
#' # Simulate data
#' mat <- matrix(rnorm(5000 * 10), nrow = 5000, ncol = 10)
#' grp <- rep(c("A", "B"), each = 5)
#' geneSets = list()
#' for (i in 0:10) {
#'   genes <- ((30 * i) + 1):(30 * (i + 1))
#'   mat[genes, grp == "A"] <- mat[genes, grp == "A"] + rnorm(1)
#'   geneSets[[paste("Set", i)]] <- genes
#' }
#' 
#' # Create eigengene matrix
#' eg_mat <- eigengenes(mat, geneSets)
#'
#' @export
#' @importFrom limma getEAWP
#'

eigengenes <- function(dat,
                       geneSets) {
  
  # Preliminaries
  if (nrow(dat) < 3L) {
    stop('dat must have at least three probes.')
  }
  if (is.null(rownames(dat))) {
    rownames(dat) <- seq_len(nrow(dat))
  }
  if (is.null(colnames(dat))) {
    colnames(dat) <- paste0('Sample', seq_len(ncol(dat)))
  }
  if (!is.list(geneSets)) {
    stop('geneSets must be a list.')
  }
  if (is.null(names(geneSets))) {
    stop('geneSets must be a named list.')
  }
  overlap <- sapply(geneSets, function(g) sum(g %in% rownames(dat)))
  geneSets <- geneSets[overlap > 1L]
  if (length(geneSets) == 0L) {
    stop('No overlap detected between the genes in dat and those in geneSets.')
  }
  
  # Tidy data
  dat <- getEAWP(dat)$expr
  out <- matrix(nrow = length(geneSets), ncol = ncol(dat), 
                dimnames = list(names(geneSets), colnames(dat)))
  for (m in seq_along(geneSets)) {
    mat <- dat[rownames(dat) %in% geneSets[[m]], ]
    pca <- prcomp(t(mat))
    out[m, ] <- pca$x[, 1]
  }
  
  # Export
  return(out)
  
}

  
  