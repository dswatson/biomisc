#' Preprocess DNA methylation data
#'
#' This function runs a QN.BMIQ pipeline on raw .idat files.
#'
#' @param targets A data frame of sample information or a vector of barcodes.
#' @param idats Path to the directory containing the relevant .idat files.
#'
#' @details
#' \code{qn.bmiq} provides wrappers for functions from several packages that
#' collectively form a complete preprocessing pipeline on raw DNA methylation data.
#' It begins by filtering autosomal probes by detection \emph{p}-value, beadcount,
#' alignment redundancy, and known SNPs using the default settings of
#' \code{ChAMP::champ.load}. It then runs color bias adjustment, background correction,
#' and quantile normalization, using functions from the \code{lumi} package. Finally,
#' data are normalized using beta-mixture quantile dilation as implemented by
#' \code{wateRmelon::BMIQ}. Missing values are imputed using \emph{k}-nearest
#' neighbors with the default setting of \emph{k} = 10. Zeros are replaced with a very
#' small value to facilitate logit transformation.
#'
#' Normalization procedures for methylation data are the subject of much active
#' research. Over a dozen methods are implemented in various Bioconductor packages.
#' The QN.BMIQ pipeline has been found to be optimal in terms of both technical
#' reproducibility (Marabita et al., 2013) and biological comparison with gold
#' standard whole genome bisulfite sequencing (Wang et al., 2015).
#'
#' @return A matrix of filtered and normalized beta values.
#'
#' @references
#' Morris, T.J., Butcher, L.M., Teschendorff, A.E., Chakravarthy, A.R., Wojdacz, T.K.
#' & Beck, S. (2014). "ChAMP: 450k Chip Analysis Methylation Pipeline."
#' \emph{Bioinformatics}, \emph{30}(3): 428-430.
#' \url{http://doi.org/10.1093/bioinformatics/btt684}
#'
#' Bolstad, B.M., Irizarry, R.A. Astrand, M. & Speed, T.P. (2003). "A comparison of
#' normalization methods for high density oligonucleotide array data based on
#' variance and bias." \emph{Bioinformatics}, \emph{19}(2): 185-193.
#' \url{https://www.ncbi.nlm.nih.gov/pubmed/12538238}
#'
#' Teschendorff, A.E., Marabita, F., Lechner, M. Bartlett, T. Tegner, J. Gomez-Cabrero,
#' D. & Beck, S. (2012). "A beta-mixture quantile normalization method for correcting
#' probe design bias in Illumina Infinium 450k DNA methylation data."
#' \emph{Bioinformatics}, \emph{29}(2): 189-196.
#' \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3546795/}
#'
#' Marabita, F. et al. (2013). "An evaluation of analysis pipelines for DNA methylation
#' profiling using the Illumina HumanMethylation450 BeadChip platform."
#' \emph{Epigenetics}, \emph{8}(3): 333-346.
#' \url{https://www.ncbi.nlm.nih.gov/pubmed/23422812}
#'
#' Wang, T. et al. (2015). "A systematic study of normalization methods for Infinium
#' 450k methylation data using whole-genome bisulfite sequencing data."
#' \emph{Epigenetics}, \emph{10}(7): 662-669.
#' \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4623491/}
#'
#' @examples
#' data()
#'
#'
#' @export
#' @import lumi
#' @import wateRmelon
#' @import ChAMP
#' @import ChAMPdata
#' @import impute
#' @import dplyr
#'

qn.bmiq <- function(targets, idats) {

  # Filter
  dat <- champ.load(idats)
  keep <- rownames(dat$beta)
  targ <- read.csv(targets, stringsAsFactors = FALSE)
  dat <- importMethyIDAT(targ, idats)
  dat <- dat[keep, ]

  # QN
  dat <- lumiMethyC(obj, method = 'quantile')    # Adjust for color bias
  dat <- lumiMethyB(dat, method = 'bgAdjust2C')  # Background correct
  dat <- lumiMethyN(dat, method = 'quantile')    # Quantile normalize
  betas <- ilogit2(exprs(dat))
  rm(dat)

  # BMIQ
  data(probeInfoALL.lv)                          # Annotate probes
  probeInfo.lv <- lapply(probeInfoALL.lv, function(x) {
    x[match(rownames(betas), probeInfoALL.lv[[5]])]
  })
  design.v <- as.numeric(probeInfo.lv[[2]])
  mat <- matrix(nrow = nrow(betas), ncol = ncol(betas),
                dimnames = list(rownames(betas), colnames(betas)))
  for (j in seq_len(ncol(betas))) {            # BMIQ normalise
    bmiq.o <- BMIQ(betas[, j], design.v, plots = FALSE,
                   sampleID = colnames(betas)[j])
    mat[, j] <- bmiq.o$nbeta
  }
  mat <- impute.knn(mat)$data                    # Impute missing values
  mat[mat == 0] <- .0000001                      # Replace zeros

  return(mat)

}


