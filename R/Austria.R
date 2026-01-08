#' STR Allele Frequencies from Austria
#'
#' @description
#' Population allele frequency data for 16 autosomal Short Tandem Repeat
#' (STR) markers from the Austrian population. Focused on core forensic
#' markers used in European laboratories.
#'
#' @format A data frame with 66 rows (alleles) and 17 columns:
#'   \describe{
#'     \item{Allele}{Allele designation (numeric repeat number)}
#'     \item{D1S1656, D2S1338, D2S441, ...}{Allele frequencies for each marker
#'       (values between 0 and 1, summing to 1 per marker)}
#'   }
#'
#' @details
#' This dataset contains allele frequencies for the following 16 STR markers:
#' D1S1656, D2S1338, D2S441, D3S1358, D8S1179, D10S1248, D12S391, D16S539,
#' D18S51, D19S433, D21S11, D22S1045, FGA, SE33, TH01, VWA.
#'
#' These markers correspond to the European Standard Set (ESS) of forensic
#' STR loci plus commonly used additional markers.
#'
#' @usage data(Austria)
#'
#' @source
#' Austrian population frequency data. Format compatible with \pkg{pedtools}
#' and \pkg{forrel} packages.
#'
#' @references
#' Parson W, et al. (2008). "The EDNAP standardization of the NGM
#' amplification kit." \emph{Forensic Science International: Genetics
#' Supplement Series}, 1(1), 183-184. \doi{10.1016/j.fsigss.2007.10.062}
#'
#' @seealso
#' \code{\link{get_allele_freqs}} for extracting frequencies,
#' \code{\link{sim_lr_genetic}} for LR simulations.
#'
#' Other European databases: \code{\link{Europe}}, \code{\link{BosniaHerz}}
#'
#' @examples
#' # Load the dataset
#' data(Austria)
#'
#' # View structure
#' head(Austria)
#'
#' # Compare with Bosnia-Herzegovina (same marker set)
#' data(BosniaHerz)
#' identical(names(Austria), names(BosniaHerz))  # TRUE
"Austria"
