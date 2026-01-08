#' STR Allele Frequencies from China
#'
#' @description
#' Comprehensive population allele frequency data for 70 autosomal Short
#' Tandem Repeat (STR) markers from Chinese populations. This is one of
#' the most extensive STR frequency databases available.
#'
#' @format A data frame with 67 rows (alleles) and 71 columns:
#'   \describe{
#'     \item{Allele}{Allele designation (numeric repeat number)}
#'     \item{CSF1PO, D1S1656, D2S441, ...}{Allele frequencies for each marker
#'       (values between 0 and 1, summing to 1 per marker)}
#'   }
#'
#' @details
#' This comprehensive dataset contains allele frequencies for 70 STR markers,
#' including all standard forensic core loci plus an extensive set of
#' additional markers. This enables very high discrimination power for
#' identification purposes.
#'
#' Core forensic markers included: CSF1PO, D1S1656, D2S441, D2S1338,
#' D3S1358, D5S818, D7S820, D8S1179, D10S1248, D12S391, D13S317, D16S539,
#' D18S51, D19S433, D21S11, D22S1045, FGA, TH01, TPOX, vWA, SE33.
#'
#' Extended markers include: PENTA D, PENTA E, D6S1043, D4S2408, D9S1122,
#' and many others for enhanced discrimination.
#'
#' @usage data(China)
#'
#' @source
#' Chinese population frequency data. Format compatible with \pkg{pedtools}
#' and \pkg{forrel} packages.
#'
#' @references
#' Wang L, et al. (2019). "Genetic polymorphisms of 21 autosomal STR loci
#' in the Han population from Shanghai, China." \emph{Forensic Science
#' International: Genetics}, 42, e4-e6. \doi{10.1016/j.fsigen.2019.05.008}
#'
#' @seealso
#' \code{\link{get_allele_freqs}} for extracting frequencies,
#' \code{\link{sim_lr_genetic}} for LR simulations.
#'
#' Other Asian databases: \code{\link{Asia}}, \code{\link{Japan}}
#'
#' @examples
#' # Load the dataset
#' data(China)
#'
#' # This is one of the most comprehensive databases
#' ncol(China) - 1  # 70 markers
#'
#' # Check common markers with other databases
#' common <- intersect(names(China), names(Japan))
#' length(common)  # Many shared markers
"China"
