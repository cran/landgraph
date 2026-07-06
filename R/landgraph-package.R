#' landgraph: Graphs and Covariance for Landscape Genetics
#'
#' \if{html}{\figure{logo.png}{options: style='float: right' alt='logo' width='120'}}
#'
#' Shared, dependency-light building blocks for landscape-genetic network methods.
#' A lightweight deme/landscape graph (\code{\link{deme_graph}}); genetic
#' covariance and distance from biallelic or multiallelic data
#' (\code{\link{cov_from_biallelic}}, \code{\link{cov_from_genetic_data}},
#' \code{\link{fst_from_biallelic}}, \code{\link{dist_from_cov}},
#' \code{\link{dist_from_biallelic}}); and antisymmetric per-edge directional
#' covariate builders (\code{\link{edge_gradient}}, \code{\link{edge_flow}}). Used
#' by the symmetric (terradish) and asymmetric (dragonflow) network methods.
#'
#' @keywords internal
#' @importFrom stats setNames
"_PACKAGE"
NULL
