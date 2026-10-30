#' F_ST from biallelic allele counts
#'
#' Estimates pairwise F\eqn{_{ST}} from counts of the derived allele across
#' biallelic markers using the ratio-of-averages estimator of Bhatia et al.
#' (2013).
#'
#' @param Y Numeric matrix of derived-allele counts with rows as populations
#'   and columns as loci.
#' @param N Haploid sample sizes in any form accepted by
#'   \code{\link{cov_from_biallelic}}; every entry must exceed 1. For diploid
#'   data without missing genotypes, \code{N = 2 * n_individuals} for each
#'   cell.
#'
#' @details
#' The estimator computes, for each locus, the squared allele-frequency
#' difference between populations, corrects for within-population
#' heterozygosity, and then forms a ratio of locus-averages (Bhatia et al.
#' 2013).  Diagonal entries are set to 0.
#'
#' The result is not guaranteed to lie in \eqn{[0, 1]} for every pair; small
#' negative values can arise from sampling noise in very similar populations.
#'
#' To use F\eqn{_{ST}} as the response in \code{terradish}, pass the
#' matrix directly to the formula left-hand side with
#' \code{terradish::mlpe()}.
#' This ratio estimator is not a Wishart response and has no Wishart degrees
#' of freedom. Use covariance-derived squared distances for Wishart models.
#'
#' @references
#' Bhatia G, Patterson N, Sankararaman S, Price AL. 2013. Estimating and
#' interpreting F\eqn{_{ST}}: the impact of rare variants. Genome Research
#' 23(9):1514-1521. \doi{10.1101/gr.154831.113}
#'
#' @return A symmetric matrix named by \code{rownames(Y)} holding pairwise
#'   Hudson F\eqn{_{ST}}, with a zero diagonal.
#'
#' @seealso \code{\link{cov_from_biallelic}}, \code{\link{dist_from_biallelic}}
#'
#' @examples
#' Y <- matrix(c(2, 1, 0,
#'               1, 1, 1,
#'               0, 1, 2), nrow = 3, byrow = TRUE)
#' N <- matrix(2, nrow = 3, ncol = 3)
#' fst_from_biallelic(Y, N)
#'
#' # A scalar or a locus-length vector also works
#' fst_from_biallelic(Y, 2)
#'
#' @export

fst_from_biallelic <- function(Y, N)
{
  if (is.data.frame(Y))
    Y <- as.matrix(Y)
  if (!is.matrix(Y) || !is.numeric(Y) || nrow(Y) < 2L || ncol(Y) < 1L)
    stop("`Y` must be a numeric matrix or data frame with at least two rows and one column",
         call. = FALSE)
  storage.mode(Y) <- "double"

  # Expand scalar and vector `N` so a locus-length vector aligns with columns.
  N <- .biallelic_sample_size_matrix(N, Y, ploidy = 2)

  if (anyNA(Y) || anyNA(N))
    stop("missing values are not currently supported")
  if (any(Y<0) || any(N<0))
    stop("Cannot have negative counts")
  if (any(Y>N))
    stop("Allele counts cannot exceed number of haploids sampled")
  if (any(N <= 1))
    stop("All sample sizes in `N` must exceed 1")

  Fr <- Y/N
  # Mean squared frequency difference over loci, without an n x n x L array.
  r <- rowSums(Fr^2)
  f2 <- (outer(r, r, "+") - 2 * tcrossprod(Fr)) / ncol(Fr)
  cornum <- rowMeans(Fr*(1-Fr)/(N-1))
  corden <- rowMeans(Fr*(1-Fr)*N/(N-1))
  fst <- (f2 - c(outer(cornum, cornum, "+")))/(f2 - c(outer(cornum, cornum, "+")) + c(outer(corden, corden, "+")))
  fst <- matrix(fst, nrow(Y), nrow(Y),
                dimnames = list(rownames(Y), rownames(Y)))
  diag(fst) <- 0
  fst
}

#' Genetic covariance from biallelic allele counts
#'
#' Computes a sample covariance matrix from counts of the derived allele across
#' biallelic (SNP) markers.  Rows may be individuals, populations, or any other
#' sampled genetic unit.
#'
#' @param Y Numeric matrix of derived-allele counts with sampled units in rows
#'   and loci in columns.  For diploid individuals genotyped as 0/1/2, this is
#'   the standard dosage matrix.
#' @param N Optional haploid sample-size specification.  Controls how many
#'   chromosomes were sampled at each locus for each unit:
#'   \itemize{
#'     \item \code{NULL}: use \code{ploidy} for all cells.
#'     \item A scalar: recycled to all cells.
#'     \item A vector of length \code{nrow(Y)}: row-specific sample sizes
#'       (same at every locus for that row).
#'     \item A vector of length \code{ncol(Y)}: locus-specific sample sizes
#'       (same across all rows).
#'     \item A matrix matching \code{dim(Y)}: fully general cell-level sizes.
#'   }
#' @param ploidy Haploid chromosome count per unit, used when \code{N = NULL}.
#'   The default \code{2} corresponds to diploid individuals with dosage coding
#'   0/1/2 (so each individual contributes two chromosomes at each locus).
#' @param monomorphic How to handle loci with pooled allele frequency 0 or 1
#'   (fixed across all rows).  \code{"drop"} (default) removes them and emits a
#'   warning; \code{"error"} stops with an error to preserve historical strict
#'   behavior.
#' @param tol Tolerance for identifying monomorphic loci.  A locus is treated
#'   as monomorphic when its pooled frequency is within \code{tol} of 0 or 1.
#'
#' @details
#' Let \eqn{p_\ell = \sum_i Y_{i\ell} / \sum_i N_{i\ell}} be the pooled
#' derived-allele frequency at locus \eqn{\ell}.  The normalized allele dosage
#' is:
#'
#' \deqn{Z_{i\ell} = \frac{Y_{i\ell} - N_{i\ell} p_\ell}
#'                        {\sqrt{N_{i\ell} p_\ell (1 - p_\ell)}}}
#'
#' The function returns \eqn{Z Z^\top / L}, where \eqn{L} is the number of
#' retained (polymorphic) loci.  The resulting matrix is the sample covariance
#' in normalized allele-frequency space, analogous to the genomic relationship
#' matrix of Yang et al. (2010).
#'
#' The function standardizes counts, not frequencies, so entries grow with the
#' haploid sample size: \eqn{Z_{i\ell} = \sqrt{N_{i\ell}}(\hat p_{i\ell} -
#' p_\ell)/\sqrt{p_\ell(1-p_\ell)}}. For individual 0/1/2 dosage this is the
#' genomic relationship matrix scale. For population counts with equal sample
#' sizes, \code{cov_from_biallelic(Y, N) / N[1]} gives the standardized
#' allele-frequency covariance that terradish fits. Do not use that shortcut
#' with unequal sizes; rows then carry unequal sampling variance that a common
#' nugget does not absorb. \code{Y} and \code{N} cannot contain \code{NA}. For
#' population counts, record missing genotypes by lowering the matching cells
#' of \code{N}. Every cell of \code{N} must be positive, so for individual
#' dosage you must filter loci or individuals first. Counting either allele
#' gives the same covariance.
#'
#' Per-locus standardization gives greater weight to rare variants, including
#' their contribution to the diagonal. Filter loci by minor allele frequency
#' before constructing the covariance, report the filtering rule, and check
#' that conclusions hold across thresholds. Unequal chromosome counts also
#' imply unequal sampling variance. With unequal counts across rows,
#' pooled-frequency centering does not generally produce zero row sums in the
#' covariance.
#'
#' To convert the covariance to a pairwise distance matrix, use
#' \code{\link{dist_from_cov}}.  For non-biallelic or multiallelic data, use
#' \code{\link{cov_from_genetic_data}}.
#'
#' @references
#' Yang J, Benyamin B, McEvoy BP, Gordon S, Henders AK, Nyholt DR,
#' Madden PA, Heath AC, Martin NG, Montgomery GW, Goddard ME, Visscher PM.
#' 2010. Common SNPs explain a large proportion of the heritability for human
#' height. Nature Genetics 42(7):565-569. \doi{10.1038/ng.608}
#'
#' @return A symmetric positive-semidefinite \code{nrow(Y)} by \code{nrow(Y)}
#'   matrix of pairwise genetic covariances on the count-standardized scale
#'   (see Details).  Row and column names match the row names of \code{Y}
#'   when present. The \code{diagonal} attribute is \code{"normalized_dosage"}.
#'   The \code{centered} attribute is \code{"sites"} when chromosome counts
#'   are equal across rows at every retained locus; otherwise it is
#'   \code{"pooled_allele_frequency"}. These labels describe the construction
#'   without changing the returned numerical values.
#'
#' @seealso \code{\link{cov_from_genetic_data}}, \code{\link{fst_from_biallelic}},
#'   \code{\link{dist_from_biallelic}}, \code{\link{dist_from_cov}},
#'   \code{terradish::wishart_covariance()}
#'
#' @examples
#' # Individual diploid dosage matrix: 3 individuals, 3 SNPs
#' Y <- matrix(c(2, 1, 0,
#'               1, 1, 1,
#'               0, 1, 2), nrow = 3, byrow = TRUE)
#' cov_from_biallelic(Y)          # ploidy = 2 by default
#'
#' # Population-level counts with unequal sampling
#' Y2 <- matrix(c(10, 4, 1,
#'                5,  5, 3), nrow = 2, byrow = TRUE)
#' N2 <- matrix(c(20, 20, 10,
#'                10, 10, 6), nrow = 2, byrow = TRUE)
#' cov_from_biallelic(Y2, N = N2)
#'
#' @export

cov_from_biallelic <- function(Y,
                               N = NULL,
                               ploidy = 2,
                               monomorphic = c("drop", "error"),
                               tol = sqrt(.Machine$double.eps))
{
  monomorphic <- match.arg(monomorphic)

  if (is.data.frame(Y))
    Y <- as.matrix(Y)
  if (!is.matrix(Y) || !is.numeric(Y) || nrow(Y) < 2L || ncol(Y) < 1L)
    stop("`Y` must be a numeric matrix or data frame with at least two rows and one column",
         call. = FALSE)
  storage.mode(Y) <- "double"

  N <- .biallelic_sample_size_matrix(N, Y, ploidy)

  if (anyNA(Y) || anyNA(N))
    stop("missing values are not currently supported")
  if (any(!is.finite(Y)) || any(!is.finite(N)))
    stop("`Y` and `N` must contain only finite values", call. = FALSE)
  if (any(Y<0) || any(N<0))
    stop("Cannot have negative counts")
  if (any(Y>N))
    stop("Allele counts cannot exceed number of haploids sampled")
  if (any(N <= 0))
    stop("All sample sizes in `N` must be positive")
  if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol < 0)
    stop("`tol` must be a finite nonnegative number", call. = FALSE)

  Fr <- colSums(Y) / colSums(N)
  variable <- is.finite(Fr) & Fr > tol & Fr < (1 - tol)
  if (!all(variable))
  {
    if (identical(monomorphic, "error"))
      stop("Each retained locus must have pooled allele frequency strictly between 0 and 1",
           call. = FALSE)
    if (!any(variable))
      stop("No variable biallelic loci remain after dropping monomorphic loci",
           call. = FALSE)
    dropped <- colnames(Y)[!variable]
    if (is.null(dropped))
      dropped <- which(!variable)
    warning("Dropping monomorphic biallelic loci: ",
            paste(dropped, collapse = ", "),
            call. = FALSE)
    Y <- Y[, variable, drop = FALSE]
    N <- N[, variable, drop = FALSE]
    Fr <- Fr[variable]
  }

  # Fr is a per-locus (length-ncol) vector; broadcast it across columns so each
  # locus frequency aligns with its column. (Bare `N*Fr` would recycle Fr
  # column-major, down the rows, misaligning frequencies when nrow(Y) > 1.)
  Frm <- matrix(Fr, nrow(Y), ncol(Y), byrow = TRUE)
  Y  <- (Y - N * Frm) / sqrt(N * Frm * (1 - Frm))

  covariance <- Y %*% t(Y) / ncol(Y)
  attr(covariance, "diagonal") <- "normalized_dosage"
  # The denominator preserves site centering only for equal counts per locus.
  equal_counts <- all(vapply(seq_len(ncol(N)), function(j)
    all(N[, j] == N[1L, j]), logical(1)))
  attr(covariance, "centered") <- if (equal_counts) "sites" else
    "pooled_allele_frequency"
  covariance
}

.biallelic_sample_size_matrix <- function(N, Y, ploidy)
{
  if (is.null(N))
  {
    if (!is.numeric(ploidy) || length(ploidy) != 1L ||
        !is.finite(ploidy) || ploidy <= 0)
      stop("`ploidy` must be one finite positive number", call. = FALSE)
    N <- ploidy
  }

  if (is.data.frame(N))
    N <- as.matrix(N)
  if (is.matrix(N))
  {
    if (!all(dim(Y) == dim(N)))
      stop("Dimension mismatch", call. = FALSE)
    if (!is.numeric(N))
      stop("`N` must be numeric", call. = FALSE)
    storage.mode(N) <- "double"
    return(N)
  }
  if (!is.numeric(N))
    stop("`N` must be numeric", call. = FALSE)
  if (length(N) == 1L)
    return(matrix(N, nrow(Y), ncol(Y), dimnames = dimnames(Y)))
  if (length(N) == nrow(Y))
    return(matrix(N, nrow(Y), ncol(Y), dimnames = dimnames(Y)))
  if (length(N) == ncol(Y))
    return(matrix(N, nrow(Y), ncol(Y), byrow = TRUE, dimnames = dimnames(Y)))

  stop("`N` must be a scalar, a matrix matching `Y`, or a row- or locus-length vector",
       call. = FALSE)
}

#' Covariance from multivariate genetic data
#'
#' Calculates individual- or population-level genetic covariance from
#' multivariate genetic data. This is useful for microsatellite, multiallelic,
#' SNP dosage, PCA score, or other numeric encodings.
#'
#' @param x Genetic data. For \code{input = "features"}, a numeric matrix or data
#'   frame with individuals in rows and genetic features in columns. For
#'   \code{input = "allele_calls"}, a matrix or data frame of allele calls with
#'   one column per allele copy.
#' @param groups Optional factor, character, or integer vector assigning each row
#'   of \code{x} to a population. If \code{NULL}, each row is treated as an
#'   individual sampled unit and an individual-level covariance matrix is
#'   returned.
#' @param input Input format. \code{"features"} treats \code{x} as an already
#'   numeric feature matrix. \code{"allele_calls"} converts allele calls to
#'   per-locus allele dosage columns before calculating covariance.
#' @param loci Required for \code{input = "allele_calls"}; a vector with one
#'   entry per column of \code{x} identifying the locus for each allele-copy
#'   column. For a diploid microsatellite matrix with two columns per locus, the
#'   two columns for each locus should have the same \code{loci} value.
#' @param center Should feature columns be centered by their global mean before
#'   covariance calculation? Default \code{TRUE}. The covariance depends only
#'   on differences among units, so centering changes only the
#'   \code{centroids}, \code{unit_features}, and \code{feature_center}
#'   attributes.
#' @param scale Should feature columns be divided by their global standard
#'   deviation? Default \code{TRUE}. Scaling divides each feature by its
#'   standard deviation, so rare alleles get the same variance as common ones
#'   and therefore more relative weight. Constant features are dropped with a
#'   warning.
#' @param tol Tolerance used to identify constant features when
#'   \code{scale = TRUE}.
#' @param diagonal How to set the returned covariance diagonal.
#'   \code{"gower"} (default) keeps the double-centered diagonal, so every
#'   entry shares one scale; use it for Wishart models. \code{"within"}
#'   substitutes summed within-population variances, as in population graphs;
#'   terradish warns when it receives one. \code{"auto"} restores the default
#'   used before 0.0.3: \code{"within"} for grouped data with replication and
#'   \code{"gower"} otherwise.
#' @param normalize Should covariance scale be normalized? \code{"none"} returns
#'   sums over retained features. \code{"features"} divides the covariance,
#'   squared distances, and within-population variances by the number of retained
#'   features. The divisor counts expanded allele features, not loci, and is
#'   not \code{nu}.
#'
#' @details
#' The function first obtains a numeric feature matrix. Microsatellite data can
#' be supplied as allele calls by setting \code{input = "allele_calls"}; these are
#' converted to allele dosage columns, one column per locus-allele combination.
#'
#' Wishart models in terradish, \code{terradish::wishart_covariance()} and
#' \code{terradish::generalized_wishart()}, need \code{nu}, an
#' effective-information value that you supply. This function does not compute
#' it, and the number of SNPs, loci, or allele features does not supply it:
#' linked loci share population history, so they carry far less independent
#' information than their count suggests. Hold \code{nu} fixed across compared
#' models and report sensitivity with \code{terradish::terradish_rescale_nu()}.
#'
#' If \code{groups = NULL}, rows of the processed feature matrix are used
#' directly and squared Euclidean distances among individuals are transformed to
#' a covariance matrix by Gower double-centering. This produces a positive
#' semidefinite row-level covariance matrix up to numerical tolerance.
#'
#' If \code{groups} is supplied, population centroids are calculated in feature
#' space before the Gower transform. With \code{diagonal = "within"}, the
#' diagonal is replaced by the sum of within-population feature variances. This
#' matches the covariance construction used before population-graph
#' partial-correlation filtering in Dyer-style population graph workflows.
#'
#' A numeric feature matrix cannot contain \code{NA}. For allele calls, each
#' missing copy becomes the most frequent allele at its locus; the
#' \code{imputed_allele_calls} attribute records the count. Modal imputation
#' shrinks distances, so report it. \code{terradish::pca_dist()} imputes allele
#' means instead. Code missing calls as \code{NA}: the function reads codes such
#' as 0 or -9 as real alleles.
#'
#' The within-population diagonal is on a different scale from the covariance
#' among population centroids and is unsuitable for a Wishart likelihood.
#' Use \code{diagonal = "gower"} for that purpose. Replacing the diagonal
#' also removes the guarantee of positive semidefiniteness and zero row sums.
#' Unequal group sizes imply unequal sampling variance. Standardizing allele
#' features gives greater weight to rare variants; filter by minor allele
#' frequency before constructing the covariance and report the filtering rule.
#'
#' @references
#' Gower JC. 1966. Some distance properties of latent root and vector methods
#' used in multivariate analysis. Biometrika 53(3-4):325-338.
#' \doi{10.1093/biomet/53.3-4.325}
#'
#' Dyer RJ, Nason JD. 2004. Population graphs: the graph theoretic shape of
#' genetic structure. Molecular Ecology 13(7):1713-1727.
#' \doi{10.1111/j.1365-294X.2004.02177.x}
#'
#' Dyer RJ. 2015. Population graphs and landscape genetics. Annual Review of
#' Ecology, Evolution, and Systematics 46:327-342.
#' \doi{10.1146/annurev-ecolsys-112414-054150}
#'
#' @return A symmetric numeric matrix with one row and column per individual
#'   (when \code{groups = NULL}) or population (when \code{groups} is
#'   supplied).  The matrix is positive semi-definite up to numerical
#'   tolerance, except when \code{diagonal = "within"} is used (see Details).
#'   The following attributes are attached. For a result \code{S}, use
#'   \code{S[, ]} to view the matrix without attributes.
#' \describe{
#'   \item{\code{centroids}, \code{unit_features}}{Feature matrix of
#'     population centroids (or individual feature vectors when ungrouped),
#'     one row per returned unit. The two attributes hold the same matrix.}
#'   \item{\code{centroid_distance2}, \code{unit_distance2}}{Squared Euclidean
#'     distance matrix among the rows of \code{centroids}, divided by
#'     \code{normalizer}. The two attributes hold the same matrix.}
#'   \item{\code{within_variance}}{Named vector of summed within-population
#'     feature variances, divided by \code{normalizer}. Attached whenever
#'     \code{groups} is supplied, with 0 for single-member groups. It fills
#'     the diagonal only when \code{diagonal = "within"}.}
#'   \item{\code{unit_size}}{Named vector of group sizes: numeric when
#'     \code{groups} is supplied, integer (all 1) when ungrouped.}
#'   \item{\code{feature_center}, \code{feature_scale}}{Centering and scaling
#'     constants applied to each retained feature.}
#'   \item{\code{retained_features}}{Names of the features kept after constant
#'     features were dropped.}
#'   \item{\code{input}}{The \code{input} format used.}
#'   \item{\code{imputed_allele_calls}, \code{imputed_modal_alleles}}{Only for
#'     \code{input = "allele_calls"}: the number of imputed calls per locus
#'     and the modal allele used to fill them, for loci with missing calls
#'     (zero-length when no call was missing).}
#'   \item{\code{level}}{Either \code{"individual"} or \code{"population"}.}
#'   \item{\code{diagonal}}{The \code{diagonal} method actually used.}
#'   \item{\code{centered}}{\code{"sites"} for the Gower covariance, or
#'     \code{"sites_before_diagonal_replacement"} when the within-population
#'     diagonal replaces the centered covariance diagonal.}
#'   \item{\code{normalize}}{The \code{normalize} option used.}
#'   \item{\code{normalizer}}{Divisor applied to the covariance and distances
#'     (\code{1} for \code{normalize = "none"}, number of retained features for
#'     \code{normalize = "features"}).}
#' }
#'
#' @seealso \code{terradish::wishart_covariance()},
#'   \code{\link{cov_from_biallelic}}, \code{\link{dist_from_cov}},
#'   \code{terradish::simulate_covariance_response()}
#'
#' @examples
#' # Numeric allele dosage / feature matrix
#' x <- matrix(c(0, 1,
#'               1, 1,
#'               2, 0,
#'               2, 1,
#'               0, 2,
#'               1, 2),
#'             ncol = 2, byrow = TRUE)
#' cov_from_genetic_data(x)
#'
#' groups <- rep(c("pop1", "pop2", "pop3"), each = 2)
#' cov_from_genetic_data(x, groups = groups)
#'
#' # Microsatellite-style allele-call columns: two allele copies per locus
#' alleles <- data.frame(
#'   loc1_a = c(100, 100, 102, 102, 104, 104),
#'   loc1_b = c(100, 102, 102, 104, 104, 100),
#'   loc2_a = c(200, 202, 200, 202, 204, 204),
#'   loc2_b = c(202, 202, 204, 204, 204, 200)
#' )
#' cov_from_genetic_data(
#'   alleles,
#'   groups = groups,
#'   input = "allele_calls",
#'   loci = c("loc1", "loc1", "loc2", "loc2")
#' )
#'
#' @export
cov_from_genetic_data <- function(x,
                                  groups = NULL,
                                  input = c("features", "allele_calls"),
                                  loci = NULL,
                                  center = TRUE,
                                  scale = TRUE,
                                  tol = sqrt(.Machine$double.eps),
                                  diagonal = c("gower", "within", "auto"),
                                  normalize = c("none", "features"))
{
  input <- match.arg(input)
  diagonal <- match.arg(diagonal)
  normalize <- match.arg(normalize)

  if (is.data.frame(x))
    x <- as.matrix(x)
  if (!is.matrix(x) || nrow(x) < 2L || ncol(x) < 1L)
    stop("`x` must be a matrix or data frame with at least two rows and one column",
         call. = FALSE)
  if (!is.null(groups))
  {
    if (length(groups) != nrow(x))
      stop("`groups` must have one entry per row of `x`", call. = FALSE)
    if (anyNA(groups))
      stop("missing values are not supported in `groups`", call. = FALSE)
  }
  if (!is.logical(center) || length(center) != 1L || is.na(center))
    stop("`center` must be TRUE or FALSE", call. = FALSE)
  if (!is.logical(scale) || length(scale) != 1L || is.na(scale))
    stop("`scale` must be TRUE or FALSE", call. = FALSE)
  if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol < 0)
    stop("`tol` must be a finite nonnegative number", call. = FALSE)

  features <- if (identical(input, "allele_calls"))
    .allele_call_feature_matrix(x, loci)
  else
    .numeric_genetic_feature_matrix(x)
  imputed_allele_calls <- attr(features, "imputed_allele_calls")
  imputed_modal_alleles <- attr(features, "imputed_modal_alleles")

  feature_center <- if (isTRUE(center))
    colMeans(features)
  else
    setNames(rep(0, ncol(features)), colnames(features))
  features <- sweep(features, 2L, feature_center, "-")

  feature_scale <- setNames(rep(1, ncol(features)), colnames(features))
  retained <- rep(TRUE, ncol(features))
  if (isTRUE(scale))
  {
    feature_scale <- apply(features, 2L, stats::sd)
    retained <- is.finite(feature_scale) & feature_scale > tol
    if (!any(retained))
      stop("No variable genetic features remain after scaling", call. = FALSE)
    if (any(!retained))
      warning("Dropping constant or near-constant genetic features: ",
              paste(colnames(features)[!retained], collapse = ", "),
              call. = FALSE)
    features <- sweep(features[, retained, drop = FALSE],
                      2L, feature_scale[retained], "/")
    feature_center <- feature_center[retained]
    feature_scale <- feature_scale[retained]
  }

  normalizer <- if (identical(normalize, "features")) ncol(features) else 1

  if (is.null(groups))
  {
    unit_features <- features
    unit_names <- rownames(features)
    if (is.null(unit_names) || anyNA(unit_names) || any(!nzchar(unit_names)))
      unit_names <- paste0("sample", seq_len(nrow(features)))
    else if (anyDuplicated(unit_names))
      unit_names <- make.unique(unit_names)
    rownames(unit_features) <- unit_names
    unit_size <- setNames(rep(1L, nrow(unit_features)), unit_names)
    within_variance <- NULL
    level <- "individual"
    if (identical(diagonal, "within"))
      stop("`diagonal = \"within\"` requires population groups; use `diagonal = \"gower\"` for individual-level covariance",
           call. = FALSE)
    diagonal <- "gower"
  }
  else
  {
    groups <- if (is.factor(groups))
      droplevels(groups)
    else
      factor(groups, levels = unique(groups))
    if (nlevels(groups) < 2L)
      stop("`groups` must contain at least two populations", call. = FALSE)
    group_size <- as.numeric(table(groups))
    names(group_size) <- levels(groups)

    unit_features <- rowsum(features, groups, reorder = TRUE)
    unit_features <- sweep(unit_features, 1L, group_size, "/")
    rownames(unit_features) <- levels(groups)
    unit_size <- group_size
    level <- if (all(group_size == 1L)) "individual" else "population"

    if (identical(diagonal, "auto"))
      diagonal <- if (identical(level, "individual")) "gower" else "within"
    if (identical(diagonal, "within"))
    {
      if (all(group_size < 2L))
        stop("`diagonal = \"within\"` requires at least one population with two or more samples",
             call. = FALSE)
      if (any(group_size < 2L))
        warning("Some populations contain fewer than two individuals; their within-population variance is set to zero.",
                call. = FALSE)
    }

    within_variance <- vapply(levels(groups), function(g)
    {
      xg <- features[groups == g, , drop = FALSE]
      if (nrow(xg) < 2L)
        return(0)
      sum(apply(xg, 2L, stats::var))
    }, numeric(1))
  }

  centroid_dist2 <- as.matrix(stats::dist(unit_features))^2 / normalizer
  row_mean <- rowMeans(centroid_dist2)
  col_mean <- colMeans(centroid_dist2)
  grand_mean <- mean(centroid_dist2)
  covariance <- -0.5 * (
    centroid_dist2 -
      outer(row_mean, rep(1, length(row_mean))) -
      outer(rep(1, length(col_mean)), col_mean) +
      grand_mean
  )

  if (!is.null(within_variance))
    within_variance <- within_variance / normalizer
  if (identical(diagonal, "within"))
    diag(covariance) <- within_variance

  covariance <- (covariance + t(covariance)) / 2
  rownames(covariance) <- colnames(covariance) <- rownames(unit_features)
  attr(covariance, "centroids") <- unit_features
  attr(covariance, "unit_features") <- unit_features
  if (!is.null(within_variance))
    attr(covariance, "within_variance") <- within_variance
  attr(covariance, "centroid_distance2") <- centroid_dist2
  attr(covariance, "unit_distance2") <- centroid_dist2
  attr(covariance, "unit_size") <- unit_size
  attr(covariance, "feature_center") <- feature_center
  attr(covariance, "feature_scale") <- feature_scale
  attr(covariance, "retained_features") <- colnames(features)
  attr(covariance, "input") <- input
  if (!is.null(imputed_allele_calls))
  {
    attr(covariance, "imputed_allele_calls") <- imputed_allele_calls
    attr(covariance, "imputed_modal_alleles") <- imputed_modal_alleles
  }
  attr(covariance, "level") <- level
  attr(covariance, "diagonal") <- diagonal
  attr(covariance, "centered") <- if (identical(diagonal, "gower"))
    "sites" else "sites_before_diagonal_replacement"
  attr(covariance, "normalize") <- normalize
  attr(covariance, "normalizer") <- normalizer

  covariance
}

.numeric_genetic_feature_matrix <- function(x)
{
  if (!is.numeric(x))
    stop("`x` must be numeric when `input = \"features\"`", call. = FALSE)
  storage.mode(x) <- "double"
  if (anyNA(x))
    stop("missing values are not currently supported in `x`", call. = FALSE)
  if (is.null(colnames(x)))
    colnames(x) <- paste0("feature", seq_len(ncol(x)))
  x
}

.allele_call_feature_matrix <- function(x, loci)
{
  if (is.null(loci))
    stop("`loci` is required when `input = \"allele_calls\"`", call. = FALSE)
  if (length(loci) != ncol(x))
    stop("`loci` must have one entry per allele-call column", call. = FALSE)
  if (anyNA(loci) || any(!nzchar(as.character(loci))))
    stop("`loci` values must be non-missing and non-empty", call. = FALSE)

  x_chr <- matrix(as.character(x), nrow = nrow(x), ncol = ncol(x))
  loci <- as.character(loci)
  locus_names <- unique(loci)
  out <- vector("list", length(locus_names))
  names(out) <- locus_names
  imputed <- integer(length(locus_names))
  names(imputed) <- locus_names
  imputed_allele <- setNames(rep(NA_character_, length(locus_names)),
                             locus_names)

  for (loc in locus_names)
  {
    loc_values <- x_chr[, loci == loc, drop = FALSE]
    missing <- is.na(loc_values)
    if (any(missing))
    {
      observed <- loc_values[!missing]
      if (!length(observed))
        stop("Cannot impute allele calls for locus `", loc,
             "` because all calls are missing", call. = FALSE)
      modal <- .modal_allele(observed)
      loc_values[missing] <- modal
      imputed[[loc]] <- sum(missing)
      imputed_allele[[loc]] <- modal
    }
    alleles <- sort(unique(c(loc_values)))
    loc_counts <- vapply(alleles, function(allele)
      rowSums(loc_values == allele),
      numeric(nrow(loc_values)))
    loc_counts <- as.matrix(loc_counts)
    colnames(loc_counts) <- paste0(loc, ":", alleles)
    out[[loc]] <- loc_counts
  }

  if (sum(imputed) > 0L)
    message("Imputed ", sum(imputed),
            " missing allele call(s) to the modal allele within each locus.")

  features <- do.call(cbind, out)
  attr(features, "imputed_allele_calls") <- imputed[imputed > 0L]
  attr(features, "imputed_modal_alleles") <- imputed_allele[imputed > 0L]
  features
}

.modal_allele <- function(x)
{
  counts <- table(x)
  modes <- names(counts)[counts == max(counts)]
  sort(modes)[1]
}

#' Squared-distance matrix from a covariance matrix
#'
#' Converts a covariance to the squared Euclidean distances it implies. For
#' \code{\link{cov_from_genetic_data}} output with \code{diagonal = "gower"},
#' the result equals its \code{centroid_distance2} attribute; with
#' \code{"within"}, it does not.
#'
#' @param Cov A square numeric covariance matrix.  Does not need to be
#'   full-rank (e.g. the output of \code{\link{cov_from_biallelic}} is
#'   positive semi-definite).
#'
#' @details
#' The squared distance between units \eqn{i} and \eqn{j} is:
#'
#' \deqn{D_{ij} = C_{ii} + C_{jj} - 2 C_{ij}}
#'
#' which is the law-of-cosines identity relating a covariance matrix to its
#' implied squared Euclidean distances.  The diagonal of the returned matrix
#' is zero.
#'
#' Squared distances from a coherent covariance suit
#' \code{terradish::generalized_wishart()}, which gives the same likelihood as
#' \code{terradish::wishart_covariance()} and needs the same \code{nu}.
#' \code{terradish::mlpe()} accepts any pairwise distance and needs no
#' \code{nu}.
#'
#' The result keeps the \code{diagonal} and \code{centered} attributes of
#' \code{Cov}, so terradish can still warn about a \code{"within"} diagonal.
#' Other attributes are dropped.
#'
#' @return A symmetric matrix with the dimnames of \code{Cov} and a zero
#'   diagonal. Off-diagonal entries are nonnegative when \code{Cov} is
#'   positive semidefinite.
#'
#' @seealso \code{\link{cov_from_biallelic}}, \code{\link{cov_from_genetic_data}},
#'   \code{\link{dist_from_biallelic}}
#'
#' @examples
#' Cov <- matrix(c(2.0, 1.2, 0.8,
#'                 1.2, 1.5, 0.7,
#'                 0.8, 0.7, 1.1), nrow = 3, byrow = TRUE)
#' dist_from_cov(Cov)
#'
#' # Round-trip: covariance -> distance -> (verify non-negative, zero diag)
#' Y <- matrix(c(2, 1, 0, 1, 1, 1, 0, 1, 2), nrow = 3, byrow = TRUE)
#' S_cov  <- cov_from_biallelic(Y)
#' S_dist <- dist_from_cov(S_cov)
#' diag(S_dist)   # all zeros
#'
#' @export

dist_from_cov <- function(Cov)
{
  Cov <- as.matrix(Cov)
  if (!is.numeric(Cov) || nrow(Cov) != ncol(Cov))
    stop("`Cov` must be a square numeric matrix", call. = FALSE)
  d <- diag(Cov)
  D <- outer(d, d, "+") - 2 * Cov[, , drop = FALSE]  # drops other attributes
  dimnames(D) <- dimnames(Cov)
  # terradish warns when attr(S, "diagonal") is "within"; keep that label.
  attr(D, "diagonal") <- attr(Cov, "diagonal", exact = TRUE)
  attr(D, "centered") <- attr(Cov, "centered", exact = TRUE)
  D
}

#' Genetic distance matrix from biallelic allele counts
#'
#' Computes pairwise squared genetic distances from biallelic (SNP) allele
#' counts.  This is a convenience wrapper for
#' \code{\link{dist_from_cov}(\link{cov_from_biallelic}(Y, N))}.
#'
#' @param Y Numeric matrix of derived-allele counts (populations/individuals in
#'   rows, loci in columns).
#' @param N Optional haploid sample sizes in any form accepted by
#'   \code{\link{cov_from_biallelic}}. The default \code{NULL} uses
#'   \code{ploidy} (2 unless set through \code{...}) for every cell.
#' @param ... Further arguments passed to \code{\link{cov_from_biallelic}},
#'   such as \code{ploidy}, \code{monomorphic}, or \code{tol}.
#'
#' @return A symmetric matrix of pairwise squared distances with a zero
#'   diagonal, on the count scale of \code{cov_from_biallelic()} and named by
#'   \code{rownames(Y)}. With equal sample sizes, the terradish response
#'   \code{dist_from_cov(cov_from_biallelic(Y, N) / N[1])} equals this matrix
#'   divided by \code{N[1]}. Suitable for \code{terradish::mlpe()} or
#'   \code{terradish::generalized_wishart()}.
#'
#' @seealso \code{\link{cov_from_biallelic}}, \code{\link{dist_from_cov}},
#'   \code{\link{fst_from_biallelic}}
#'
#' @examples
#' Y <- matrix(c(2, 1, 0,
#'               1, 1, 1,
#'               0, 1, 2), nrow = 3, byrow = TRUE)
#' N <- matrix(2, nrow = 3, ncol = 3)
#' dist_from_biallelic(Y, N)
#' dist_from_biallelic(Y)   # N = NULL uses ploidy = 2 for every cell
#'
#' @export

dist_from_biallelic <- function(Y, N = NULL, ...)
{
  dist_from_cov(cov_from_biallelic(Y, N, ...))
}
