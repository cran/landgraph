# landgraph 0.0.4

* Documentation now matches terradish. `nu` is an effective-information value
  that the user supplies, not a SNP, locus, or feature count. The help for
  `cov_from_biallelic()` explains the count-standardized scale and why
  terradish divides by `N[1]` only with equal sample sizes. The help and the
  vignette describe the `"gower"` default for `diagonal` and how missing data
  and `NA` coding are handled.
* `dist_from_cov()` keeps the dimnames of its input and drops every attribute
  except `diagonal` and `centered`, so terradish can still flag a `"within"`
  diagonal.
* `fst_from_biallelic()` names its result by `rownames(Y)` and accepts `N` in
  any form accepted by `cov_from_biallelic()`. A locus-length `N` vector was
  previously recycled down the columns and misaligned; that case now gives
  different (correct) values. No other results change.
* `dist_from_biallelic()` gains `N = NULL` and passes `...` to
  `cov_from_biallelic()`.
* `deme_graph()` requires `k` to be a single positive whole number for
  `neighbours = "knn"`, and documents its `coords` element.
* `deme_graph()` objects are no longer described as interchangeable with
  terradish fitting graphs, and `edge_gradient()` no longer cites a terradish
  directed model that does not exist.
* DESCRIPTION cites the implemented methods, and `citation("landgraph")`
  reports the installed version.

# landgraph 0.0.3

* Default to the Gower diagonal for grouped genetic covariance. This preserves
  site centering and avoids mixing covariance and within-population variance
  scales. Explicit within and legacy auto choices remain available.

# landgraph 0.0.2

* Record covariance diagonal construction and centering metadata without
  changing numerical values. Distinguish pooled allele-frequency centering
  with unequal sample sizes and centering before diagonal replacement.
* Clarify rare-variant weighting, unequal sampling variance, and the use of
  Gower covariances for Wishart models. FST ratio estimates are not Wishart
  responses.

# landgraph 0.0.1

* Initial release. Shared landscape-genetic primitives extracted from terradish so
  that landscape-genetic network methods stand on a common, dependency-light
  base.
* Graph: `deme_graph()` builds a lightweight deme/landscape graph (vertex
  coordinates + undirected edge list) from coordinates, with Delaunay, k-nearest-
  neighbor, or lattice (rook/queen) adjacency. The result is class
  `c("landgraph", "terradish_graph")`.
* Genetic covariance / distance: `cov_from_biallelic()`, `cov_from_genetic_data()`,
  `fst_from_biallelic()`, `dist_from_cov()`, `dist_from_biallelic()`.
* Directional edge covariates: `edge_gradient()` (gradient of a scalar potential)
  and `edge_flow()` (projection of a vector flow field, carrying a curl component).
* `cov_from_biallelic()` applies the per-locus pooled frequency correctly per column
  (fixes a column-major recycling error that affected multi-deme covariances).
