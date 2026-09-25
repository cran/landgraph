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
  neighbour, or lattice (rook/queen) adjacency. The result is class
  `c("landgraph", "terradish_graph")`.
* Genetic covariance / distance: `cov_from_biallelic()`, `cov_from_genetic_data()`,
  `fst_from_biallelic()`, `dist_from_cov()`, `dist_from_biallelic()`.
* Directional edge covariates: `edge_gradient()` (gradient of a scalar potential)
  and `edge_flow()` (projection of a vector flow field, carrying a curl component).
* `cov_from_biallelic()` applies the per-locus pooled frequency correctly per column
  (fixes a column-major recycling error that affected multi-deme covariances).
