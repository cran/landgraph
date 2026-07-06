# landgraph 0.0.1

* Initial release. Shared landscape-genetic primitives extracted from terradish so
  that the symmetric (terradish) and asymmetric (dragonflow) network methods stand
  on a common, dependency-light base.
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
