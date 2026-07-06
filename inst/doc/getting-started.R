## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment  = "#>",
  fig.width  = 6,
  fig.height = 4.25
)

## ----load---------------------------------------------------------------------
library(landgraph)

## ----graph--------------------------------------------------------------------
coords <- as.matrix(expand.grid(x = 0:2, y = 0:1))   # 6 demes on a 3 x 2 grid
rownames(coords) <- paste0("deme", seq_len(nrow(coords)))
coords

g <- deme_graph(coords, neighbours = "lattice")
g

## ----graph-parts--------------------------------------------------------------
str(g)
g$edge_pairs

## ----graph-plot, fig.alt = "Six demes on a grid joined by rook-adjacency edges."----
op <- par(no.readonly = TRUE)
par(mar = c(4, 4, 1, 1))
plot(g$vertex_coordinates, type = "n", asp = 1, xlab = "x", ylab = "y")
ep <- g$edge_pairs
segments(g$vertex_coordinates[ep[, 1], 1], g$vertex_coordinates[ep[, 1], 2],
         g$vertex_coordinates[ep[, 2], 1], g$vertex_coordinates[ep[, 2], 2],
         col = "grey60")
points(g$vertex_coordinates, pch = 21, bg = "steelblue", cex = 3)
text(g$vertex_coordinates, labels = seq_len(nrow(coords)), col = "white")
par(op)

## ----graph-knn----------------------------------------------------------------
g_knn <- deme_graph(coords, neighbours = "knn", k = 2)
nrow(g_knn$edge_pairs)   # edge count under 2-nearest-neighbour adjacency

## ----snp-data-----------------------------------------------------------------
set.seed(42)
n_demes <- nrow(coords)
n_snp   <- 8
freqs   <- runif(n_snp, 0.1, 0.9)                 # a true frequency per SNP
Y <- vapply(freqs, function(p) rbinom(n_demes, size = 40, prob = p),
            numeric(n_demes))
rownames(Y) <- rownames(coords)
colnames(Y) <- paste0("snp", seq_len(n_snp))
Y

## ----snp-cov------------------------------------------------------------------
S <- cov_from_biallelic(Y, N = 40)
round(S, 3)

## ----fst----------------------------------------------------------------------
Nmat <- matrix(40, n_demes, n_snp)
fst <- fst_from_biallelic(Y, Nmat)
round(fst, 4)

## ----multi-features-----------------------------------------------------------
x <- matrix(c(0, 1,
              1, 1,
              2, 0,
              2, 1,
              0, 2,
              1, 2),
            ncol = 2, byrow = TRUE)
groups <- rep(c("pop1", "pop2", "pop3"), each = 2)

Sg <- cov_from_genetic_data(x, groups = groups)
round(Sg, 3)

## ----multi-attrs--------------------------------------------------------------
attr(Sg, "level")            # "population" once groups have replication
attr(Sg, "diagonal")         # the diagonal rule actually applied
attr(Sg, "within_variance")  # the within-population variances on the diagonal
attr(Sg, "unit_size")        # number of individuals per population

## ----multi-msat---------------------------------------------------------------
alleles <- data.frame(
  loc1_a = c(100, 100, 102, 102, 104, 104),
  loc1_b = c(100, 102, 102, 104, 104, 100),
  loc2_a = c(200, 202, 200, 202, 204, 204),
  loc2_b = c(202, 202, 204, 204, 204, 200)
)

Sm <- cov_from_genetic_data(
  alleles,
  groups = groups,
  input  = "allele_calls",
  loci   = c("loc1", "loc1", "loc2", "loc2")
)
round(Sm, 3)

## ----dist---------------------------------------------------------------------
D <- dist_from_cov(S)
round(D, 3)

## ----dist-wrap----------------------------------------------------------------
D2 <- dist_from_biallelic(Y, N = 40)
all.equal(D, D2)   # same as the two-step route above

## ----grad---------------------------------------------------------------------
elevation <- coords[, "x"] + coords[, "y"]   # one value per deme
elevation

eg <- edge_gradient(elevation, g)
str(eg)

## ----flow---------------------------------------------------------------------
cen <- colMeans(g$vertex_coordinates)
rotational <- function(xy) cbind(-(xy[, 2] - cen[2]), xy[, 1] - cen[1])

circ <- edge_flow(rotational, g)
round(circ, 3)

## ----flow-plot, fig.alt = "A counter-clockwise vector field drawn as arrows at each deme over the graph edges."----
op <- par(no.readonly = TRUE)
par(mar = c(4, 4, 1, 1))
vc  <- g$vertex_coordinates
fld <- rotational(vc)
plot(vc, type = "n", asp = 1, xlab = "x", ylab = "y",
     xlim = range(vc[, 1]) + c(-0.4, 0.4),
     ylim = range(vc[, 2]) + c(-0.4, 0.4))
segments(vc[ep[, 1], 1], vc[ep[, 1], 2],
         vc[ep[, 2], 1], vc[ep[, 2], 2], col = "grey80")
arrows(vc[, 1], vc[, 2],
       vc[, 1] + 0.3 * fld[, 1], vc[, 2] + 0.3 * fld[, 2],
       length = 0.08, col = "firebrick")
points(vc, pch = 21, bg = "steelblue", cex = 2.5)
par(op)

## ----quick-ref, eval = FALSE--------------------------------------------------
# library(landgraph)
# 
# # 1. Build the spatial graph from deme coordinates
# g <- deme_graph(coords, neighbours = "lattice")    # or "knn" / "delaunay"
# 
# # 2. Genetic covariance and distance from biallelic (SNP) counts
# S <- cov_from_biallelic(Y, N = 40)                 # deme covariance
# D <- dist_from_cov(S)                              # squared genetic distance
# D <- dist_from_biallelic(Y, N = 40)                # the two steps in one call
# fst <- fst_from_biallelic(Y, Nmat)                 # pairwise F_ST (N as a matrix)
# 
# # 2b. Covariance from microsatellite / multivariate data
# Sg <- cov_from_genetic_data(x, groups = groups)    # population covariance
# Sg <- cov_from_genetic_data(alleles, groups = groups,
#                             input = "allele_calls",
#                             loci = c("loc1", "loc1", "loc2", "loc2"))
# 
# # 3. Directional edge covariates for directed models
# eg   <- edge_gradient(elevation, g)                # downhill drop; eg$d per directed edge
# circ <- edge_flow(rotational, g)                   # circulation per undirected edge
# 
# # 4. Hand off downstream
# #    terradish::wishart_covariance(S = S, ...)      # symmetric resistance
# #    dragonflow::dragon(..., circulation = circ)    # asymmetric gene flow

