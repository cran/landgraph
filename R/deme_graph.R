#' Lightweight deme / landscape graph from coordinates
#'
#' Builds the minimal spatial graph that DRAGON-style network methods consume: a
#' set of vertex coordinates and an undirected edge list. For deme-scale problems
#' (dozens to a few hundred nodes) this avoids the heavier raster / conductance
#' graph builders, and the returned object is drop-in compatible with the graph
#' consumed by \code{dragon} (in dragonflow) and with a \code{terradish_graph}.
#'
#' @param coords A two-column numeric matrix or data frame of node (deme)
#'   coordinates, one row per node (columns x, y).
#' @param neighbours Edge construction. \code{"delaunay"} (the Delaunay
#'   triangulation; needs the \pkg{deldir} package), \code{"knn"} (symmetric
#'   k-nearest-neighbour adjacency), or \code{"lattice"} (rook or queen adjacency
#'   for points lying on an integer grid).
#' @param k Number of neighbours for \code{neighbours = "knn"}.
#' @param queen For \code{neighbours = "lattice"}, also connect diagonal (queen)
#'   neighbours; the default rook adjacency reproduces
#'   \code{conductance_surface(directions = 4)}.
#' @return An object of class \code{c("landgraph", "terradish_graph")} with
#'   \code{vertex_coordinates} (an \code{n x 2} matrix), \code{edge_pairs} (an
#'   \code{m x 2} integer matrix of 1-based undirected edges, each pair once,
#'   \code{a < b}), and \code{n_vertices}.
#' @seealso \code{\link{edge_gradient}}, \code{\link{edge_flow}}
#' @examples
#' coords <- as.matrix(expand.grid(x = 0:4, y = 0:4))
#' g <- deme_graph(coords, neighbours = "lattice")
#' nrow(g$edge_pairs)
#' @export
deme_graph <- function(coords, neighbours = c("delaunay", "knn", "lattice"),
                       k = 6L, queen = FALSE)
{
  neighbours <- match.arg(neighbours)
  coords <- as.matrix(coords)
  if (ncol(coords) != 2L)
    stop("`coords` must have two columns (x, y).", call. = FALSE)
  storage.mode(coords) <- "double"
  if (nrow(coords) < 2L)
    stop("`coords` needs at least two nodes.", call. = FALSE)

  ep <- switch(neighbours,
    delaunay = .deme_edges_delaunay(coords),
    knn      = .deme_edges_knn(coords, k),
    lattice  = .deme_edges_lattice(coords, queen))

  ep <- matrix(as.integer(ep), ncol = 2L)
  ep <- ep[ep[, 1] != ep[, 2], , drop = FALSE]
  ep <- t(apply(ep, 1L, function(e) if (e[1] < e[2]) e else e[2:1]))  # a < b
  ep <- ep[!duplicated(ep), , drop = FALSE]
  ep <- ep[order(ep[, 1], ep[, 2]), , drop = FALSE]
  storage.mode(ep) <- "integer"

  structure(list(vertex_coordinates = coords, edge_pairs = ep,
                 coords = coords, n_vertices = nrow(coords)),
            class = c("landgraph", "terradish_graph"))
}

#' @export
print.landgraph <- function(x, ...)
{
  cat("landgraph: ", nrow(x$vertex_coordinates), " vertices, ",
      nrow(x$edge_pairs), " undirected edges\n", sep = "")
  invisible(x)
}

# rook (|dx|+|dy| == 1) or queen (also the four diagonals) adjacency for points
# on an integer grid; distance-based so spacing need only be uniform.
.deme_edges_lattice <- function(coords, queen)
{
  d <- as.matrix(stats::dist(coords))
  step <- min(d[d > 1e-9])
  thr <- if (queen) sqrt(2) * step + 1e-6 else step + 1e-6
  idx <- which(d > 1e-9 & d <= thr, arr.ind = TRUE)
  idx[idx[, 1] < idx[, 2], , drop = FALSE]
}

.deme_edges_knn <- function(coords, k)
{
  n <- nrow(coords)
  d <- as.matrix(stats::dist(coords))
  diag(d) <- Inf
  k <- min(as.integer(k), n - 1L)
  do.call(rbind, lapply(seq_len(n), function(i) cbind(i, order(d[i, ])[seq_len(k)])))
}

.deme_edges_delaunay <- function(coords)
{
  if (!requireNamespace("deldir", quietly = TRUE))
    stop("neighbours = \"delaunay\" needs the 'deldir' package; ",
         "install it, or use \"knn\" / \"lattice\".", call. = FALSE)
  dd <- deldir::deldir(coords[, 1], coords[, 2], suppressMsge = TRUE)
  as.matrix(dd$delsgs[, c("ind1", "ind2")])
}
