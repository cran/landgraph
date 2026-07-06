test_that("edge_gradient builds an antisymmetric node-potential gradient", {
  g <- deme_graph(as.matrix(expand.grid(x = 0:2, y = 0:2)), neighbours = "lattice")
  elev <- g$vertex_coordinates[, 1]                 # potential = x coordinate
  eg <- edge_gradient(elev, g)
  expect_named(eg, c("edges", "d"))
  # directed edge a->b carries elev_a - elev_b; its reverse carries the negative
  m <- nrow(g$edge_pairs)
  expect_equal(eg$d[seq_len(m)], -eg$d[m + seq_len(m)])
})

test_that("edge_flow builds an antisymmetric circulation covariate from a vector field", {
  g <- deme_graph(as.matrix(expand.grid(x = 0:2, y = 0:2)), neighbours = "lattice")
  vc <- g$vertex_coordinates
  cen <- colMeans(vc)
  rot <- function(xy) cbind(-(xy[, 2] - cen[2]), xy[, 1] - cen[1])   # counter-clockwise curl

  circ <- edge_flow(rot, g)
  expect_length(circ, nrow(g$edge_pairs))
  expect_true(any(abs(circ) > 1e-8))                                 # non-trivial flow
  expect_equal(edge_flow(rot(vc), g), circ)                          # function == matrix form

  # reversing each edge's orientation negates the covariate (antisymmetry)
  g_rev <- g; g_rev$edge_pairs <- g$edge_pairs[, 2:1, drop = FALSE]
  expect_equal(edge_flow(rot, g_rev), -circ)

  expect_error(edge_flow(matrix(1, nrow(vc), 3), g), "two-column")
})
