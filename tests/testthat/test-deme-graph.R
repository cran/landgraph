test_that("deme_graph builds a valid graph with rook/queen/knn adjacency", {
  coords <- as.matrix(expand.grid(x = 0:2, y = 0:2))   # 3 x 3 integer grid

  g <- deme_graph(coords, neighbours = "lattice")
  expect_s3_class(g, "landgraph")
  expect_true(inherits(g, "terradish_graph"))
  expect_equal(nrow(g$vertex_coordinates), 9L)
  expect_equal(nrow(g$edge_pairs), 12L)                # rook edges of a 3x3 grid
  expect_true(all(g$edge_pairs[, 1] < g$edge_pairs[, 2]))
  expect_false(any(duplicated(g$edge_pairs)))

  gq <- deme_graph(coords, neighbours = "lattice", queen = TRUE)
  expect_equal(nrow(gq$edge_pairs), 20L)               # 12 rook + 8 diagonal

  gk <- deme_graph(coords, neighbours = "knn", k = 2)
  expect_true(nrow(gk$edge_pairs) > 0L)
  expect_true(all(gk$edge_pairs[, 1] < gk$edge_pairs[, 2]))

  expect_error(deme_graph(matrix(1:9, ncol = 3)), "two columns")
})
