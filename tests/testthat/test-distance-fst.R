test_that("dist_from_cov keeps names and inverts Gower", {
  x <- matrix(c(0, 1, 1, 1, 2, 0, 2, 1, 0, 2, 1, 2), ncol = 2, byrow = TRUE)
  C <- cov_from_genetic_data(x, groups = rep(c("p1", "p2", "p3"), each = 2))
  D <- dist_from_cov(C)
  expect_identical(dimnames(D), dimnames(C))
  expect_equal(D[, ], attr(C, "centroid_distance2"))
  Y <- matrix(c(2, 1, 0, 1, 1, 1, 0, 1, 2), 3, byrow = TRUE)
  expect_equal(dist_from_biallelic(Y), dist_from_cov(cov_from_biallelic(Y)))
})

test_that("dist_from_cov keeps only the construction labels", {
  x <- matrix(c(0, -1, 0, 1, 1, -1, 1, 1, 2, -1, 2, 1), ncol = 2, byrow = TRUE)
  groups <- rep(c("p1", "p2", "p3"), each = 2)
  within <- cov_from_genetic_data(x, groups, scale = FALSE, diagonal = "within")
  D <- dist_from_cov(within)
  # terradish warns on a "within" diagonal, so the label must survive.
  expect_equal(attr(D, "diagonal"), "within")
  expect_equal(attr(D, "centered"), "sites_before_diagonal_replacement")
  expect_null(attr(D, "centroid_distance2"))
  expect_false(isTRUE(all.equal(D[, ], attr(within, "centroid_distance2"))))

  plain <- dist_from_cov(matrix(c(2, 1, 1, 3), 2))
  expect_identical(names(attributes(plain)), "dim")
  expect_equal(diag(plain), c(0, 0))
  expect_error(dist_from_cov(matrix(1:6, 2)), "square numeric")
})

test_that("dist_from_biallelic passes N and further arguments through", {
  Y <- matrix(c(2, 1, 0, 1, 1, 1, 0, 1, 2), 3, byrow = TRUE,
              dimnames = list(c("a", "b", "c"), NULL))
  expect_equal(dist_from_biallelic(Y, 4), dist_from_cov(cov_from_biallelic(Y, 4)))
  expect_equal(dist_from_biallelic(Y, ploidy = 4),
               dist_from_cov(cov_from_biallelic(Y, ploidy = 4)))
  expect_identical(rownames(dist_from_biallelic(Y)), rownames(Y))
})

test_that("fst_from_biallelic matches Hudson and aligns N", {
  set.seed(1)
  Y <- matrix(rbinom(15, 16, 0.4), 3, dimnames = list(c("p1", "p2", "p3"), NULL))
  N <- matrix(20, 3, 5); p <- Y / N
  num <- (p[1, ] - p[2, ])^2 - p[1, ] * (1 - p[1, ]) / 19 - p[2, ] * (1 - p[2, ]) / 19
  den <- p[1, ] * (1 - p[2, ]) + p[2, ] * (1 - p[1, ])
  fst <- fst_from_biallelic(Y, N)
  expect_equal(fst[1, 2], mean(num) / mean(den))
  expect_identical(rownames(fst), rownames(Y))
  expect_equal(fst_from_biallelic(N - Y, N), fst)
  Nl <- c(20, 20, 18, 20, 16)
  expect_equal(fst_from_biallelic(Y, Nl),
               fst_from_biallelic(Y, matrix(Nl, 3, 5, byrow = TRUE)))
  expect_equal(fst_from_biallelic(Y, 20), fst)
  expect_error(fst_from_biallelic(Y, matrix(20, 3, 4)), "Dimension mismatch")
  expect_error(fst_from_biallelic(matrix(c(0, 1, 1, 0), 2), 1), "must exceed 1")
})
