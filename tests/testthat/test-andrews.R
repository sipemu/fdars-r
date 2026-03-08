test_that("andrews_transform returns fdata with correct dimensions", {
  X <- matrix(rnorm(50), nrow = 10, ncol = 5)
  fd <- andrews_transform(X)

  expect_s3_class(fd, "fdata")
  expect_equal(nrow(fd$data), 10)
  expect_equal(ncol(fd$data), 200)
  expect_equal(length(fd$argvals), 200)
  expect_equal(range(fd$argvals), c(-pi, pi))
})

test_that("andrews_transform accepts custom m", {
  X <- matrix(rnorm(30), nrow = 5, ncol = 6)
  fd <- andrews_transform(X, m = 50)

  expect_equal(ncol(fd$data), 50)
  expect_equal(length(fd$argvals), 50)
})

test_that("andrews_transform accepts data.frame input", {
  df <- data.frame(a = 1:5, b = 6:10, c = 11:15)
  fd <- andrews_transform(df)

  expect_s3_class(fd, "fdata")
  expect_equal(nrow(fd$data), 5)
})

test_that("andrews_transform attaches basis and varnames attributes", {
  X <- matrix(rnorm(30), nrow = 5, ncol = 6)
  colnames(X) <- paste0("var", 1:6)
  fd <- andrews_transform(X, m = 100)

  basis <- attr(fd, "andrews_basis")
  varnames <- attr(fd, "andrews_varnames")

  expect_true(!is.null(basis))
  expect_equal(dim(basis), c(100, 6))
  expect_equal(varnames, paste0("var", 1:6))
})

test_that("andrews_transform generates default varnames when colnames NULL", {
  X <- matrix(rnorm(20), nrow = 4, ncol = 5)
  fd <- andrews_transform(X)

  expect_equal(attr(fd, "andrews_varnames"), paste0("V", 1:5))
})

test_that("andrews_transform rejects non-numeric input", {
  expect_error(andrews_transform("not a matrix"),
               "numeric matrix or data.frame")
})

test_that("andrews_transform basis has correct structure", {
  X <- matrix(rnorm(30), nrow = 5, ncol = 6)
  fd <- andrews_transform(X, m = 50)
  basis <- attr(fd, "andrews_basis")

  # First column is constant 1/sqrt(2)
  expect_equal(basis[, 1], rep(1 / sqrt(2), 50))

  # Columns alternate sin/cos with increasing frequency
  t_grid <- seq(-pi, pi, length.out = 50)
  expect_equal(basis[, 2], sin(1 * t_grid))    # j=2: sin(t)
  expect_equal(basis[, 3], cos(1 * t_grid))    # j=3: cos(t)
  expect_equal(basis[, 4], sin(2 * t_grid))    # j=4: sin(2t)
  expect_equal(basis[, 5], cos(2 * t_grid))    # j=5: cos(2t)
  expect_equal(basis[, 6], sin(3 * t_grid))    # j=6: sin(3t)
})

test_that("andrews_loadings returns correct structure", {
  X <- scale(iris[, 1:4])
  fd <- andrews_transform(X)
  fpca <- fdata2pc(fd, ncomp = 3)
  loadings <- andrews_loadings(fpca, fd)

  expect_s3_class(loadings, "andrews_loadings")
  expect_s3_class(loadings, "data.frame")
  expect_named(loadings, c("Variable", "Loading", "PC"))
  expect_equal(nrow(loadings), 4 * 3)  # 4 variables x 3 PCs
  expect_equal(unique(loadings$PC), c("PC1", "PC2", "PC3"))
  expect_equal(loadings$Variable[1:4], colnames(iris)[1:4])
})

test_that("andrews_loadings uses correct variable names", {
  X <- matrix(rnorm(30), nrow = 10, ncol = 3)
  colnames(X) <- c("alpha", "beta", "gamma")
  fd <- andrews_transform(X)
  fpca <- fdata2pc(fd, ncomp = 2)
  loadings <- andrews_loadings(fpca, fd)

  expect_equal(unique(loadings$Variable), c("alpha", "beta", "gamma"))
})

test_that("andrews_loadings ncomp parameter works", {
  X <- scale(iris[, 1:4])
  fd <- andrews_transform(X)
  fpca <- fdata2pc(fd, ncomp = 4)

  loadings2 <- andrews_loadings(fpca, fd, ncomp = 2)
  expect_equal(nrow(loadings2), 4 * 2)
  expect_equal(unique(loadings2$PC), c("PC1", "PC2"))

  loadings_all <- andrews_loadings(fpca, fd)
  expect_equal(nrow(loadings_all), 4 * 4)
})

test_that("andrews_loadings errors on missing basis attribute", {
  X <- scale(iris[, 1:4])
  fd <- andrews_transform(X)
  fpca <- fdata2pc(fd, ncomp = 2)

  # Remove the attribute
  fd_plain <- fd
  attr(fd_plain, "andrews_basis") <- NULL

  expect_error(andrews_loadings(fpca, fd_plain),
               "andrews_transform")
})

test_that("andrews_loadings errors on wrong fpca class", {
  X <- scale(iris[, 1:4])
  fd <- andrews_transform(X)

  expect_error(andrews_loadings(list(), fd),
               "fdata2pc")
})
