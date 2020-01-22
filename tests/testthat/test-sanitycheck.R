test_that("Not a matrix", {
  df <- data.frame(a = 1:5, b = 1:5)
  expect_error(sanity_check(df, TRUE))
})

test_that("Empty rows", {
  mat <- matrix(1,5,5)
  row.names(mat) <- 1:5
  mat[2,] <- 0
  expect_error(sanity_check(mat, TRUE))
})


test_that("Empty names", {
  mat <- matrix(1,5,5)
  expect_warning(sanity_check(mat, TRUE))

})
