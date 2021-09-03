test_that("SingVarRiskSummaries works", {
  # setup
  set.seed(111)
  dat <- SimData(n = 50, M = 4)
  y <- dat$y
  Z <- dat$Z
  X <- dat$X
  
  set.seed(111)
  fitkm <- kmbayes(y = y, Z = Z, X = X, verbose = FALSE)
  
  expect_snapshot_value(
    SingVarRiskSummaries(fitkm),
    style = "serialize"
  )
})
