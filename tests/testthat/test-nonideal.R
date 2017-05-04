context("nonideal")

test_that("loggam and logK values are consistent", {
  rxn1 <- subcrt(c("anhydrite", "Ca+2", "SO4-2"), c(-1, 1, 1), P=1, T=25, I=0)
  rxn2 <- subcrt(c("anhydrite", "Ca+2", "SO4-2"), c(-1, 1, 1), P=1, T=25, I=0.7)
  expect_equal(rxn1$out$logK, rxn2$out$loggam + rxn2$out$logK)
})
