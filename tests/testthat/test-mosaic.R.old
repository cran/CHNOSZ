context("mosaic")

test_that("results are consistent with affinity()", {
  basis(c("CO2", "H2O", "NH3", "O2"), c(0, 0, 0, 0))
  species(c("alanine", "glycine"))
  a <- affinity()
  # this is a degenerate case because we only allow NH3 to swap for NH3, and CO2 for CO2;
  # however it still exercises the affinity scaling and summing code
  m1 <- mosaic("NH3", "CO2", blend=TRUE)
  # this failed before we divided by loga.tot to get _relative_ abundances of basis species in mosaic.R
  expect_equal(a$values, m1$A.species$values)
  # the next call failed when which.pmax(), called by diagram(), choked on a list of length one
  m2 <- mosaic("NH3", "CO2")
  expect_equal(a$values, m2$A.species$values)
})

# TODO: test that basis specifications can be exchanged between bases and bases2 without altering output
