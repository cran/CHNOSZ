context("thermo")

# clear out any previous basis definition or database alterations
suppressPackageStartupMessages(data(thermo))

test_that("NAs in thermo$obigt propagate to subcrt()", {
  # first of all, water is in thermo$obigt but its properties
  # are actually calculated using water() so it has NAs for some parameters
  expect_equal(info(1)$a, as.numeric(NA))
  # get the existing value of c for [Ala](cr) (it's 0)
  expect_equal(c.Ala <- info(info("[Ala]", "cr"))$c, 0)
  # when we make a protein, its G depends on temperature
  expect_true(all(diff(subcrt("LYSC_CHICK", "cr")$out[[1]]$G) < 0))
  # turn the values of G and S for [Ala](cr) into NA
  mod.obigt(name="[Ala]", state="cr", G=NA, S=NA)
  # now when we make a protein(cr), its G is NA
  expect_true(all(is.na(subcrt("RNAS1_BOVIN", "cr")$out[[1]]$G)))
  # also check propagation of NA for aqueous species
  mod.obigt(name="[Ala]", state="aq", G=NA, S=NA)
  expect_true(all(is.na(subcrt("[Ala]", "aq")$out[[1]]$G)))
  # be nice and restore the database
  suppressPackageStartupMessages(data(thermo))
})

test_that("minimal usage of mod.obigt() creates usable data entries", {
  # we need at least a name and some property
  expect_error(mod.obigt(list(name="test")))
  # the default state is aq
  expect_message(itest <- mod.obigt(list(name="test", date=today())), "added test\\(aq\\)")
  # we should get NA values of G for a species with NA properties 
  expect_true(all(is.na(subcrt(itest)$out[[1]]$G)))
  # values for Cp and c1 integrate to the same values of G
  G.Cp <- subcrt(mod.obigt(list(name="test", formula="C0", G=0, S=0, Cp=100)))$out[[1]]$G
  G.c1 <- subcrt(mod.obigt(list(name="test", formula="C0", G=0, S=0, c1=100)))$out[[1]]$G
  expect_equal(G.Cp, G.c1)
})
