context("iprotein")

# clear out any prior database alterations
suppressPackageStartupMessages(data(thermo))

test_that("basic searches and conversions work as expected", {
  expect_equal(iprotein(c("LYSC_CHICK", "MYGPHYCA")), c(6, NA))

})

test_that("errors and messages occur in some circumstances", {
  expect_message(iprotein(c("LYSC_CHICK", "MYGPHYCA")), "1 protein not matched")
  expect_error(seq2aa("LYS_CHICK", "XXX"), "no characters match an amino acid")
  expect_error(add.protein(count.aa("AAA")), "not a data frame with the same columns as thermo\\$protein")
  expect_message(add.protein(ip2aa(iprotein("CYC_BOVIN"))), "added 0 of 1 proteins")
})

test_that("group additivity for proteins gives expected values", {
  # values for chicken lysozyme calculated using group additivity values
  # from Dick et al., 2006 (Biogeosciences 3, 311-336)
  G <- -4206050
  Cp <- 6415.5
  V <- 10421
  formula <- "C613H959N193O185S10"
  # use add.obigt to load the parameters for [Met] sidechain group
  # from above reference instead of the updated values
  # from LaRowe and Dick, 2012 (Geochim Cosmochim Acta 80, 70-91)
  add.obigt()
  lprop <- info(info("LYSC_CHICK"))
  expect_equal(G, lprop$G)
  expect_equal(Cp, lprop$Cp, tolerance=1e-5)
  expect_equal(V, lprop$V, tolerance=1e-4)
  expect_equal(formula, lprop$formula)
})

test_that("amino acid counts taken from a fasta file can be added",{
  ffile <- system.file("extdata/fasta/EF-Tu.aln", package="CHNOSZ")
  aa <- read.aa(ffile)
  expect_message(ip1 <- add.protein(aa), "added 8 of 8")
  expect_message(ip2 <- add.protein(aa), "added 0 of 8")
  expect_equal(ip1, ip2)
})

test_that("add.protein returns correct indices for existing names", {
  file <- system.file("data/protein.csv", package="CHNOSZ")
  aa <- read.aa(file)
  aa$protein[1] <- "TEST"
  # only the first one is added, because the others already exist
  expect_message(ip <- add.protein(aa), "added 1 of")
  expect_equal(thermo$protein$protein[ip], aa$protein)
})

# for the future... make info() faster!
# (especially the loading of ionizable groups for proteins)
#test_that("calculations for ionized proteins meet performance expectations", {
#  expect_that({
#    basis("CHNOS+")
#    i <- info(c("LYSC_CHICK","RNAS1_BOVIN"))
#    species(i)
#    a <- affinity()
#  }, takes_less_than(0.4))
#})
