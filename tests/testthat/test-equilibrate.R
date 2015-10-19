context("equilibrate")

# set up some simple systems
# one acid
basis("CHNOS+")
species("acetic acid")
aone <- suppressMessages(affinity())
# acids
species(c("formic acid", "formate", "acetate"))
aacid <- suppressMessages(affinity())
# acids plus a sulfur species
species("H2S")
aacidS <- suppressMessages(affinity())
# proteins
species(delete=TRUE)
species(c("LYSC_CHICK", "MYG_PHYCA", "RNAS1_BOVIN", "CYC_BOVIN"))
aprot <- suppressMessages(affinity())

test_that("equilibrate() gives expected messages and errors for balance calculation", {
  # the following error is triggered by equil.react, not equil.boltzmann
  expect_error(equilibrate(aone), "at least two species needed")
  expect_message(equilibrate(aacid), "coefficients are moles of CO2")
  expect_message(equilibrate(aacid), "balancing coefficients are 2 1 1 2 ")
  expect_message(equilibrate(aacid), "logarithm of total moles of CO2 is -2.221848")
  expect_message(equilibrate(aacid, loga.balance=-3), "logarithm of total moles of CO2 \\(from loga.balance\\) is -3")
  expect_error(equilibrate(aacid, balance="length"), "some species are not proteins")
  expect_message(equilibrate(aacidS), "coefficients are unity")
  expect_message(equilibrate(aacidS), "balancing coefficients are 1 1 1 1 1")
  expect_message(equilibrate(aacidS), "logarithm of total moles of species is -2.301029")
  expect_error(equilibrate(aacidS, balance="CO2"), "some species have no CO2 in the formation reaction")
  expect_message(equilibrate(aprot), "coefficients are protein length")
  expect_message(equilibrate(aprot), "balancing coefficients are 129 153 124 104")
  expect_message(equilibrate(aprot), "logarithm of total protein length is -0.292429")
  expect_message(equilibrate(aprot, normalize=TRUE), "using 'normalize' for molar formulas")
})

test_that("equilibrate() gives expected messages and errors for species selection", {
  # an error if we select no species
  expect_error(equilibrate(aacid, ispecies=numeric()), "the length of ispecies is zero")
  # an error if all affinities are NA
  aNA <- aacid
  aNA$values[1:2] <- NA
  expect_error(equilibrate(aNA, ispecies=1:2), "all species have NA affinities")
  # a message if we select only certain of the species
  expect_message(equilibrate(aacid, ispecies=1:2), "using 2 of 4 species")
})

test_that("equilibrate() keeps the same total loga.balance for normalize=TRUE or FALSE", {
  # use the proteins
  e.norm <- equilibrate(aprot, normalize=TRUE)
  e <- equilibrate(aprot)
  # the total activity of the balance in the two cases
  sumact.balance.norm <- sum(10^unlist(e.norm$loga.equil)*e.norm$m.balance)
  sumact.balance <- sum(10^unlist(e$loga.equil)*e$n.balance)
  expect_equal(sumact.balance.norm, sumact.balance)
})

test_that("equilibrate() reproduces an example from the literature", {
  # the reference values are the equilibrium logarithms of activities
  # of sulfur species at logfO2=-30 from Seewald, 2001
  # we name them here because S5O6-2 isn't on the plot at logfO2=-30, 
  # and to get them in order
  species.ref <- c("S3O6-2", "S2O6-2", "S2O4-2", "S3-2", "S2-2", "S2O3-2", "HSO3-", "SO2", "HSO4-", "H2S")
  # these values were read from the plot using g3data (http://www.frantz.fi/software/g3data.php)
  loga.ref <- c(-28.82, -24.70, -22.10, -14.19, -12.12, -11.86, -8.40, -7.40, -6.54, -1.95)
  # set up the system - see ?diagram for an example showing the entire plot
  basis("CHNOS+")
  basis(c("pH", "O2"), c(5, -30))
  # we include here all the species shown by Seewald, 2001
  species(c("H2S", "S2-2", "S3-2", "S2O3-2", "S2O4-2", "S3O6-2", "S5O6-2", "S2O6-2", "HSO3-", "SO2", "HSO4-"))
  a <- affinity(T=325, P=350)
  # loga.balance=-2 signifies 10 mmolal total sulfur
  e <- equilibrate(a, loga.balance=-2)
  # get the calculated activities of the reference species
  loga.equil <- unlist(e$loga.equil[match(species.ref, e$species$name)])
  # the test... the tolerance may seem high, but consider that the reference values
  # were read from a plot with 30 logfO2 units spanning 4 inches
  expect_true(all(abs(loga.equil-loga.ref) < 0.36))
})

#}

# references

# Seewald, J. S. (2001) 
#   Aqueous geochemistry of low molecular weight hydrocarbons at elevated temperatures and
#   pressures: Constraints from mineral buffered laboratory experiments
#   Geochim. Cosmochim. Acta 65, 1641--1664. http://dx.doi.org/10.1016/S0016-7037(01)00544-0
