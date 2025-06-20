## ----HTML, include=FALSE------------------------------------------------------
# Some frequently used HTML expressions
# Use lowercase because some of these are used as variables in the examples
h2o <- "H<sub>2</sub>O"
h2 <- "H<sub>2</sub>"
o2 <- "O<sub>2</sub>"
co2 <- "CO<sub>2</sub>"
h2s <- "H<sub>2</sub>S"
Psat <- "<i>P</i><sub>sat</sub>"
zc <- "<i>Z</i><sub>C</sub>"

## ----setup, include=FALSE---------------------------------------------------------------------------------------------------------------------------------------------------------
library(knitr)

options(width = 180)

# Function to colorize messages 20171031
# Adapted from https://gist.github.com/yihui/2629886#file-knitr-color-msg-rnw
color_block = function(color) {
  function(x, options) sprintf('<pre style="color:%s">%s</pre>', color, x)
}
# Set knitr hooks
knit_hooks$set(
  warning = color_block('magenta'),
  error = color_block('red'),
  message = color_block('blue'),
  # Use pngquant to optimize PNG images
  pngquant = hook_pngquant
)

# Define a custom hook to use smaller plot margins
knit_hooks$set(small_mar = function(before, options, envir) {
  if(before) {
    par(mar = c(5, 5, 1, 1))
  }
})

# Define variables used below
basedpi <- if(nzchar(Sys.getenv("CHNOSZ_BUILD_LARGE_VIGNETTES"))) 60 else 40
pngquant <- if(nzchar(Sys.getenv("CHNOSZ_BUILD_LARGE_VIGNETTES"))) "--speed=1 --quality=0-25" else "--speed=1 --quality=0-10"
# Avoid errors if pngquant isn't available (perhaps on R-Forge?)
if (!nzchar(Sys.which("pngquant"))) pngquant <- NULL 

# Set chunk options
opts_chunk$set(
  tidy = FALSE,
  # Figure options
  fig.margin = TRUE, fig.width = 6, fig.height = 4,
  out.width = "100%", dpi = basedpi,
  cache = TRUE, pngquant = pngquant,
  # From "Tufte Handout" example dated 2016-12-27
  # Invalidate cache when the tufte version changes
  cache.extra = packageVersion('tufte')
)

## ----library_CHNOSZ---------------------------------------------------------------------------------------------------------------------------------------------------------------
library(CHNOSZ)

## ----info_CH4---------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Get database index for aqueous methane
info("CH4")

## ----info_methane-----------------------------------------------------------------------------------------------------------------------------------------------------------------
# Two ways to lookup methane gas
info("methane")
info("CH4", "gas")

## ----info_info_CH4----------------------------------------------------------------------------------------------------------------------------------------------------------------
# Get thermodynamic properties for aqueous methane
info(info("CH4"))

## ----info_fuzzy-------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Search for ribose-related entries
info("ribose+")

## ----subcrt_CH4-------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Properties of aqueous methane at default T and P
subcrt("CH4")

## ----subcrt_TP--------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Custom T, P grid for water in supercritical region
subcrt("H2O", T = c(400, 500), P = c(250, 300))

## ----E.units----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Change energy units from Joules to calories
E.units("cal")
subcrt("CH4", T = 25)
reset()  # Restore defaults

## ----subcrt_CO2_reaction----------------------------------------------------------------------------------------------------------------------------------------------------------
# CO2 dissolution reaction
subcrt(c("CO2", "CO2"), c("gas", "aq"), c(-1, 1), T = 25)

## ----basis_subcrt-----------------------------------------------------------------------------------------------------------------------------------------------------------------
basis(c("CO2", "H2O", "H+", "e-"))
subcrt(c("CH4", "acetate"), c("aq", "aq"), c(1, -1), T = 25)

## ----basis_species----------------------------------------------------------------------------------------------------------------------------------------------------------------
basis(c("CO2", "H2O", "H+", "e-"))
species(c("CH4", "acetate"))

## ----setup2, include = FALSE, cache = FALSE---------------------------------------------------------------------------------------------------------------------------------------
# Change the default options for the following chunks to hide results and messages
opts_chunk$set(
  results="hide", message=FALSE
)

## ----diagram, fig.cap = "Eh--pH (Pourbaix) diagram for S"-------------------------------------------------------------------------------------------------------------------------
# Setup the C-H-N-O-S basis system with electron
basis("CHNOSe")
# Define aqueous sulfur species
species(c("SO4-2", "HSO4-", "HS-", "H2S"))
# Calculate affinities on an Eh-pH grid
a <- affinity(pH = c(0, 14), Eh = c(-0.8, 1.2))
# Create an Eh-pH (Pourbaix) diagram
diagram(a, col = 4, col.names = 4, italic = TRUE)
# Add legend
TC <- convert(a$T, "C")
legend <- c(describe.property("T", TC), quote(italic(P) == "1 bar"))
legend("topright", legend = legend, bty = "n")

## ----mosaic, fig.cap = "Mosaic diagram showing effect of aqueous S speciation on the relative stabilities of Cu minerals and aqueous species"-------------------------------------
# Create a mosaic diagram for Cu-S-Cl-H2O system
basis(c("Cu", "H2S", "Cl-", "H2O", "H+", "e-"))
basis("H2S", -6)
basis("Cl-", -1)
species(c("CuCl", "CuCl2-", "CuCl3-2", "CuCl+", "CuCl2", "CuCl3-", "CuCl4-2"))
species(c("chalcocite", "tenorite", "cuprite", "copper"), add = TRUE)
bases <- c("H2S", "HS-", "HSO4-", "SO4-2")
m <- mosaic(bases, pH = c(0, 12), Eh = c(-1, 1), T = 200)
a <- m$A.species
diagram(a, lwd = 2, bold = species()$state == "cr")
diagram(m$A.bases, add = TRUE, col = 4, col.names = 4, italic = TRUE)
water.lines(m$A.species, lty = 3, col = 8)
# Add legend
legend <- lTP(convert(a$T, "C"), a$P)
legend("topright", legend = legend, bty = "n")

## ----equilibrate, fig.cap="Carbonate speciation as a function of pH and temperature"----------------------------------------------------------------------------------------------
# Carbonate speciation vs pH
basis(c("CO2", "H2O", "H+", "e-"))
species(c("CO2", "HCO3-", "CO3-2"))
# 25 degrees C
a <- affinity(pH = c(0, 14))
e <- equilibrate(a)
diagram(e, alpha = TRUE)
# 100 degrees C
a <- affinity(pH = c(4, 12), T = 100)
e <- equilibrate(a)
diagram(e, alpha = TRUE, add = TRUE, col = 2, names = NA)
# Add legend
legend <- as.expression(list(lT(25), lT(100)))
legend("left", legend = legend, lty = 1, col = c(1, 2))

## ----corundum_solubility, fig.cap="Solubility of corundum (green line) and equilibrium concentrations of aqueous species (black lines)"-------------------------------------------
# Corundum solubility vs pH
basis(c("Al+3", "H2O", "H+", "e-"))
species("corundum")
iaq <- info(c("Al(OH)4-", "Al(OH)3", "Al(OH)2+", "AlOH+2", "Al+3"))
s <- solubility(iaq, pH = c(2, 10), in.terms.of = "Al+3")
diagram(s, col = "green4", lwd = 2, ylim = c(-10, -2))
diagram(s, type = "loga.equil", add = TRUE)
legend("topright", c("25 °C", "1 bar"), bty = "n")

## ----corundum_solubility_IS, fig.cap="Solubility of corundum dependent on ionic strength"-----------------------------------------------------------------------------------------
# Corundum solubility again
basis(c("Al+3", "H2O", "H+", "e-"))
species("corundum")
iaq <- info(c("Al(OH)4-", "Al(OH)3", "Al(OH)2+", "AlOH+2", "Al+3"))
# Calculate with ionic strength of 0 and 1 molal
s0 <- solubility(iaq, pH = c(2, 10))
s1 <- solubility(iaq, pH = c(2, 10), IS = 1)
diagram(s1, col = 4, lwd = 3, ylim = c(-8, -2))
diagram(s0, col = "green4", lwd = 2, add = TRUE)
legend("topleft", legend = c(1, 0), lwd = c(3, 2), col = c(4, "green4"), title = "I (mol/kg)", bty = "n")
legend("topright", c("25 °C", "1 bar"), bty = "n")

## ----refs, eval = FALSE-----------------------------------------------------------------------------------------------------------------------------------------------------------
# # Return data frame with references for one or more species
# thermo.refs(info("CH4", c("aq", "gas")))
# # View all references in a browser
# thermo.refs()

## ----quasisolubility, echo=FALSE, fig.cap="The total solubility contour (green) at log *m*Fe = -6 lies inside the quasisolubility contour (boundary between tan and blue areas), showing that the latter underestimates solubility; the largest contributions to total solubility are from Fe^+2^ and FeCl^+^ at low pH and FeO~2~^-^ and HFeO~2~^-^ at high pH", fig.margin=FALSE, fig.width=12, fig.height=4.2, cache=TRUE----
par(mfrow = c(1, 2))

# System definition
metal <- "Fe"
# The concentration to be used for the quasisolubility contour
logm_metal <- -6
T <- 150
P <- "Psat"
Stot <- 1e-3
# Approximate chloride molality and ionic strength for mNaCl = 1.9
m_Clminus <- IS <- 1.5
# Plot variables
res <- 300
pH <- c(2, 12, res)
O2 <- c(-55, -40, res)

# Setup basis species
basis(c(metal, "H2S", "Cl-", "oxygen", "H2O", "H+"))
basis("H2S", log10(Stot))
basis("Cl-", log10(m_Clminus))
# Define basis species to swap through for mosaic calculation
bases <- c("H2S", "HS-", "HSO4-", "SO4-2")

# Get minerals and aqueous species
icr <- info(c("hematite", "magnetite", "pyrite", "pyrrhotite"))
iaq <- retrieve(metal, c("Cl", "S", "O", "H"), state = "aq")

## First plot: predominance diagram (quasisolubility) overlaid with total solubility contours
species(icr)
species(iaq, logm_metal, add = TRUE)
m <- mosaic(bases, pH = pH, O2 = O2, T = T, P = P, IS = IS)
diagram(m$A.species, bold = species()$state == "cr")

# Make total solubility contours
species(icr)
s <- solubility(iaq, bases = bases, pH = pH, O2 = O2, T = T, P = P, IS = IS, in.terms.of = metal)
levels <- seq(-9, -3, 3)
diagram(s, levels = levels, contour.method = "flattest", add = TRUE, col = "green4", lwd = 2)
# Add a legend
leg_list <- c(
  lTP(T, P),
  lNaCl(1.9),
  lS(Stot)
)
leg_expr <- lex(leg_list)
legend("bottomleft", legend = leg_expr, bg = "white")
title("Total solubility is higher than quasisolubility", font.main = 1)

## Second plot: activities of aqueous species in solubility calculation
basis("O2", -46)
s <- solubility(iaq, bases = bases, pH = pH, T = T, P = P, IS = IS, in.terms.of = metal)
diagram(s, col = "green4", lwd = 2, ylim = c(-10, -6))
# Locate the stability boundary between pyrite and magnetite along the the pH axis 
ipH <- which.min(abs(s$loga.equil[[3]] - s$loga.equil[[2]]))
# Show py-mag stability boundary
pH_py_mag <- s$vals$pH[ipH]
abline(v = pH_py_mag, col = 8)
text(pH_py_mag, -6, "pyrite", adj = c(1.2, 1.5), font = 2)
text(pH_py_mag, -6, "magnetite", adj = c(-0.1, 1.5), font = 2)

# Consider each mineral separately to get activities of aqueous speices
species("pyrite")
s_py <- solubility(iaq, bases = bases, pH = pH, T = T, P = P, IS = IS, in.terms.of = metal)
species("magnetite")
s_mag <- solubility(iaq, bases = bases, pH = pH, T = T, P = P, IS = IS, in.terms.of = metal)
# Merge activities for stable regions of each mineral
loga.equil_both <- lapply(1:length(s_py$loga.equil), function(ispecies) {
  c(s_py$loga.equil[[ispecies]][1:(ipH-1)], s_mag$loga.equil[[ispecies]][ipH:res])
})
# Put merged activities into new object for plotting
s_both <- s_py
s_both$loga.equil <- loga.equil_both
diagram(s_both, type = "loga.equil", add = TRUE)

leg_expr <- expr.species("O2", "gas", value = -46, log = TRUE)
legend("bottomleft", legend = leg_expr, bg = "white")
title("Species contributing to total solubility", font.main = 1)

## ----quasisolubility, eval=FALSE, cache=FALSE-------------------------------------------------------------------------------------------------------------------------------------
# par(mfrow = c(1, 2))
# 
# # System definition
# metal <- "Fe"
# # The concentration to be used for the quasisolubility contour
# logm_metal <- -6
# T <- 150
# P <- "Psat"
# Stot <- 1e-3
# # Approximate chloride molality and ionic strength for mNaCl = 1.9
# m_Clminus <- IS <- 1.5
# # Plot variables
# res <- 300
# pH <- c(2, 12, res)
# O2 <- c(-55, -40, res)
# 
# # Setup basis species
# basis(c(metal, "H2S", "Cl-", "oxygen", "H2O", "H+"))
# basis("H2S", log10(Stot))
# basis("Cl-", log10(m_Clminus))
# # Define basis species to swap through for mosaic calculation
# bases <- c("H2S", "HS-", "HSO4-", "SO4-2")
# 
# # Get minerals and aqueous species
# icr <- info(c("hematite", "magnetite", "pyrite", "pyrrhotite"))
# iaq <- retrieve(metal, c("Cl", "S", "O", "H"), state = "aq")
# 
# ## First plot: predominance diagram (quasisolubility) overlaid with total solubility contours
# species(icr)
# species(iaq, logm_metal, add = TRUE)
# m <- mosaic(bases, pH = pH, O2 = O2, T = T, P = P, IS = IS)
# diagram(m$A.species, bold = species()$state == "cr")
# 
# # Make total solubility contours
# species(icr)
# s <- solubility(iaq, bases = bases, pH = pH, O2 = O2, T = T, P = P, IS = IS, in.terms.of = metal)
# levels <- seq(-9, -3, 3)
# diagram(s, levels = levels, contour.method = "flattest", add = TRUE, col = "green4", lwd = 2)
# # Add a legend
# leg_list <- c(
#   lTP(T, P),
#   lNaCl(1.9),
#   lS(Stot)
# )
# leg_expr <- lex(leg_list)
# legend("bottomleft", legend = leg_expr, bg = "white")
# title("Total solubility is higher than quasisolubility", font.main = 1)
# 
# ## Second plot: activities of aqueous species in solubility calculation
# basis("O2", -46)
# s <- solubility(iaq, bases = bases, pH = pH, T = T, P = P, IS = IS, in.terms.of = metal)
# diagram(s, col = "green4", lwd = 2, ylim = c(-10, -6))
# # Locate the stability boundary between pyrite and magnetite along the the pH axis
# ipH <- which.min(abs(s$loga.equil[[3]] - s$loga.equil[[2]]))
# # Show py-mag stability boundary
# pH_py_mag <- s$vals$pH[ipH]
# abline(v = pH_py_mag, col = 8)
# text(pH_py_mag, -6, "pyrite", adj = c(1.2, 1.5), font = 2)
# text(pH_py_mag, -6, "magnetite", adj = c(-0.1, 1.5), font = 2)
# 
# # Consider each mineral separately to get activities of aqueous speices
# species("pyrite")
# s_py <- solubility(iaq, bases = bases, pH = pH, T = T, P = P, IS = IS, in.terms.of = metal)
# species("magnetite")
# s_mag <- solubility(iaq, bases = bases, pH = pH, T = T, P = P, IS = IS, in.terms.of = metal)
# # Merge activities for stable regions of each mineral
# loga.equil_both <- lapply(1:length(s_py$loga.equil), function(ispecies) {
#   c(s_py$loga.equil[[ispecies]][1:(ipH-1)], s_mag$loga.equil[[ispecies]][ipH:res])
# })
# # Put merged activities into new object for plotting
# s_both <- s_py
# s_both$loga.equil <- loga.equil_both
# diagram(s_both, type = "loga.equil", add = TRUE)
# 
# leg_expr <- expr.species("O2", "gas", value = -46, log = TRUE)
# legend("bottomleft", legend = leg_expr, bg = "white")
# title("Species contributing to total solubility", font.main = 1)

## ----dissolution_logK, fig.cap="Equilibrium constants of dissolution reactions", small_mar = TRUE---------------------------------------------------------------------------------
T <- seq(0, 350, 10)
CO2_logK <- subcrt(c("CO2", "CO2"), c("gas", "aq"), c(-1, 1), T = T)$out$logK
CO_logK <- subcrt(c("CO", "CO"), c("gas", "aq"), c(-1, 1), T = T)$out$logK
CH4_logK <- subcrt(c("CH4", "CH4"), c("gas", "aq"), c(-1, 1), T = T)$out$logK
logK <- data.frame(T, CO2_logK, CO_logK, CH4_logK)
matplot(logK[, 1], logK[, -1], type = "l", col = 1, lty = 1,
        xlab = axis.label("T"), ylab = axis.label("logK"))
text(80, -1.7, expr.species("CO2"))
text(240, -2.37, expr.species("CO"))
text(300, -2.57, expr.species("CH4"))
# Add legend
legend <- describe.property("P", "Psat")
legend("topright", legend = legend, bty = "n")

## ----retrieve_Al_default, results = "show"----------------------------------------------------------------------------------------------------------------------------------------
# List aqueous Al species in the default database
iaq <- retrieve("Al", ligands = c("F", "H", "O"), state = "aq")
# Print the first few rows and columns
info(iaq)[1:3, 1:5]
# Use the species index or name in subcrt()
subcrt(iaq[1:3], T = 25)
#subcrt(names(iaq)[1:3], T = 25)  # same as above

## ----Al_diagram-------------------------------------------------------------------------------------------------------------------------------------------------------------------
basis(c("Al+3", "H2O", "F-", "H+", "e-"))
species(iaq)
species(names(iaq))  # same as above

## ----add_OBIGT_SLOP98, results = "show"-------------------------------------------------------------------------------------------------------------------------------------------
add.OBIGT("SLOP98")
iaq_all <- retrieve("Al", ligands = c("F", "H", "O"), state = "aq")
# Use difference between sets to get the added species
info(setdiff(iaq_all, iaq))

## ----corundum_solubility_2, fig.cap="Corundum solubility with species from SLOP98"------------------------------------------------------------------------------------------------
# Add superseded species from SLOP98
add.OBIGT("SLOP98")
# List all aqueous Al species
iaq <- retrieve("Al", state = "aq")
# Keep only species from Shock et al. (1997)
iaq <- iaq[grepl("SSWS97", info(iaq)$ref1)]
# Plot corundum solubility vs pH
basis(c("Al+3", "H2O", "H+", "e-"))
species("corundum")
s <- solubility(iaq, pH = c(2, 10), in.terms.of = "Al+3")
## Alternatively, we could use the species names
#s <- solubility(names(iaq), pH = c(2, 10), in.terms.of = "Al+3")
diagram(s, col = "green4", lwd = 2, ylim = c(-10, -2))
diagram(s, type = "loga.equil", add = TRUE)
legend("topright", c("25 °C", "1 bar"), bty = "n")
# Reset the database for subsequent calculations
reset()

## ----basis_Al_F_OH_1, error=TRUE--------------------------------------------------------------------------------------------------------------------------------------------------
try({
basis(c("Al+3", "F-", "OH-"))
})

## ----basis_Al_F_OH_2, error=TRUE--------------------------------------------------------------------------------------------------------------------------------------------------
try({
basis(c("Al+3", "F-", "H+", "H2O"))
})

## ----basis_Al_F_OH_3--------------------------------------------------------------------------------------------------------------------------------------------------------------
# Use "oxygen" to get oxygen gas (for log fO2 diagrams)
basis(c("Al+3", "F-", "H+", "H2O", "oxygen"))
# Use "e-" to get aqueous electron (for Eh diagrams)
basis(c("Al+3", "F-", "H+", "H2O", "e-"))

## ----species_logact---------------------------------------------------------------------------------------------------------------------------------------------------------------
basis(c("Al+3", "F-", "H+", "H2O", "e-"))
iaq <- retrieve("Al", ligands = c("F", "H", "O"), state = "aq")
# Check that the data are from the same source
stopifnot(all(info(iaq)$ref1 == "TS01"))
species(iaq, -6)

## ----Pourbaix_Mn, fig.cap = "Pourbaix diagram for Mn with two quasisolubility contours"-------------------------------------------------------------------------------------------
basis(c("Mn+2", "H+", "H2O", "e-"))
icr <- retrieve("Mn", ligands = c("H", "O"), state = "cr")
iaq <- retrieve("Mn", ligands = c("H", "O"), state = "aq")
# First layer, logact(aq) = -3
species(icr)
species(iaq, add = TRUE)
a <- affinity(pH = c(4, 14), Eh = c(-1, 1), T = 100)
# Use names = NA to avoid plotting labels twice
diagram(a, lty = 2, names = NA)
# Second layer, logact(aq) = -4
species(icr)
species(iaq, -4, add = TRUE)
a <- affinity(pH = c(4, 14), Eh = c(-1, 1), T = 100)
d <- diagram(a, bold = species()$state == "cr", add = TRUE)
# Add water stability limits
water.lines(d, lty = 3, col = 8)
# Add legends
legend("topright", legend = c(lT(100), lP("Psat")), bty = "n")
title = as.expression(quote(log~italic(a)*"Mn(aq)"))
legend("bottomleft", legend = c(-3, -4), lty = c(2, 1), title = title, bty = "n")

## ----NaCl, results = "show"-------------------------------------------------------------------------------------------------------------------------------------------------------
T <- 300
P <- 1000
pH <- 5
m_NaCl = 0.8
NaCl <- NaCl(m_NaCl = m_NaCl, T = T, P = P, pH = pH)
print(paste("mol NaCl added to 1 kg H2O:", m_NaCl))
print(paste("Ionic strength (mol/kg):", NaCl$IS))
print(paste("Chloride concentration (mol/kg):", NaCl$m_Clminus))

## ----solubility, echo=FALSE, fig.cap="Stability diagram for minerals; predominance diagram for aqueous species; composite diagram with quasisolubility contour; diagram with total solubility contours in units of log *m*", fig.margin=FALSE, fig.fullwidth=TRUE, fig.width=16, fig.height=3, cache=TRUE, dpi=1.2*basedpi----
par(mfrow = c(1, 4))

# System definition
metal <- "Fe"
# The concentration to be used for the quasisolubility contour
logm_metal <- -6
T <- 150
P <- "Psat"
Stot <- 1e-3
wt_percent_NaCl <- 10
# Plot variables
res <- 300
pH <- c(0, 14, res)
O2 <- c(-55, -35, res)

# Work out NaCl molality from weight percent
# Convert to permil to get parts out of 1000 g (1 kg) of solution
wt_permil_NaCl <- wt_percent_NaCl * 10
# Divide by molecular weight to get moles of NaCl in 1000 g of solution
moles_NaCl <- wt_permil_NaCl / mass("NaCl")
# Subtract from 1000 g to get mass of H2O
grams_H2O <- 1000 - wt_permil_NaCl
# This gives the moles of NaCl added to 1 kg of H2O
m_NaCl <- 1000 * moles_NaCl / grams_H2O
# Now calculate ionic strength and molality of Cl-
NaCl_res <- NaCl(m_NaCl, T = T, P = P)
IS = NaCl_res$IS
m_Clminus = NaCl_res$m_Clminus

# Setup basis species
basis(c(metal, "H2S", "Cl-", "oxygen", "H2O", "H+"))
basis("H2S", log10(Stot))
basis("Cl-", log10(m_Clminus))
# Define basis species to swap through for mosaic calculation
bases <- c("H2S", "HS-", "HSO4-", "SO4-2")

# Get minerals and aqueous species
icr <- retrieve(metal, c("Cl", "S", "O", "H"), state = "cr")
iaq <- retrieve(metal, c("Cl", "S", "O", "H"), state = "aq")

# Make diagram for minerals
species(icr)
mcr <- mosaic(bases, pH = pH, O2 = O2, T = T, P = P, IS = IS)
diagram(mcr$A.species, bold = TRUE, fill = "#FAEBD788")
diagram(mcr$A.bases, lty = 2, col = 4, col.names = 4, italic = TRUE, add = TRUE)
title(paste(metal, "minerals"), font.main = 1)
# Add a legend
leg_list <- c(
  lTP(T, P),
  lNaCl(m_NaCl),
  lS(Stot)
)
leg_expr <- lex(leg_list)
legend("topright", legend = leg_expr, bty = "n")

# Make diagram for aqueous species
species(iaq)
maq <- mosaic(bases, pH = pH, O2 = O2, T = T, P = P, IS = IS)
diagram(maq$A.species, fill = "#F0F8FF88")
title(paste(metal, "aqueous species"), font.main = 1)

# Make diagram for minerals and aqueous species
species(icr)
species(iaq, logm_metal, add = TRUE)
m <- mosaic(bases, pH = pH, O2 = O2, T = T, P = P, IS = IS)
diagram(m$A.species, bold = species()$state == "cr")
diagram(m$A.bases, lty = 2, col = 4, col.names = 4, italic = TRUE, add = TRUE)
main = bquote("Quasisolubility contour for" ~ log ~ italic(m)*.(metal)*"(aq)" == .(logm_metal))
title(main, font.main = 1)

# Make total solubility contours
species(icr)
s <- solubility(iaq, bases = bases, pH = pH, O2 = O2, T = T, P = P, IS = IS, in.terms.of = metal)
levels <- seq(-12, 9, 3)
diagram(s, levels = levels, contour.method = "flattest", col = "green4")
title("Total solubility contours", font.main = 1)

## ----solubility, eval=FALSE, cache=FALSE------------------------------------------------------------------------------------------------------------------------------------------
# par(mfrow = c(1, 4))
# 
# # System definition
# metal <- "Fe"
# # The concentration to be used for the quasisolubility contour
# logm_metal <- -6
# T <- 150
# P <- "Psat"
# Stot <- 1e-3
# wt_percent_NaCl <- 10
# # Plot variables
# res <- 300
# pH <- c(0, 14, res)
# O2 <- c(-55, -35, res)
# 
# # Work out NaCl molality from weight percent
# # Convert to permil to get parts out of 1000 g (1 kg) of solution
# wt_permil_NaCl <- wt_percent_NaCl * 10
# # Divide by molecular weight to get moles of NaCl in 1000 g of solution
# moles_NaCl <- wt_permil_NaCl / mass("NaCl")
# # Subtract from 1000 g to get mass of H2O
# grams_H2O <- 1000 - wt_permil_NaCl
# # This gives the moles of NaCl added to 1 kg of H2O
# m_NaCl <- 1000 * moles_NaCl / grams_H2O
# # Now calculate ionic strength and molality of Cl-
# NaCl_res <- NaCl(m_NaCl, T = T, P = P)
# IS = NaCl_res$IS
# m_Clminus = NaCl_res$m_Clminus
# 
# # Setup basis species
# basis(c(metal, "H2S", "Cl-", "oxygen", "H2O", "H+"))
# basis("H2S", log10(Stot))
# basis("Cl-", log10(m_Clminus))
# # Define basis species to swap through for mosaic calculation
# bases <- c("H2S", "HS-", "HSO4-", "SO4-2")
# 
# # Get minerals and aqueous species
# icr <- retrieve(metal, c("Cl", "S", "O", "H"), state = "cr")
# iaq <- retrieve(metal, c("Cl", "S", "O", "H"), state = "aq")
# 
# # Make diagram for minerals
# species(icr)
# mcr <- mosaic(bases, pH = pH, O2 = O2, T = T, P = P, IS = IS)
# diagram(mcr$A.species, bold = TRUE, fill = "#FAEBD788")
# diagram(mcr$A.bases, lty = 2, col = 4, col.names = 4, italic = TRUE, add = TRUE)
# title(paste(metal, "minerals"), font.main = 1)
# # Add a legend
# leg_list <- c(
#   lTP(T, P),
#   lNaCl(m_NaCl),
#   lS(Stot)
# )
# leg_expr <- lex(leg_list)
# legend("topright", legend = leg_expr, bty = "n")
# 
# # Make diagram for aqueous species
# species(iaq)
# maq <- mosaic(bases, pH = pH, O2 = O2, T = T, P = P, IS = IS)
# diagram(maq$A.species, fill = "#F0F8FF88")
# title(paste(metal, "aqueous species"), font.main = 1)
# 
# # Make diagram for minerals and aqueous species
# species(icr)
# species(iaq, logm_metal, add = TRUE)
# m <- mosaic(bases, pH = pH, O2 = O2, T = T, P = P, IS = IS)
# diagram(m$A.species, bold = species()$state == "cr")
# diagram(m$A.bases, lty = 2, col = 4, col.names = 4, italic = TRUE, add = TRUE)
# main = bquote("Quasisolubility contour for" ~ log ~ italic(m)*.(metal)*"(aq)" == .(logm_metal))
# title(main, font.main = 1)
# 
# # Make total solubility contours
# species(icr)
# s <- solubility(iaq, bases = bases, pH = pH, O2 = O2, T = T, P = P, IS = IS, in.terms.of = metal)
# levels <- seq(-12, 9, 3)
# diagram(s, levels = levels, contour.method = "flattest", col = "green4")
# title("Total solubility contours", font.main = 1)

## ----convert----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Convert 100 degrees C to value in Kelvin
convert(100, "K")
# Convert 273.15 K to value in degrees C
convert(273.15, "C")
# Convert 1 bar to value in MPa
convert(1, "MPa")
# Convert 100 MPa to value in bar
convert(100, "bar")
# Convert 1 cal/mol to value in J/mol
convert(1, "J")
# Convert 10 J/mol to value in cal/mol
convert(10, "cal")

## ----convert_ppm, fig.cap = "Solubility in units of ppm"--------------------------------------------------------------------------------------------------------------------------
sppm <- convert(s, "ppm")
levels <- c(1e-6, 1e-3, 1e0, 1e3)
diagram(sppm, levels = levels, col = "green4")

## ----transect_setup---------------------------------------------------------------------------------------------------------------------------------------------------------------
basis("CHNOSe")
species(c("NO3-", "NO2-", "NH4+", "NH3"))
a <- affinity(pH = c(6, 8, 6, 8), Eh = c(0.5, 0.5, -0.5, -0.5))

## ----transect_results, results="show", cache=FALSE--------------------------------------------------------------------------------------------------------------------------------
# Print pH and Eh values used for calculation
a$vals
# Print affinity values calculated for each species
a$values

## ----rainbow, echo=FALSE, fig.cap="Affinities of organic synthesis reactions per mole of C, H<sub>2</sub>, or formed species", fig.margin=FALSE, fig.fullwidth=TRUE, fig.width=12, fig.height=3, cache=TRUE, dpi=1.2*basedpi----
basis(c("CO2", "H2", "NH4+", "H2O", "H2S", "H+"))
# Constant activity of CH4 is a simplification of the model
species("CH4", -3)
# Add other organic species with lower activity
species(c("propanoic acid", "adenine", "leucine", "deoxyribose"), -6, add = TRUE)
# Read file with variable conditions
file <- system.file("extdata/misc/SC10_Rainbow.csv", package = "CHNOSZ")
# Use `check.names = FALSE` to preserve the `NH4+` column name
df <- read.csv(file, check.names = FALSE)
# Reverse the rows to get increasing T
df <- df[rev(rownames(df)), ]
# Define six variables on a transect
a <- affinity(T = df$T, CO2 = df$CO2, H2 = df$H2,
              `NH4+` = df$`NH4+`, H2S = df$H2S, pH = df$pH)
# Get T in Kelvin to make unit conversions
T <- convert(a$vals[[1]], "K")
# Convert dimensionless affinity to Gibbs energy in units of J/mol
# nb. this is the same as converting logK to ΔG of reaction
a$values <- lapply(a$values, convert, "G", T)
# Convert J/mol to cal/mol
a$values <- lapply(a$values, convert, "cal")
# Convert cal/mol to kcal/mol
a$values <- lapply(a$values, `*`, -0.001)
# Setup figure
par(mfrow = c(1, 3))
par(cex = 0.9)
# First plot: balance on C
ylab = quote(italic(A)~"(kcal/mol) per C")
diagram(a, ylim = c(-20, 40), ylab  = ylab, col = 1:5)
abline(h = 0, lty = 2, lwd = 2)
# Second plot: balance on H2
ylab = quote(italic(A)~"(kcal/mol) per H2")
diagram(a, balance = "H2", ylim = c(-15, 10), ylab = ylab, col = 1:5)
abline(h = 0, lty = 2, lwd = 2)
# Third plot: no balancing
ylab <- quote(italic(A)~"(kcal/mol) per species")
diagram(a, balance = 1, ylim = c(-100, 100), ylab = ylab, col = 1:5)
abline(h = 0, lty = 2, lwd = 2)

## ----rainbow, eval=FALSE, cache=FALSE---------------------------------------------------------------------------------------------------------------------------------------------
# basis(c("CO2", "H2", "NH4+", "H2O", "H2S", "H+"))
# # Constant activity of CH4 is a simplification of the model
# species("CH4", -3)
# # Add other organic species with lower activity
# species(c("propanoic acid", "adenine", "leucine", "deoxyribose"), -6, add = TRUE)
# # Read file with variable conditions
# file <- system.file("extdata/misc/SC10_Rainbow.csv", package = "CHNOSZ")
# # Use `check.names = FALSE` to preserve the `NH4+` column name
# df <- read.csv(file, check.names = FALSE)
# # Reverse the rows to get increasing T
# df <- df[rev(rownames(df)), ]
# # Define six variables on a transect
# a <- affinity(T = df$T, CO2 = df$CO2, H2 = df$H2,
#               `NH4+` = df$`NH4+`, H2S = df$H2S, pH = df$pH)
# # Get T in Kelvin to make unit conversions
# T <- convert(a$vals[[1]], "K")
# # Convert dimensionless affinity to Gibbs energy in units of J/mol
# # nb. this is the same as converting logK to ΔG of reaction
# a$values <- lapply(a$values, convert, "G", T)
# # Convert J/mol to cal/mol
# a$values <- lapply(a$values, convert, "cal")
# # Convert cal/mol to kcal/mol
# a$values <- lapply(a$values, `*`, -0.001)
# # Setup figure
# par(mfrow = c(1, 3))
# par(cex = 0.9)
# # First plot: balance on C
# ylab = quote(italic(A)~"(kcal/mol) per C")
# diagram(a, ylim = c(-20, 40), ylab  = ylab, col = 1:5)
# abline(h = 0, lty = 2, lwd = 2)
# # Second plot: balance on H2
# ylab = quote(italic(A)~"(kcal/mol) per H2")
# diagram(a, balance = "H2", ylim = c(-15, 10), ylab = ylab, col = 1:5)
# abline(h = 0, lty = 2, lwd = 2)
# # Third plot: no balancing
# ylab <- quote(italic(A)~"(kcal/mol) per species")
# diagram(a, balance = 1, ylim = c(-100, 100), ylab = ylab, col = 1:5)
# abline(h = 0, lty = 2, lwd = 2)

## ----subcrt_nonideal, results = "show"--------------------------------------------------------------------------------------------------------------------------------------------
# Set the nonideality model
old_ni <- nonideal("Alberty")
# Calculate standard and adjusted Gibbs energy at 25 °C
species <- c("MgH2ATP", "MgHATP-", "MgATP-2")
subcrt(species, T = 25, IS = c(0, 0.25), property = "G")$out
# Restore the original nonideality model
nonideal(old_ni)

## ----subcrt_affinity, results = "show"--------------------------------------------------------------------------------------------------------------------------------------------
species <- c("Mg+2", "ATP-4", "MgATP-2")
coeffs <- c(-1, -1, 1)
T <- c(25, 50, 100)
# Drop some columns for more compact output
idrop <- c(2:4, 6:9)
# Unit activity: affinity is opposite of standard Gibbs energy
subcrt(species, coeffs, T = T, logact = c(0, 0, 0))$out[, -idrop]
# Non-unit activity: affinity is opposite of non-standard Gibbs energy
subcrt(species, coeffs, T = T, logact = c(-5, -4, -3))$out[, -idrop]

## ----organic, echo=FALSE, fig.cap="Non-standard Gibbs energies of organic reactions as a function of CO<sub>2</sub> fugacity", cache=TRUE-----------------------------------------
# Format reactions with describe.reaction()
basis(c("CO2", "H2", "H2O", "H+"))
reactions <- c(
  # hydrogenotrophic methanogenesis
  describe.reaction(subcrt("CH4", 1)$reaction),
  # acetate oxidation
  describe.reaction(subcrt("acetate", -1)$reaction),
  # acetoclastic methanogenesis
  describe.reaction(subcrt(c("acetate", "CH4"), c(-1, 1))$reaction)
)
# Add space before "CO2" to avoid getting cut off to "O2" by the plot limits
reactions[[1]][2][[1]][[2]][[2]][[2]][[3]] <- "       CO"

# System 1: one C species in basis; two C species in formation reactions
basis(c("CO2", "H2", "H2O", "H+"))
# Change CO2 and H2 to gas
basis(c("CO2", "H2"), "gas")
# Set H2 (log fugacity) and pH
basis(c("H2", "pH"), c(-3.92, 7.3))
# Write reactions to form methane (gas) and acetate (aq) at given activity and fugacity
species(c("methane", "acetate"), c(-0.18, -3.4))
# Calculate affinities of formation reactions with varying log fugacity of CO2
a1 <- affinity(CO2 = c(-3, 3), T = 55, P = 50)

# System 2: two C species in basis; three C species in formation reactions
# Swap H2 for CH4 and set fugacity
swap.basis("H2", "CH4")
basis("CH4", "gas")
basis("CH4", -0.18)
# Remove methane as formed species - we only want acetate now
species(1, delete = TRUE)
a2 <- affinity(CO2 = c(-3, 3), T = 55, P = 50)

# Convert from dimensionless affinity to Gibbs energy in kJ/mol
TK <- convert(55, "K")
a1$values[[1]] <- convert(a1$values[[1]], "G", T = TK) / 1000
# Negate to reverse reaction (acetate as a reactant)
a1$values[[2]] <- -convert(a1$values[[2]], "G", T = TK) / 1000
a2$values[[1]] <- -convert(a2$values[[1]], "G", T = TK) / 1000

# Make diagram without balancing constraint (no normalization by C)
diagram(a1, balance = 1, lty = 1, lwd = 2, col = c(4, 2), names = reactions[1:2], ylab = axis.label("DG", prefix = "k"), srt = c(-13.5, 26))
diagram(a2, balance = 1, lwd = 2, col = 3, names = reactions[3], srt = 14, add = TRUE)
abline(h = 0, lty = 2, col = 8)

## ----organic, eval=FALSE, cache=FALSE---------------------------------------------------------------------------------------------------------------------------------------------
# # Format reactions with describe.reaction()
# basis(c("CO2", "H2", "H2O", "H+"))
# reactions <- c(
#   # hydrogenotrophic methanogenesis
#   describe.reaction(subcrt("CH4", 1)$reaction),
#   # acetate oxidation
#   describe.reaction(subcrt("acetate", -1)$reaction),
#   # acetoclastic methanogenesis
#   describe.reaction(subcrt(c("acetate", "CH4"), c(-1, 1))$reaction)
# )
# # Add space before "CO2" to avoid getting cut off to "O2" by the plot limits
# reactions[[1]][2][[1]][[2]][[2]][[2]][[3]] <- "       CO"
# 
# # System 1: one C species in basis; two C species in formation reactions
# basis(c("CO2", "H2", "H2O", "H+"))
# # Change CO2 and H2 to gas
# basis(c("CO2", "H2"), "gas")
# # Set H2 (log fugacity) and pH
# basis(c("H2", "pH"), c(-3.92, 7.3))
# # Write reactions to form methane (gas) and acetate (aq) at given activity and fugacity
# species(c("methane", "acetate"), c(-0.18, -3.4))
# # Calculate affinities of formation reactions with varying log fugacity of CO2
# a1 <- affinity(CO2 = c(-3, 3), T = 55, P = 50)
# 
# # System 2: two C species in basis; three C species in formation reactions
# # Swap H2 for CH4 and set fugacity
# swap.basis("H2", "CH4")
# basis("CH4", "gas")
# basis("CH4", -0.18)
# # Remove methane as formed species - we only want acetate now
# species(1, delete = TRUE)
# a2 <- affinity(CO2 = c(-3, 3), T = 55, P = 50)
# 
# # Convert from dimensionless affinity to Gibbs energy in kJ/mol
# TK <- convert(55, "K")
# a1$values[[1]] <- convert(a1$values[[1]], "G", T = TK) / 1000
# # Negate to reverse reaction (acetate as a reactant)
# a1$values[[2]] <- -convert(a1$values[[2]], "G", T = TK) / 1000
# a2$values[[1]] <- -convert(a2$values[[1]], "G", T = TK) / 1000
# 
# # Make diagram without balancing constraint (no normalization by C)
# diagram(a1, balance = 1, lty = 1, lwd = 2, col = c(4, 2), names = reactions[1:2], ylab = axis.label("DG", prefix = "k"), srt = c(-13.5, 26))
# diagram(a2, balance = 1, lwd = 2, col = 3, names = reactions[3], srt = 14, add = TRUE)
# abline(h = 0, lty = 2, col = 8)

## ----organic_reactions, results="show"--------------------------------------------------------------------------------------------------------------------------------------------
a1$species
a2$species

## ----metastable_sulfur, fig.cap = "Eh--pH diagram for metastable S species"-------------------------------------------------------------------------------------------------------
basis("CHNOSe")
iaq <- retrieve("S", c("H", "O"), state = "aq")
species(iaq)
a <- affinity(pH = c(0, 14), Eh = c(-0.8, 1.2))
# Save results but don't plot the first diagram
d <- diagram(a, plot.it = FALSE)
# Remove the stable species
istable <- unique(as.numeric(d$predominant))
species(istable, delete = TRUE)
# Make a diagram for the metastable species
a <- affinity(pH = c(0, 14), Eh = c(-0.8, 1.2))
d <- diagram(a, col = 4, col.names = 4, italic = TRUE)

## ----makeup, results="show"-------------------------------------------------------------------------------------------------------------------------------------------------------
# Element count of a species in the database
iCli <- info("carrollite")
makeup(iCli)
# Sum the elements of formulas supplied as character strings
summed_elements <- makeup(c("CO2", "CH4"), sum = TRUE)
# Use the result to write a new chemical formula
as.chemical.formula(summed_elements)

## ----thermo, results="show"-------------------------------------------------------------------------------------------------------------------------------------------------------
str(thermo(), max.level = 1)
str(thermo()$OBIGT)

## ----thermo_T.units, results="show"-----------------------------------------------------------------------------------------------------------------------------------------------
# This has the same effect as T.units("K")
thermo("opt$T.units" = "K")
# Return to default
thermo("opt$T.units" = "C")

## ----buffer_1---------------------------------------------------------------------------------------------------------------------------------------------------------------------
# View available buffers
thermo()$buffer
# Define a new buffer or modify an existing one
mod.buffer("PPM", c("pyrite", "pyrrhotite", "magnetite"), "cr", 0)
# Buffer made of one species (acetic acid with log activity = -10)
mod.buffer("AC", "acetic acid", "aq", -10)

## ----buffer_2---------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Specify basis species
basis(c("FeCHNOS"))
# Associate O2 with the PPM buffer
basis("O2", "PPM")
# Calculate and retrieve the buffered fugacity of O2
a <- affinity(T = 200, P = 2000, return.buffer = TRUE)
# Access the buffered O2 fugacity (log fO2)
log_fO2 <- a$O2  # -44.28
# Calculate buffered fugacities across temperature range
a <- affinity(T = c(200, 400, 11), P = 2000, return.buffer = TRUE)

## ----buffer_3---------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Setup basis species
basis(c("FeS2", "H2S", "O2", "H2O"), c(0, -10, -50, 0))
# Associate both H2S and O2 with the PPM buffer
basis(c("H2S", "O2"), c("PPM", "PPM"))
# Retrieve values for both buffered species
buffer_values <- affinity(T = 300, P = 2000, return.buffer = TRUE)
log_aH2S <- buffer_values$H2S  #  -2.57
log_fO2 <- buffer_values$O2    # -37.16

## ----buffer_4, fig.cap = "Equilibrium log&nbsp;H<sub>2</sub> fugacity for 10^-6^ activity of HCN or formaldehyde with water, 1 bar of N<sub>2</sub> and 10 bar of CO<sub>2</sub>"----
# Setup a system in terms of gases and liquid water
basis(c("hydrogen", "carbon dioxide", "nitrogen", "water"))
# Use 10 bar of CO2 and 1 bar of other gases (default)
basis("CO2", 1)
# Load aqueous species with given log activity
species(c("HCN", "formaldehyde"), c(-6, -6))
# Calculate affinities to form aqueous species from basis species
a <- affinity(T = c(0, 350), P = 300)
# Create diagram showing H2 activity where affinity = 0
d <- diagram(a, type = "H2", lty = c(2, 3))
legend("bottomright", c("HCN", "formaldehyde"), lty = c(2, 3))

## ----buffer_5, fig.cap = "Gold solubility at 300 °C with PPM buffer for *f*O<sub>2</sub> and *a*H<sub>2</sub>S"-------------------------------------------------------------------
# Setup system for gold solubility calculation
basis(c("Au", "Fe", "H2S", "H2O", "oxygen", "H+"))
# Apply PPM buffer for fO2 and H2S
basis("O2", "PPM")
basis("H2S", "PPM")
# Calculate gold solubility using the buffered values
species("Au")
iaq <- info(c("Au(HS)2-", "AuHS", "AuOH"))
s <- solubility(iaq, pH = c(2, 8), T = 300, P = 1000)
# Create solubility diagram
diagram(s, ylim = c(-10, -5), col = "green4", lwd = 2)
col <- c("#ED4037", "#F58645", "#0F9DE2")  # Au(HS)2-, AuHS, AuOH
diagram(s, type = "loga.equil", add = TRUE, col = col)

## ----buffer_6---------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Set basis species
basis(c("Fe", "SiO2", "CO3-2", "H2O", "oxygen", "H+"))
# Calculate log fO2 in QFM buffer across temperature range
basis("O2", "QFM")
T <- seq(600, 1000, 100)
buf <- affinity(T = T, P = 5000, return.buffer = TRUE)
# Use buffered fO2 values in downstream calculations
species(c("CH4", "CO2", "HCO3-", "CO3-2"))
# Set values of pH for transect
pH <- seq(3.8, 4.3, length.out = length(T))
# Adding 2 log units below QFM buffer
a <- affinity(T = T, O2 = buf$O2 - 2, pH = pH, P = 5000)

## ----buffer_7---------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Calculate neutral pH at 300°C and 1000 bar
T <- 300
P <- 1000
neutral_pH <- -subcrt(c("H2O", "H+", "OH-"), c(-1, 1, 1), T = T, P = P)$out$logK/2

## ----buffer_8---------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Define the KMQ buffer
mod.buffer("KMQ", c("quartz", "muscovite", "K-feldspar"), "cr", 0)
# Setup the system
basis(c("Al2O3", "quartz", "K+", "H2O", "oxygen", "H+"))
# Set K+ molality for KCl solution
basis("K+", log10(0.5))
# Associate H+ with the KMQ buffer
basis("H+", "KMQ")
# Calculate buffered pH
pH_KMQ <- -affinity(T = 300, P = 1000, return.buffer = TRUE)$`H+`

## ----buffer_gold, echo=FALSE, fig.cap="Effects of different buffers on gold solubility", fig.margin=FALSE, fig.fullwidth=TRUE, fig.width=12, fig.height=3, cache=TRUE, dpi=1.2*basedpi----
# Setup the system
basis(c("Au", "Al2O3", "quartz", "Fe", "K+", "Cl-", "H2S", "H2O", "oxygen", "H+"))
# Set log activities for 0.5 m NaCl + 0.1 m KCl solution (assuming ideality for now)
basis("Cl-", log10(0.6))
basis("K+", log10(0.1))
# Set 0.01 m H2S
basis("H2S", -2)
# Setup the KMQ pH buffer
mod.buffer("KMQ", c("quartz", "muscovite", "K-feldspar"), "cr", 0)
basis("H+", "KMQ")
# Define colors for gold species
col <- c("#ED4037", "#F58645", "#0F9DE2", "#22CC88")  # Au(HS)2-, AuHS, AuOH, AuCl2-
# Create temperature transect with HM and PPM redox buffers
temps <- seq(200, 500, 10)
P <- 1000
# Initialize plot layout
par(mfrow = c(1, 4), mar = c(4, 4, 3, 1), font.main = 1, cex.main = 1.1)

# 1. Compare gold species distribution under HM buffer
basis("O2", "HM")
species("Au")
iaq <- info(c("Au(HS)2-", "AuHS", "AuOH", "AuCl2-"))
s_HM <- solubility(iaq, T = temps, P = P)
diagram(s_HM, xlab = "Temperature (°C)", ylab = "log activity", ylim = c(-10, -4),
        col = "green4", lwd = 2, main = "Solubility vs T (KMQ: pH, HM: fO2)")
diagram(s_HM, type = "loga.equil", col = col, add = TRUE)

# 2. Compare gold species distribution under PPM buffer
basis("O2", "PPM")
basis("H2S", "PPM")
s_PPM <- solubility(iaq, T = temps, P = P)
diagram(s_PPM, xlab = "Temperature (°C)", ylab = "log activity", ylim = c(-10, -4),
        col = "green4", lwd = 2, main = "Solubility vs T (KMQ: pH, PPM: fO2 and aH2S)")
diagram(s_PPM, type = "loga.equil", col = col, add = TRUE)

# 3. Gold solubility as function of pH at 350°C under PPM buffer
# Remove KMQ buffer
basis("H+", 0)
T_fixed <- 350
# Calculate solubility across pH range
s_pH <- solubility(iaq, pH = c(2, 8), T = T_fixed, P = P)
diagram(s_pH, xlab = "pH", ylab = "log activity", ylim = c(-10, -4),
        col = "green4", lwd = 2, main = "Solubility vs pH (350 °C, PPM: fO2 and aH2S)")
diagram(s_pH, type = "loga.equil", col = col, add = TRUE)
# Get neutral pH at this temperature
neutral_pH <- -subcrt(c("H2O", "H+", "OH-"), c(-1, 1, 1), T = T_fixed, P = P)$out$logK/2
# Add neutral pH line
abline(v = neutral_pH, lty = 2)
text(neutral_pH + 0.2, -5.5, "neutral pH", srt = 90, cex = 0.8)
# Add KMQ pH buffer line
basis("H+", "KMQ")
pH_KMQ <- -affinity(T = T_fixed, P = P, return.buffer = TRUE)$`H+`
abline(v = pH_KMQ, lty = 3, col = 4)
text(pH_KMQ - 0.2, -5.5, "KMQ", srt = 90, col = 4, cex = 0.8)

# 4. log fO2-pH diagram with dominant gold species at 350°C
# Remove KMQ pH buffer
basis("H+", 0)
# Remove PPM fO2 buffer
basis("O2", 0)
# Use constant aH2S from PPM buffer
a <- affinity(pH = c(2, 8), T = T, P = P, return.buffer = TRUE)
log_aH2S <- unique(a$H2S)
basis("H2S", log_aH2S)
# Load aqueous species for relative stability calculation
species(iaq)
# Define basis species to swap through for mosaic calculation
bases <- c("H2S", "HS-", "SO4-2", "HSO4-")
m <- mosaic(bases = bases, pH = c(2, 8), O2 = c(-32, -24), T = T_fixed, P = P)
a <- m$A.species
fill <- adjustcolor(col, alpha.f = 0.5)
d <- diagram(a, fill = fill, main = "Predominance fO2 vs pH (350°C, PPM: aH2S)")
water.lines(d, lty = 3)
# Add mineral buffer lines
# PPM buffer line
basis("O2", "PPM")
O2_PPM <- affinity(T = T_fixed, P = P, return.buffer = TRUE)$O2
abline(h = O2_PPM, lty = 2, col = "darkred")
text(2.5, O2_PPM, "PPM", adj = c(0, -0.5), col = "darkred", cex = 0.8)
# HM buffer line
basis("O2", "HM")
O2_HM <- affinity(T = T_fixed, P = P, return.buffer = TRUE)$O2
abline(h = O2_HM, lty = 2, col = "darkred")
text(2.5, O2_HM, "HM", adj = c(0, -0.5), col = "darkred", cex = 0.8)
# Add KMQ pH buffer line and neutral pH
abline(v = neutral_pH, lty = 2)
text(neutral_pH + 0.2, -28, "neutral pH", srt = 90, cex = 0.8)
abline(v = pH_KMQ, lty = 2, col = 4)
text(pH_KMQ - 0.2, -28, "KMQ", srt = 90, col = 4, cex = 0.8)

## ----buffer_gold, eval=FALSE, cache=FALSE-----------------------------------------------------------------------------------------------------------------------------------------
# # Setup the system
# basis(c("Au", "Al2O3", "quartz", "Fe", "K+", "Cl-", "H2S", "H2O", "oxygen", "H+"))
# # Set log activities for 0.5 m NaCl + 0.1 m KCl solution (assuming ideality for now)
# basis("Cl-", log10(0.6))
# basis("K+", log10(0.1))
# # Set 0.01 m H2S
# basis("H2S", -2)
# # Setup the KMQ pH buffer
# mod.buffer("KMQ", c("quartz", "muscovite", "K-feldspar"), "cr", 0)
# basis("H+", "KMQ")
# # Define colors for gold species
# col <- c("#ED4037", "#F58645", "#0F9DE2", "#22CC88")  # Au(HS)2-, AuHS, AuOH, AuCl2-
# # Create temperature transect with HM and PPM redox buffers
# temps <- seq(200, 500, 10)
# P <- 1000
# # Initialize plot layout
# par(mfrow = c(1, 4), mar = c(4, 4, 3, 1), font.main = 1, cex.main = 1.1)
# 
# # 1. Compare gold species distribution under HM buffer
# basis("O2", "HM")
# species("Au")
# iaq <- info(c("Au(HS)2-", "AuHS", "AuOH", "AuCl2-"))
# s_HM <- solubility(iaq, T = temps, P = P)
# diagram(s_HM, xlab = "Temperature (°C)", ylab = "log activity", ylim = c(-10, -4),
#         col = "green4", lwd = 2, main = "Solubility vs T (KMQ: pH, HM: fO2)")
# diagram(s_HM, type = "loga.equil", col = col, add = TRUE)
# 
# # 2. Compare gold species distribution under PPM buffer
# basis("O2", "PPM")
# basis("H2S", "PPM")
# s_PPM <- solubility(iaq, T = temps, P = P)
# diagram(s_PPM, xlab = "Temperature (°C)", ylab = "log activity", ylim = c(-10, -4),
#         col = "green4", lwd = 2, main = "Solubility vs T (KMQ: pH, PPM: fO2 and aH2S)")
# diagram(s_PPM, type = "loga.equil", col = col, add = TRUE)
# 
# # 3. Gold solubility as function of pH at 350°C under PPM buffer
# # Remove KMQ buffer
# basis("H+", 0)
# T_fixed <- 350
# # Calculate solubility across pH range
# s_pH <- solubility(iaq, pH = c(2, 8), T = T_fixed, P = P)
# diagram(s_pH, xlab = "pH", ylab = "log activity", ylim = c(-10, -4),
#         col = "green4", lwd = 2, main = "Solubility vs pH (350 °C, PPM: fO2 and aH2S)")
# diagram(s_pH, type = "loga.equil", col = col, add = TRUE)
# # Get neutral pH at this temperature
# neutral_pH <- -subcrt(c("H2O", "H+", "OH-"), c(-1, 1, 1), T = T_fixed, P = P)$out$logK/2
# # Add neutral pH line
# abline(v = neutral_pH, lty = 2)
# text(neutral_pH + 0.2, -5.5, "neutral pH", srt = 90, cex = 0.8)
# # Add KMQ pH buffer line
# basis("H+", "KMQ")
# pH_KMQ <- -affinity(T = T_fixed, P = P, return.buffer = TRUE)$`H+`
# abline(v = pH_KMQ, lty = 3, col = 4)
# text(pH_KMQ - 0.2, -5.5, "KMQ", srt = 90, col = 4, cex = 0.8)
# 
# # 4. log fO2-pH diagram with dominant gold species at 350°C
# # Remove KMQ pH buffer
# basis("H+", 0)
# # Remove PPM fO2 buffer
# basis("O2", 0)
# # Use constant aH2S from PPM buffer
# a <- affinity(pH = c(2, 8), T = T, P = P, return.buffer = TRUE)
# log_aH2S <- unique(a$H2S)
# basis("H2S", log_aH2S)
# # Load aqueous species for relative stability calculation
# species(iaq)
# # Define basis species to swap through for mosaic calculation
# bases <- c("H2S", "HS-", "SO4-2", "HSO4-")
# m <- mosaic(bases = bases, pH = c(2, 8), O2 = c(-32, -24), T = T_fixed, P = P)
# a <- m$A.species
# fill <- adjustcolor(col, alpha.f = 0.5)
# d <- diagram(a, fill = fill, main = "Predominance fO2 vs pH (350°C, PPM: aH2S)")
# water.lines(d, lty = 3)
# # Add mineral buffer lines
# # PPM buffer line
# basis("O2", "PPM")
# O2_PPM <- affinity(T = T_fixed, P = P, return.buffer = TRUE)$O2
# abline(h = O2_PPM, lty = 2, col = "darkred")
# text(2.5, O2_PPM, "PPM", adj = c(0, -0.5), col = "darkred", cex = 0.8)
# # HM buffer line
# basis("O2", "HM")
# O2_HM <- affinity(T = T_fixed, P = P, return.buffer = TRUE)$O2
# abline(h = O2_HM, lty = 2, col = "darkred")
# text(2.5, O2_HM, "HM", adj = c(0, -0.5), col = "darkred", cex = 0.8)
# # Add KMQ pH buffer line and neutral pH
# abline(v = neutral_pH, lty = 2)
# text(neutral_pH + 0.2, -28, "neutral pH", srt = 90, cex = 0.8)
# abline(v = pH_KMQ, lty = 2, col = 4)
# text(pH_KMQ - 0.2, -28, "KMQ", srt = 90, col = 4, cex = 0.8)

## ----protein_1--------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Search by name in thermo()$protein
ip1 <- pinfo("LYSC_CHICK")       # Using protein_organism format
ip2 <- pinfo("LYSC", "CHICK")    # Using separate arguments
# Search for the same protein in different organisms
ip3 <- pinfo("MYG", c("HORSE", "PHYCA"))

## ----protein_2--------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Reading amino acid compositions from a CSV file
file <- system.file("extdata/protein/POLG.csv", package = "CHNOSZ")
aa <- read.csv(file)
# Add the proteins to CHNOSZ and get their indices
iprotein <- add.protein(aa)

# For FASTA files, use the canprot package
fasta_file <- system.file("extdata/protein/rubisco.fasta", package = "CHNOSZ")
aa <- canprot::read_fasta(fasta_file)
iprotein <- add.protein(aa)

## ----protein_3, results = "show"--------------------------------------------------------------------------------------------------------------------------------------------------
# Get protein length (number of amino acids)
pl <- protein.length("LYSC_CHICK")
# Get chemical formula
pf <- protein.formula("LYSC_CHICK")
# Display results
list(length = pl, formula = pf)

## ----protein_4--------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Calculate ZC for a protein
ZC_value <- ZC(protein.formula("LYSC_CHICK"))
# For multiple proteins
proteins <- c("LYSC_CHICK", "MYG_HORSE", "RNAS1_BOVIN")
ZC_values <- ZC(protein.formula(proteins))

## ----protein_5, results = "show"--------------------------------------------------------------------------------------------------------------------------------------------------
# Properties of non-ionized protein
subcrt("LYSC_CHICK")$out[[1]][1:6, ]

## ----protein_6--------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Calculate ionization properties
aa <- pinfo(pinfo("LYSC_CHICK"))
charge <- ionize.aa(aa, pH = c(4, 7, 9))
# Calculate heat capacity of ionization
Cp_ionization <- ionize.aa(aa, property = "Cp", pH = 7, T = c(25, 100))

## ----protein_7--------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Define the basis species with H+
basis("CHNOS+")
# Add proteins to the system
species(c("LYSC_CHICK", "MYG_HORSE", "RNAS1_BOVIN"))

## ----protein_8--------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Calculate affinity as a function of pH
basis("CHNOS+")
species(c("LYSC_CHICK", "MYG_HORSE", "RNAS1_BOVIN"))
a1 <- affinity(pH = c(0, 14))

## ----protein_9--------------------------------------------------------------------------------------------------------------------------------------------------------------------
species(delete = TRUE)
ip <- pinfo(c("LYSC_CHICK", "MYG_HORSE", "RNAS1_BOVIN"))
# Set logarithm of activity with loga.protein
a2 <- affinity(pH = c(0, 14), iprotein = ip, loga.protein = -3)

# Check that both methods produce equivalent results
for(i in 1:3) stopifnot(all.equal(a1$values[[i]], a2$values[[i]]))

## ----protein_10, fig.cap = "Metastable equilibrium of proteins from hot-spring metagenomes"---------------------------------------------------------------------------------------
# Calculate equilibrium distribution as a function of Eh
basis("CHNOSe")
proteins <- paste("overall", paste0("bison", c("N", "S", "R", "Q", "P")), sep = "_")
ip <- pinfo(proteins)
a <- affinity(Eh = c(-0.35, -0.15), iprotein = ip)
# Normalize by protein length to get residue-equivalent distribution
# Set loga.balance to distribute proteins with total activity of residues equal to 1
e <- equilibrate(a, normalize = TRUE, loga.balance = 0)
col <- c("darkred", "red", "darkgray", "blue", "darkblue")
diagram(e, ylim = c(-4, -2.5), col = col, lwd = 2)
legend("bottomleft", c("High-T reducing", NA, NA, NA, "Low-T oxidizing"),
  lty = 1:5, col = col, title = "Environmental Conditions", inset = c(0.1, 0))

## ----protein_11-------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Sample protein data from YeastGFP study
protein <- c("YDL195W", "YHR098C", "YLR208W", "YNL049C", "YPL085W")
abundance <- c(1840, 12200, 21400, 1720, 358)
# Find protein indices in CHNOSZ database
ip <- match(protein, thermo()$protein$protein)
# Get protein lengths
pl <- protein.length(ip)
# Scale protein abundance so total abundance of residues is unity
scaled_abundance <- 10^unitize(log10(abundance), pl)
# Check that sum for residues is unity
stopifnot(all.equal(sum(scaled_abundance * pl), 1))

## ----protein_12, fig.cap = "Unoptimized predicted abundances of proteins compared to experimental abundances, both scaled to unit total activity of residues; dashed line is 1:1 line", small_mar = TRUE----
basis("CHNOSe")
# Make a guess for Eh
basis("Eh", -0.5)
a <- affinity(iprotein = ip)
e <- equilibrate(a, normalize = TRUE, loga.balance = 0)
# Check for unit total abundance of residues
predicted_abundance <- 10^unlist(e$loga.equil)
stopifnot(all.equal(sum(predicted_abundance * pl), 1))
plot(log10(scaled_abundance), log10(predicted_abundance), pch = 19)
# Show 1:1 line
lims <- range(par("usr"))
lines(c(lims[1], lims[2]), c(lims[1], lims[2]), lty = 2)
MAE <- mean(abs(log10(scaled_abundance) - log10(predicted_abundance)))
legend("topleft", paste("MAE =", round(MAE, 2)), bty = "n")

## ----protein_13-------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Load canprot package
library(canprot)
# Get amino acid compositions
ip <- pinfo(c("LYSC_CHICK", "MYG_HORSE", "RNAS1_BOVIN"))
aa <- pinfo(ip)
# canprot has Zc(); CHNOSZ has ZC()
Zc(aa)
# Stoichiometric oxygen and water content
nO2(aa)
nH2O(aa)
# Isoelectric point and GRAVY
pI(aa)
GRAVY(aa)

## ----protein_optimization, echo=FALSE, fig.cap = "Optimizing redox potential to fit experimental protein abundances", fig.margin=FALSE, fig.fullwidth=TRUE, fig.width=12, fig.height=3, cache=TRUE, dpi=1.2*basedpi----
# Protein data from YeastGFP study
protein <- c("YDL195W", "YHR098C", "YLR208W", "YNL049C", "YPL085W")
abundance <- c(1840, 12200, 21400, 1720, 358)
# Find protein indices in CHNOSZ database
ip <- match(protein, thermo()$protein$protein)
# Calculate protein lengths
pl <- protein.length(ip)
# Scale experimental abundances to unit total abundance of residues
scaled_logabundance <- unitize(log10(abundance), pl)

# Assign colors from low to high oxidation state (red to blue)
col_in <- c("darkred", "red", "darkgray", "blue", "darkblue")
aa <- pinfo(ip)
Zc_values <- canprot::Zc(aa)
col <- col_in[rank(Zc_values)]

# Setup chemical system with H+ and e-
basis("CHNOSe")
# Calculate affinities as a function of redox potential
a <- affinity(Eh = c(-0.6, 0), iprotein = ip)

# Calculate equilibrium distributions without normalization
e <- equilibrate(a, loga.balance = 0)
# Plot results
par(mfrow = c(1, 4))
diagram(e, ylim = c(-5, -2), col = col, lwd = 2, main = "Without normalization")

# With normalization
e_norm <- equilibrate(a, normalize = TRUE, loga.balance = 0)
diagram(e_norm, ylim = c(-5, -2.5), col = col, lwd = 2, main = "With normalization")

# Calculate and plot mean absolute error as a function of Eh
predicted_logabundances <- sapply(e_norm$loga.equil, c)
MAE <- apply(abs(t(predicted_logabundances) - scaled_logabundance), 2, mean)
par(xaxs = "r", yaxs = "r")
plot(a$vals$Eh, MAE)
imin <- which.min(MAE)
abline(v = a$vals$Eh[imin], lty = 2)
title("Optimized model Eh")

# Show calculated and experimental abundances
optimized_logabundance <- predicted_logabundances[imin, ]
plot(scaled_logabundance, optimized_logabundance, pch = 19, col = col)
# Show 1:1 line
lims <- range(par("usr"))
lines(c(lims[1], lims[2]), c(lims[1], lims[2]), lty = 2)
title("Predicted vs experimental abundance")
# Add legends
legend("topleft", paste("MAE =", round(MAE[imin], 2)), bty = "n")
legend("bottomright", c("High Zc", NA, NA, NA, "Low Zc"), pch = 19, col = rev(col_in))

## ----protein_optimization, eval=FALSE, cache=FALSE--------------------------------------------------------------------------------------------------------------------------------
# # Protein data from YeastGFP study
# protein <- c("YDL195W", "YHR098C", "YLR208W", "YNL049C", "YPL085W")
# abundance <- c(1840, 12200, 21400, 1720, 358)
# # Find protein indices in CHNOSZ database
# ip <- match(protein, thermo()$protein$protein)
# # Calculate protein lengths
# pl <- protein.length(ip)
# # Scale experimental abundances to unit total abundance of residues
# scaled_logabundance <- unitize(log10(abundance), pl)
# 
# # Assign colors from low to high oxidation state (red to blue)
# col_in <- c("darkred", "red", "darkgray", "blue", "darkblue")
# aa <- pinfo(ip)
# Zc_values <- canprot::Zc(aa)
# col <- col_in[rank(Zc_values)]
# 
# # Setup chemical system with H+ and e-
# basis("CHNOSe")
# # Calculate affinities as a function of redox potential
# a <- affinity(Eh = c(-0.6, 0), iprotein = ip)
# 
# # Calculate equilibrium distributions without normalization
# e <- equilibrate(a, loga.balance = 0)
# # Plot results
# par(mfrow = c(1, 4))
# diagram(e, ylim = c(-5, -2), col = col, lwd = 2, main = "Without normalization")
# 
# # With normalization
# e_norm <- equilibrate(a, normalize = TRUE, loga.balance = 0)
# diagram(e_norm, ylim = c(-5, -2.5), col = col, lwd = 2, main = "With normalization")
# 
# # Calculate and plot mean absolute error as a function of Eh
# predicted_logabundances <- sapply(e_norm$loga.equil, c)
# MAE <- apply(abs(t(predicted_logabundances) - scaled_logabundance), 2, mean)
# par(xaxs = "r", yaxs = "r")
# plot(a$vals$Eh, MAE)
# imin <- which.min(MAE)
# abline(v = a$vals$Eh[imin], lty = 2)
# title("Optimized model Eh")
# 
# # Show calculated and experimental abundances
# optimized_logabundance <- predicted_logabundances[imin, ]
# plot(scaled_logabundance, optimized_logabundance, pch = 19, col = col)
# # Show 1:1 line
# lims <- range(par("usr"))
# lines(c(lims[1], lims[2]), c(lims[1], lims[2]), lty = 2)
# title("Predicted vs experimental abundance")
# # Add legends
# legend("topleft", paste("MAE =", round(MAE[imin], 2)), bty = "n")
# legend("bottomright", c("High Zc", NA, NA, NA, "Low Zc"), pch = 19, col = rev(col_in))

## ----Cas_Zc, echo=FALSE, fig.cap = "Carbon oxidation state and size of CRISPR-Cas effector modules", cache=TRUE, small_mar = TRUE-------------------------------------------------
# Read data table
file <- system.file("extdata/protein/Cas/Cas_uniprot.csv", package = "CHNOSZ")
dat <- read.csv(file)
# Use UniProt ID as the file name
ID <- dat$UniProt
# In case UniProt ID is missing, use alternate ID
ID[ID == ""] <- dat$Protein[ID == ""]
# Store ID in data frame
dat$ID <- ID
# Remove missing IDs
dat <- subset(dat, ID != "")
# Keep proteins in effector complexes - the functional units that cut DNA
dat <- subset(dat, Effector)

# Read amino acid compositions
aafile <- system.file("extdata/protein/Cas/Cas_aa.csv", package = "CHNOSZ")
aa <- read.csv(aafile)

# Loop over evolutionary subtypes (I-A, I-B etc.)
subtypes <- unique(dat$Subtype)
effector_aa_list <- lapply(subtypes, function(subtype) {
  # Get the IDs for this evolutionary subtype
  idat <- dat$Subtype == subtype
  ID <- dat$ID[idat]
  # Combine amino acid compositions of all effector proteins in this subtype
  iaa <- aa$ref %in% ID
  all_aa <- aa[iaa, ]
  summed_aa <- canprot::sum_aa(all_aa)
  # Store gene names and IDs for reference
  summed_aa$protein <- paste(all_aa$protein, collapse = ",")
  summed_aa$ref <- paste(all_aa$ref, collapse = ",")
  # Label with the evolutionary subtype
  summed_aa$abbrv <- subtype
  # Return the combined amino acid composition for this evolutionary branch
  summed_aa
})
# Create data frame of amino acid compositions for each evolutionary branch
effector_aa <- do.call(rbind, effector_aa_list)

# Identify the six major evolutionary types (I to VI)
type_names <- c("I", "II", "III", "IV", "V", "VI")
# Map subtypes to their parent types
subtype_type <- sapply(strsplit(effector_aa$abbrv, "-"), "[", 1)
# Assign evolutionary classes (1 or 2) based on system architecture
class <- ifelse(subtype_type %in% c("I", "III", "IV"), 1, 2)

# Color scheme to distinguish evolutionary classes
col1 <- hcl.colors(5, "Peach")[1:3]
col2 <- hcl.colors(5, "Purp")[1:3]
type_col <- c(
  # Class 1 - multi-protein complexes
  I = col1[1], III = col1[2], IV = col1[3],
  # Class 2 - single-protein effectors
  II = col2[1], V = col2[2], VI = col2[3]
)
# Point sizes to show types in each class
type_cex <- c(
  I = 1, III = 1.5, IV = 2,
  II = 1, V = 1.5, VI = 2
)
# Apply colors and sizes to subtypes
subtype_col <- type_col[subtype_type]
subtype_cex <- type_cex[subtype_type]

# Calculate chemical features of effector complexes
Zc <- canprot::Zc(effector_aa)
naa <- protein.length(effector_aa)
# Plot chemical features
plot(naa, Zc, pch = class + 20, cex = subtype_cex, col = subtype_col, bg = subtype_col, xlab = "Number of amino acids", ylab = quote(italic(Z)[C]))
legend("topright", c("I               ", "III", "IV"), pch = 21, pt.cex = type_cex[1:3], col = col1, pt.bg = col1, title = "Class 1               ", cex = 0.95)
legend("topright", c("II", "V", "VI"), pch = 22, pt.cex = type_cex[4:6], col = col2, pt.bg = col2, title = "Class 2", bty = "n", cex = 0.95)

## ----Cas_Zc, eval=FALSE, cache=FALSE----------------------------------------------------------------------------------------------------------------------------------------------
# # Read data table
# file <- system.file("extdata/protein/Cas/Cas_uniprot.csv", package = "CHNOSZ")
# dat <- read.csv(file)
# # Use UniProt ID as the file name
# ID <- dat$UniProt
# # In case UniProt ID is missing, use alternate ID
# ID[ID == ""] <- dat$Protein[ID == ""]
# # Store ID in data frame
# dat$ID <- ID
# # Remove missing IDs
# dat <- subset(dat, ID != "")
# # Keep proteins in effector complexes - the functional units that cut DNA
# dat <- subset(dat, Effector)
# 
# # Read amino acid compositions
# aafile <- system.file("extdata/protein/Cas/Cas_aa.csv", package = "CHNOSZ")
# aa <- read.csv(aafile)
# 
# # Loop over evolutionary subtypes (I-A, I-B etc.)
# subtypes <- unique(dat$Subtype)
# effector_aa_list <- lapply(subtypes, function(subtype) {
#   # Get the IDs for this evolutionary subtype
#   idat <- dat$Subtype == subtype
#   ID <- dat$ID[idat]
#   # Combine amino acid compositions of all effector proteins in this subtype
#   iaa <- aa$ref %in% ID
#   all_aa <- aa[iaa, ]
#   summed_aa <- canprot::sum_aa(all_aa)
#   # Store gene names and IDs for reference
#   summed_aa$protein <- paste(all_aa$protein, collapse = ",")
#   summed_aa$ref <- paste(all_aa$ref, collapse = ",")
#   # Label with the evolutionary subtype
#   summed_aa$abbrv <- subtype
#   # Return the combined amino acid composition for this evolutionary branch
#   summed_aa
# })
# # Create data frame of amino acid compositions for each evolutionary branch
# effector_aa <- do.call(rbind, effector_aa_list)
# 
# # Identify the six major evolutionary types (I to VI)
# type_names <- c("I", "II", "III", "IV", "V", "VI")
# # Map subtypes to their parent types
# subtype_type <- sapply(strsplit(effector_aa$abbrv, "-"), "[", 1)
# # Assign evolutionary classes (1 or 2) based on system architecture
# class <- ifelse(subtype_type %in% c("I", "III", "IV"), 1, 2)
# 
# # Color scheme to distinguish evolutionary classes
# col1 <- hcl.colors(5, "Peach")[1:3]
# col2 <- hcl.colors(5, "Purp")[1:3]
# type_col <- c(
#   # Class 1 - multi-protein complexes
#   I = col1[1], III = col1[2], IV = col1[3],
#   # Class 2 - single-protein effectors
#   II = col2[1], V = col2[2], VI = col2[3]
# )
# # Point sizes to show types in each class
# type_cex <- c(
#   I = 1, III = 1.5, IV = 2,
#   II = 1, V = 1.5, VI = 2
# )
# # Apply colors and sizes to subtypes
# subtype_col <- type_col[subtype_type]
# subtype_cex <- type_cex[subtype_type]
# 
# # Calculate chemical features of effector complexes
# Zc <- canprot::Zc(effector_aa)
# naa <- protein.length(effector_aa)
# # Plot chemical features
# plot(naa, Zc, pch = class + 20, cex = subtype_cex, col = subtype_col, bg = subtype_col, xlab = "Number of amino acids", ylab = quote(italic(Z)[C]))
# legend("topright", c("I               ", "III", "IV"), pch = 21, pt.cex = type_cex[1:3], col = col1, pt.bg = col1, title = "Class 1               ", cex = 0.95)
# legend("topright", c("II", "V", "VI"), pch = 22, pt.cex = type_cex[4:6], col = col2, pt.bg = col2, title = "Class 2", bty = "n", cex = 0.95)

## ----Cas_stability, echo=FALSE, fig.cap = "Groupwise relative stabilities of effector modules in different types of CRISPR-Cas systems as a function of Eh and pH at 25 °C; dashed lines are water stability limits", cache=TRUE----
# Setup figure
par(mfrow = c(1, 2))

# Setup basis species to represent chemical variables
basis("QEC+")
swap.basis("O2", "e-")

# Load amino acid compositions normalized to one residue equivalent
iprotein <- add.protein(effector_aa, as.residue = TRUE)

# Calculate formation affinities for Class 1 effectors across environmental gradients
aout <- affinity(pH = c(0, 14), Eh = c(-0.8, 0.8), iprotein = iprotein[class == 1])
# Group proteins by evolutionary type for comparative analysis
class_1_groups <- sapply(c("I", "III", "IV"), `==`, subtype_type[class == 1], simplify = FALSE)
# Calculate average rank of formation affinity for each evolutionary group
arank <- rank.affinity(aout, groups = class_1_groups)
# Create stability diagram showing thermodynamically favored types
d <- diagram(arank, fill = col1)
water.lines(d, lty = 2)
title("Class 1", font.main = 1)

# Generate second diagram for Class 2
aout <- affinity(pH = c(0, 14), Eh = c(-0.8, 0.8), iprotein = iprotein[class == 2])
class_2_groups <- sapply(c("II", "V", "VI"), `==`, subtype_type[class == 2], simplify = FALSE)
arank <- rank.affinity(aout, groups = class_2_groups)
d <- diagram(arank, fill = col2)
water.lines(d, lty = 2)
title("Class 2", font.main = 1)

## ----Cas_stability, eval=FALSE, cache=FALSE---------------------------------------------------------------------------------------------------------------------------------------
# # Setup figure
# par(mfrow = c(1, 2))
# 
# # Setup basis species to represent chemical variables
# basis("QEC+")
# swap.basis("O2", "e-")
# 
# # Load amino acid compositions normalized to one residue equivalent
# iprotein <- add.protein(effector_aa, as.residue = TRUE)
# 
# # Calculate formation affinities for Class 1 effectors across environmental gradients
# aout <- affinity(pH = c(0, 14), Eh = c(-0.8, 0.8), iprotein = iprotein[class == 1])
# # Group proteins by evolutionary type for comparative analysis
# class_1_groups <- sapply(c("I", "III", "IV"), `==`, subtype_type[class == 1], simplify = FALSE)
# # Calculate average rank of formation affinity for each evolutionary group
# arank <- rank.affinity(aout, groups = class_1_groups)
# # Create stability diagram showing thermodynamically favored types
# d <- diagram(arank, fill = col1)
# water.lines(d, lty = 2)
# title("Class 1", font.main = 1)
# 
# # Generate second diagram for Class 2
# aout <- affinity(pH = c(0, 14), Eh = c(-0.8, 0.8), iprotein = iprotein[class == 2])
# class_2_groups <- sapply(c("II", "V", "VI"), `==`, subtype_type[class == 2], simplify = FALSE)
# arank <- rank.affinity(aout, groups = class_2_groups)
# d <- diagram(arank, fill = col2)
# water.lines(d, lty = 2)
# title("Class 2", font.main = 1)

## ----the_end----------------------------------------------------------------------------------------------------------------------------------------------------------------------
   ######    ##   ##    ##   ##    ######     #####  #####
 ##         ##===##    ## \\##   ##    ##     \\       //
 ######    ##   ##    ##   ##    ######    #####      #####

## ----sessionInfo, echo = FALSE----------------------------------------------------------------------------------------------------------------------------------------------------
## For finding what versions of packages are on R-Forge and winbuilder
#sessionInfo()

