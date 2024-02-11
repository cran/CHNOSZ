## ----setup, include=FALSE-----------------------------------------------------
options(width = 80)
## Use pngquant to optimize PNG images
library(knitr)
knit_hooks$set(pngquant = hook_pngquant)
pngquant <- "--speed=1 --quality=0-25"
if (!nzchar(Sys.which("pngquant"))) pngquant <- NULL

## Resolution settings
# Change this to TRUE to make high-resolution plots
# (default is FALSE to save time in CRAN checks)
hires <- nzchar(Sys.getenv("CHNOSZ_BUILD_LARGE_VIGNETTES"))
res1.lo <- 150
res1.hi <- 256
res1 <- if(hires) res1.hi else res1.lo
res2.lo <- 200
res2.hi <- 400
res2 <- if(hires) res2.hi else res2.lo

## logK with a thin space 20200627
logK <- "log&thinsp;<i>K</i>"

## Set dpi 20231129
knitr::opts_chunk$set(
  dpi = if(nzchar(Sys.getenv("CHNOSZ_BUILD_LARGE_VIGNETTES"))) 100 else 72
)

## ----CHNOSZ_reset, include=FALSE----------------------------------------------
library(CHNOSZ)
reset()

## ----res, results = "asis", echo = FALSE--------------------------------------
cat("```")
cat("\n")
cat(paste0(if(hires) "# " else "", "res1 <- ", res1.lo))
cat("\n")
cat(paste0(if(hires) "" else "# ", "res1 <- ", res1.hi))
cat("\n")
cat(paste0(if(hires) "# " else "", "res2 <- ", res2.lo))
cat("\n")
cat(paste0(if(hires) "" else "# ", "res2 <- ", res2.hi))
cat("\n")
cat("```")

## ----mash, echo = 1:8, eval = FALSE-------------------------------------------
#  par(mfrow = c(1, 2))
#  basis("CHNOS+")
#  species(c("CH4", "CO2", "HCO3-", "CO3-2"))
#  aC <- affinity(pH = c(0, 14), O2 = c(-75, -60))
#  dC <- diagram(aC, dx = c(0, 1, 0, 0), dy = c(0, 1, 0, 0))
#  species(c("H2S", "HS-", "HSO4-", "SO4-2"))
#  aS <- affinity(aC)  # argument recall
#  dS <- diagram(aS, add = TRUE, col = 4, col.names = 4, dx = c(0, -0.5, 0, 0))
#  aCS <- mash(dC, dS)
#  srt <- c(0, 0, 90, 0, 0, 0, 90, 0, 0, 0)
#  cex.names <- c(1, 1, 0.8, 1, 1, 1, 1, 1, 1, 1)
#  dy <- c(0, 0, 0, -0.2, 0, 0, 0, 0, 0, 0)
#  diagram(aCS, srt = srt, cex.names = cex.names, dy = dy)
#  legend("topright", legend = lTP(25, 1), bty = "n")

## ----mash, echo = 9:14,  results = "hide", message = FALSE, fig.width = 10, fig.height = 5, out.width = "100%"----
par(mfrow = c(1, 2))
basis("CHNOS+")
species(c("CH4", "CO2", "HCO3-", "CO3-2"))
aC <- affinity(pH = c(0, 14), O2 = c(-75, -60))
dC <- diagram(aC, dx = c(0, 1, 0, 0), dy = c(0, 1, 0, 0))
species(c("H2S", "HS-", "HSO4-", "SO4-2"))
aS <- affinity(aC)  # argument recall
dS <- diagram(aS, add = TRUE, col = 4, col.names = 4, dx = c(0, -0.5, 0, 0))
aCS <- mash(dC, dS)
srt <- c(0, 0, 90, 0, 0, 0, 90, 0, 0, 0)
cex.names <- c(1, 1, 0.8, 1, 1, 1, 1, 1, 1, 1)
dy <- c(0, 0, 0, -0.2, 0, 0, 0, 0, 0, 0)
diagram(aCS, srt = srt, cex.names = cex.names, dy = dy)
legend("topright", legend = lTP(25, 1), bty = "n")

## ----materials, message = FALSE, results = "hide"-----------------------------
## Formation energies (eV / atom) for solids from Materials API, e.g.
# from pymatgen import MPRester
# m = MPRester("USER_API_KEY")
# m.query(criteria={"task_id": "mp-1279742"}, properties=["formation_energy_per_atom"])
# mp-13, mp-1279742, mp-19306, mp-19770
#Fe.cr <- c(Fe = 0, FeO = -1.72803, Fe3O4 = -1.85868, Fe2O3 = -1.90736)  # 20201109
Fe.cr <- c(Fe = 0, FeO = -1.72768, Fe3O4 = -1.85838, Fe2O3 = -1.90708)  # 20210219
# mp-146, mp-18937, mp-1275946, mp-19094, mp-756395, mp-25279
#V.cr <- c(V = 0, V2O3 = -2.52849, V3O5 = -2.52574, VO2 = -2.48496, V3O7 = -2.32836, V2O5 = -2.29524)  # 20201109
V.cr <- c(V = 0, V2O3 = -2.52787, V3O5 = -2.52516, VO2 = -2.48447, V3O7 = -2.32789, V2O5 = -2.29480)  # 20210219

# Convert formation energies from eV / atom to eV / molecule
natom.Fe <- sapply(makeup(names(Fe.cr)), sum)
Fe.cr <- Fe.cr * natom.Fe
natom.V <- sapply(makeup(names(V.cr)), sum)
V.cr <- V.cr * natom.V

# Convert formation energies from eV / molecule to J / mol
eVtoJ <- function(eV) eV * 1.602176634e-19 * 6.02214076e23
Fe.cr <- eVtoJ(Fe.cr)
V.cr <- eVtoJ(V.cr)

# Gibbs energies of formation (J / mol) for aqueous species
# Most are from Wagman et al., 1982
Fe.aq <- 1000 * c("Fe+2" = -78.90, "Fe+3" = -4.7, "FeO2-2" = -295.3,
  "FeOH+" = -277.4, "FeOH+2" = -229.41, "HFeO2-" = -377.7,
  "Fe(OH)2+" = -438.0, "Fe(OH)3" = -659.3,
  "FeO2-" = -368.2, # SSWS97
  "FeO4-2" = -322.2 # Mis73
)

V.aq <- 1000 * c("VO+2" = -446.4, "VO2+" = -587.0, "VO3-" = -783.6, "VO4-3" = -899.0,
  "V2O7-4" = -1719, "HVO4" = -745.1, "HVO4-2" = -974.8,
  "VOH2O2+3" = -523.4, "VO2H2O2+" = -746.3, "V2HO7-3" = 1792.2, "V2H3O7-" = -1863.8,
  "HV10O28-5" = -7702, "H2V10O28-4" = -7723
)

# Gibbs energies of formation (J / mol) for solids from Wagman et al., 1982
Fe3O4 <- 1000 * -1015.4 # magnetite
V3O5 <- 1000 * -1803

# Calculate correction for difference between reference and DFT energies (Persson et al., 2012)
Fe.corr <- (Fe.cr["Fe3O4"] - Fe3O4) / 3
V.corr <- (V.cr["V3O5"] - V3O5) / 2

# Apply correction to standard Gibbs energies of aqueous species (Persson et al., 2012)
nFe <- sapply(makeup(names(Fe.aq)), "[", "Fe")
Fe.aq <- Fe.aq + nFe * Fe.corr
nV <- sapply(makeup(names(V.aq)), "[", "V")
V.aq <- V.aq + nV * V.corr

# Add energies to OBIGT
# This function modifies OBIGT and returns the species indices of the affected species
modfun <- function(x, state, model = NULL) sapply(seq_along(x), function(i) {
  # We explicitly set the units to Joules (this is the default as of CHNOSZ 2.0.0)
  if(is.null(model)) mod.OBIGT(names(x)[i], formula = names(x)[i], state = state, E_units = "J", G = x[i])
  else mod.OBIGT(names(x)[i], formula = names(x)[i], state = state, model = model, E_units = "J", G = x[i])
})
# We need model = "CGL" to override the Berman model for some minerals 20220919
iFe.cr <- modfun(Fe.cr, "cr", model = "CGL")
iFe.aq <- modfun(Fe.aq, "aq")
iV.cr <- modfun(V.cr, "cr", model = "CGL")
iV.aq <- modfun(V.aq, "aq")

# Formation energies (eV / atom) for bimetallic solids from Materials API
# mp-1335, mp-1079399, mp-866134, mp-558525, mp-504509 (triclinic FeVO4)
#FeV.cr <- c(FeV = -0.12928, FeV3 = -0.17128, Fe3V = -0.12854, Fe2V4O13 = -2.23833, FeVO4 = -1.75611)  # 20201109
FeV.cr <- c(FeV = -0.12815, FeV3 = -0.17114, Fe3V = -0.12832, Fe2V4O13 = -2.23967, FeVO4 = -1.75573)  # 20210219
# Convert energies and add to OBIGT
natom.FeV <- sapply(makeup(names(FeV.cr)), sum)
FeV.cr <- FeV.cr * natom.FeV
FeV.cr <- eVtoJ(FeV.cr)
iFeV.cr <- modfun(FeV.cr, "cr")

## ----mixing1, eval = FALSE, echo = 1:14---------------------------------------
#  par(mfrow = c(1, 3))
#  loga.Fe <- -5
#  loga.V <- -5
#  # Fe-O-H diagram
#  basis(c("VO+2", "Fe+2", "H2O", "e-", "H+"))
#  species(c(iFe.aq, iFe.cr))
#  species(1:length(iFe.aq), loga.Fe)
#  aFe <- affinity(pH = c(4, 10, res1), Eh = c(-1.5, 0, res1))
#  dFe <- diagram(aFe, plot.it = FALSE)
#  # V-O-H diagram
#  species(c(iV.aq, iV.cr))
#  species(1:length(iV.aq), loga.V)
#  aV <- affinity(aFe)  # argument recall
#  dV <- diagram(aV, plot.it = FALSE)
#  
#  # Calculate affinities for bimetallic species
#  species(iFeV.cr)
#  aFeV <- affinity(aFe)  # argument recall
#  dFeV <- diagram(aFeV, plot.it = FALSE, bold = TRUE)
#  # 1:1 mixture (Fe:V)
#  a11 <- mix(dFe, dV, dFeV, c(1, 1))
#  # Adjust labels 20210219
#  iV2O3 <- info("V2O3")
#  iFeO <- info("FeO", "cr")
#  iFe3V <- info("Fe3V")
#  srt <- rep(0, nrow(a11$species))
#  srt[a11$species$ispecies == paste(iFeO, iV2O3, sep = ",")] <- 90
#  srt[a11$species$ispecies == paste(iV2O3, iFe3V, sep = ",")] <- -13
#  diagram(a11, min.area = 0.01, srt = srt)
#  title("Fe:V = 1:1")
#  label.figure(lTP(25, 1), xfrac = 0.12)
#  # 1:3 mixture
#  a13 <- mix(dFe, dV, dFeV, c(1, 3))
#  srt <- rep(0, nrow(a13$species))
#  srt[a13$species$ispecies == paste(iFeO, iV2O3, sep = ",")] <- 90
#  srt[a13$species$ispecies == paste(iV2O3, iFe3V, sep = ",")] <- -13
#  diagram(a13, min.area = 0.01, srt = srt)
#  title("Fe:V = 1:3")
#  # 1:5 mixture
#  a15 <- mix(dFe, dV, dFeV, c(1, 5))
#  iFeV3 <- info("FeV3")
#  srt <- rep(0, nrow(a15$species))
#  srt[a15$species$ispecies == paste(iFeO, iV2O3, sep = ",")] <- 90
#  srt[a15$species$ispecies == paste(iV2O3, iFe3V, sep = ",")] <- -13
#  srt[a15$species$ispecies == paste(iV2O3, iFeV3, sep = ",")] <- -13
#  diagram(a15, min.area = 0.01, srt = srt)
#  title("Fe:V = 1:5")

## ----mixing1, echo = 16:47, message = FALSE, results = "hide", fig.width = 9, fig.height = 3, out.width = "100%", out.extra='class="full-width"', pngquant = pngquant----
par(mfrow = c(1, 3))
loga.Fe <- -5
loga.V <- -5
# Fe-O-H diagram
basis(c("VO+2", "Fe+2", "H2O", "e-", "H+"))
species(c(iFe.aq, iFe.cr))
species(1:length(iFe.aq), loga.Fe)
aFe <- affinity(pH = c(4, 10, res1), Eh = c(-1.5, 0, res1))
dFe <- diagram(aFe, plot.it = FALSE)
# V-O-H diagram
species(c(iV.aq, iV.cr))
species(1:length(iV.aq), loga.V)
aV <- affinity(aFe)  # argument recall
dV <- diagram(aV, plot.it = FALSE)

# Calculate affinities for bimetallic species
species(iFeV.cr)
aFeV <- affinity(aFe)  # argument recall
dFeV <- diagram(aFeV, plot.it = FALSE, bold = TRUE)
# 1:1 mixture (Fe:V)
a11 <- mix(dFe, dV, dFeV, c(1, 1))
# Adjust labels 20210219
iV2O3 <- info("V2O3")
iFeO <- info("FeO", "cr")
iFe3V <- info("Fe3V")
srt <- rep(0, nrow(a11$species))
srt[a11$species$ispecies == paste(iFeO, iV2O3, sep = ",")] <- 90
srt[a11$species$ispecies == paste(iV2O3, iFe3V, sep = ",")] <- -13
diagram(a11, min.area = 0.01, srt = srt)
title("Fe:V = 1:1")
label.figure(lTP(25, 1), xfrac = 0.12)
# 1:3 mixture
a13 <- mix(dFe, dV, dFeV, c(1, 3))
srt <- rep(0, nrow(a13$species))
srt[a13$species$ispecies == paste(iFeO, iV2O3, sep = ",")] <- 90
srt[a13$species$ispecies == paste(iV2O3, iFe3V, sep = ",")] <- -13
diagram(a13, min.area = 0.01, srt = srt)
title("Fe:V = 1:3")
# 1:5 mixture
a15 <- mix(dFe, dV, dFeV, c(1, 5))
iFeV3 <- info("FeV3")
srt <- rep(0, nrow(a15$species))
srt[a15$species$ispecies == paste(iFeO, iV2O3, sep = ",")] <- 90
srt[a15$species$ispecies == paste(iV2O3, iFe3V, sep = ",")] <- -13
srt[a15$species$ispecies == paste(iV2O3, iFeV3, sep = ",")] <- -13
diagram(a15, min.area = 0.01, srt = srt)
title("Fe:V = 1:5")

## ----FeVO4, eval = FALSE, echo = 1:29-----------------------------------------
#  layout(t(matrix(1:3)), widths = c(1, 1, 0.2))
#  par(cex = 1)
#  # Fe-bearing species
#  basis(c("VO+2", "Fe+2", "H2O", "e-", "H+"))
#  species(c(iFe.aq, iFe.cr))$name
#  species(1:length(iFe.aq), loga.Fe)
#  aFe <- affinity(pH = c(0, 14, res2), Eh = c(-1.5, 2, res2))
#  dFe <- diagram(aFe, plot.it = FALSE)
#  # V-bearing species
#  species(c(iV.aq, iV.cr))$name
#  species(1:length(iV.aq), loga.V)
#  aV <- affinity(aFe)  # argument recall
#  dV <- diagram(aV, plot.it = FALSE)
#  # Bimetallic species
#  species(iFeV.cr)
#  aFeV <- affinity(aFe)  # argument recall
#  dFeV <- diagram(aFeV, plot.it = FALSE, bold = TRUE)
#  # 1:1 mixture (Fe:V)
#  a11 <- mix(dFe, dV, dFeV, c(1, 1))
#  # Adjust labels 20210219
#  iV2O3 <- info("V2O3")
#  iFe3V <- info("Fe3V")
#  iVO4m3 <- info("VO4-3")
#  iFe2O3 <- info("Fe2O3")
#  srt <- rep(0, nrow(a11$species))
#  srt[a11$species$ispecies == paste(iV2O3, iFe3V, sep = ",")] <- -13
#  srt[a11$species$ispecies == paste(iFe2O3, iVO4m3, sep = ",")] <- 90
#  d11 <- diagram(a11, min.area = 0.01, srt = srt)
#  water.lines(d11, col = "orangered")
#  
#  # Calculate affinity of FeVO4
#  species("FeVO4")
#  aFeVO4 <- affinity(aFe)  # argument recall
#  # Calculate difference from stable species
#  aFeVO4_vs_stable <- aFeVO4$values[[1]] - d11$predominant.values
#  # Overlay lines from diagram on color map
#  diagram(a11, fill = NA, names = FALSE, limit.water = FALSE)
#  opar <- par(usr = c(0, 1, 0, 1))
#  col <- rev(topo.colors(128)) # No hcl.colors() in R < 3.6.0
#  if(getRversion() >= "3.6.0") col <- rev(hcl.colors(128, palette = "YlGnBu", alpha = 0.8))
#  image(aFeVO4_vs_stable, col = col, add = TRUE)
#  par(opar)
#  diagram(a11, fill = NA, add = TRUE, names = FALSE)
#  water.lines(d11, col = "orangered")
#  thermo.axis()
#  
#  imax <- arrayInd(which.max(aFeVO4_vs_stable), dim(aFeVO4_vs_stable))
#  pH <- d11$vals$pH[imax[1]]
#  Eh <- d11$vals$Eh[imax[2]]
#  points(pH, Eh, pch = 10, cex = 2, lwd = 2, col = "gold")
#  stable <- d11$names[d11$predominant[imax]]
#  text(pH, Eh, stable, adj = c(0.5, -1), cex = 1.2, col = "gold")
#  
#  # Make color scale 20210228
#  par(mar = c(3, 0, 2.5, 2.7))
#  plot.new()
#  levels <- 1:length(col)
#  plot.window(xlim = c(0, 1), ylim = range(levels), xaxs = "i", yaxs = "i")
#  rect(0, levels[-length(levels)], 1, levels[-1L], col = rev(col), border = NA)
#  box()
#  # To get the limits, convert range of affinities to eV/atom
#  arange <- rev(range(aFeVO4_vs_stable))
#  # This gets us to J/mol
#  Jrange <- convert(arange, "G")
#  # And to eV/atom
#  eVrange <- Jrange / 1.602176634e-19 / 6.02214076e23 / 6
#  ylim <- formatC(eVrange, digits = 3, format = "f")
#  axis(4, at = range(levels), labels = ylim)
#  mtext(quote(Delta*italic(G)[pbx]*", eV/atom"), side = 4, las = 0, line = 1)

## ----FeVO4, echo = 31:44, message = FALSE, results = "hide", fig.width = 11, fig.height = 5, out.width = "100%", pngquant = pngquant----
layout(t(matrix(1:3)), widths = c(1, 1, 0.2))
par(cex = 1)
# Fe-bearing species
basis(c("VO+2", "Fe+2", "H2O", "e-", "H+"))
species(c(iFe.aq, iFe.cr))$name
species(1:length(iFe.aq), loga.Fe)
aFe <- affinity(pH = c(0, 14, res2), Eh = c(-1.5, 2, res2))
dFe <- diagram(aFe, plot.it = FALSE)
# V-bearing species
species(c(iV.aq, iV.cr))$name
species(1:length(iV.aq), loga.V)
aV <- affinity(aFe)  # argument recall
dV <- diagram(aV, plot.it = FALSE)
# Bimetallic species
species(iFeV.cr)
aFeV <- affinity(aFe)  # argument recall
dFeV <- diagram(aFeV, plot.it = FALSE, bold = TRUE)
# 1:1 mixture (Fe:V)
a11 <- mix(dFe, dV, dFeV, c(1, 1))
# Adjust labels 20210219
iV2O3 <- info("V2O3")
iFe3V <- info("Fe3V")
iVO4m3 <- info("VO4-3")
iFe2O3 <- info("Fe2O3")
srt <- rep(0, nrow(a11$species))
srt[a11$species$ispecies == paste(iV2O3, iFe3V, sep = ",")] <- -13
srt[a11$species$ispecies == paste(iFe2O3, iVO4m3, sep = ",")] <- 90
d11 <- diagram(a11, min.area = 0.01, srt = srt)
water.lines(d11, col = "orangered")

# Calculate affinity of FeVO4
species("FeVO4")
aFeVO4 <- affinity(aFe)  # argument recall
# Calculate difference from stable species
aFeVO4_vs_stable <- aFeVO4$values[[1]] - d11$predominant.values
# Overlay lines from diagram on color map
diagram(a11, fill = NA, names = FALSE, limit.water = FALSE)
opar <- par(usr = c(0, 1, 0, 1))
col <- rev(topo.colors(128)) # No hcl.colors() in R < 3.6.0
if(getRversion() >= "3.6.0") col <- rev(hcl.colors(128, palette = "YlGnBu", alpha = 0.8))
image(aFeVO4_vs_stable, col = col, add = TRUE)
par(opar)
diagram(a11, fill = NA, add = TRUE, names = FALSE)
water.lines(d11, col = "orangered")
thermo.axis()

imax <- arrayInd(which.max(aFeVO4_vs_stable), dim(aFeVO4_vs_stable))
pH <- d11$vals$pH[imax[1]]
Eh <- d11$vals$Eh[imax[2]]
points(pH, Eh, pch = 10, cex = 2, lwd = 2, col = "gold")
stable <- d11$names[d11$predominant[imax]]
text(pH, Eh, stable, adj = c(0.5, -1), cex = 1.2, col = "gold")

# Make color scale 20210228
par(mar = c(3, 0, 2.5, 2.7))
plot.new()
levels <- 1:length(col)
plot.window(xlim = c(0, 1), ylim = range(levels), xaxs = "i", yaxs = "i")
rect(0, levels[-length(levels)], 1, levels[-1L], col = rev(col), border = NA)
box()
# To get the limits, convert range of affinities to eV/atom
arange <- rev(range(aFeVO4_vs_stable))
# This gets us to J/mol
Jrange <- convert(arange, "G")
# And to eV/atom
eVrange <- Jrange / 1.602176634e-19 / 6.02214076e23 / 6
ylim <- formatC(eVrange, digits = 3, format = "f")
axis(4, at = range(levels), labels = ylim)
mtext(quote(Delta*italic(G)[pbx]*", eV/atom"), side = 4, las = 0, line = 1)

## ----Gpbx_min, echo = 2:8, message = FALSE, fig.keep = "none"-----------------
plot(1:10) # so we can run "points" in this chunk
imax <- arrayInd(which.max(aFeVO4_vs_stable), dim(aFeVO4_vs_stable))
pH <- d11$vals$pH[imax[1]]
Eh <- d11$vals$Eh[imax[2]]
points(pH, Eh, pch = 10, cex = 2, lwd = 2, col = "gold")
stable <- d11$names[d11$predominant[imax]]
text(pH, Eh, stable, adj = c(0.3, 2), cex = 1.2, col = "gold")
(Apbx <- range(aFeVO4_vs_stable[d11$predominant == d11$predominant[imax]]))

## ----hull, message = FALSE----------------------------------------------------
b <- basis(c("Fe2O3", "Fe2V4O13", "O2"))
J_mol <- subcrt("FeVO4", 1, T = 25)$out$G
stopifnot(all.equal(rep(convert(J_mol, "logK"), 2), Apbx))
eV_mol <- J_mol / 1.602176634e-19
eV_atom <- eV_mol / 6.02214076e23 / 6
round(eV_atom, 3)
stopifnot(round(eV_atom, 3) == 0.415)

## ----reset, message = FALSE---------------------------------------------------
reset()

## ----stack1_1, results = "hide", message = FALSE------------------------------
logaH2S <- -2
T <- 200
pH <- c(0, 14, res2)
O2 <- c(-48, -33, res2)
basis(c("Cu+", "pyrite", "H2S", "oxygen", "H2O", "H+"))
basis("H2S", logaH2S)
S.aq <- c("H2S", "HS-", "HSO4-", "SO4-2")
Fe.cr <- c("pyrite", "pyrrhotite", "magnetite", "hematite")
Fe.abbrv <- c("Py", "Po", "Mag", "Hem")

## ----stack1_2, eval = FALSE, echo = 1:4---------------------------------------
#  species(Fe.cr)
#  mFe <- mosaic(S.aq, pH = pH, O2 = O2, T = T)
#  diagram(mFe$A.bases, lty = 2, col = 4, col.names = 4, italic = TRUE, dx = c(0, 1, 0, 0), dy = c(-1.5, 0, 1, 0))
#  dFe <- diagram(mFe$A.species, add = TRUE, lwd = 2, names = Fe.abbrv, dx = c(0, 0.5, 0, 0), dy = c(-1, 0, 0.5, 0))
#  FeCu.cr <- c("chalcopyrite", "bornite")
#  Cu.cr <- c("copper", "cuprite", "tenorite", "chalcocite", "covellite")
#  FeCu.abbrv <- c("Ccp", "Bn", "Cu", "Cpr", "Tnr", "Cct", "Cv")
#  species(c(FeCu.cr, Cu.cr))
#  mFeCu <- mosaic(list(S.aq, Fe.cr), pH = pH, O2 = O2,
#                T = T, stable = list(NULL, dFe$predominant))
#  col <- c("#FF8C00", rep(2, 6))
#  lwd <- c(2, 1, 1, 1, 1, 1, 1)
#  dy = c(0, 0, 0, 0, 0, 1, 0)
#  diagram(mFeCu$A.species, add = TRUE, col = col, lwd = lwd, col.names = col, bold = TRUE, names = FeCu.abbrv, dy = dy)
#  TPS <- c(describe.property(c("T", "P"), c(T, "Psat")), expression(sum(S) == 0.01*m))
#  legend("topright", TPS, bty = "n")
#  title("Cu-Fe-S-O-H (minerals only)", font.main = 1)

## ----stack1_2, echo=5:17, results = "hide", message = FALSE, fig.width = 6, fig.height = 5, out.width = "75%", fig.align = "center", pngquant = pngquant----
species(Fe.cr)
mFe <- mosaic(S.aq, pH = pH, O2 = O2, T = T)
diagram(mFe$A.bases, lty = 2, col = 4, col.names = 4, italic = TRUE, dx = c(0, 1, 0, 0), dy = c(-1.5, 0, 1, 0))
dFe <- diagram(mFe$A.species, add = TRUE, lwd = 2, names = Fe.abbrv, dx = c(0, 0.5, 0, 0), dy = c(-1, 0, 0.5, 0))
FeCu.cr <- c("chalcopyrite", "bornite")
Cu.cr <- c("copper", "cuprite", "tenorite", "chalcocite", "covellite")
FeCu.abbrv <- c("Ccp", "Bn", "Cu", "Cpr", "Tnr", "Cct", "Cv")
species(c(FeCu.cr, Cu.cr))
mFeCu <- mosaic(list(S.aq, Fe.cr), pH = pH, O2 = O2,
              T = T, stable = list(NULL, dFe$predominant))
col <- c("#FF8C00", rep(2, 6))
lwd <- c(2, 1, 1, 1, 1, 1, 1)
dy = c(0, 0, 0, 0, 0, 1, 0)
diagram(mFeCu$A.species, add = TRUE, col = col, lwd = lwd, col.names = col, bold = TRUE, names = FeCu.abbrv, dy = dy)
TPS <- c(describe.property(c("T", "P"), c(T, "Psat")), expression(sum(S) == 0.01*m))
legend("topright", TPS, bty = "n")
title("Cu-Fe-S-O-H (minerals only)", font.main = 1)

## ----stack2, eval = FALSE-----------------------------------------------------
#  # Define system
#  pH <- c(0, 14, res2)
#  O2 <- c(-48, -33, res2)
#  T <- 200
#  logmS <- -2
#  m_NaCl <- 0.1
#  logm_aq <- -6 # for both Fe- and Cu-bearing aq species
#  # Basis species
#  S.aq <- c("H2S", "HS-", "HSO4-", "SO4-2")
#  # Minerals
#  Fe.cr <- c("pyrite", "pyrrhotite", "magnetite", "hematite")
#  Fe.abbrv <- c("Py", "Po", "Mag", "Hem")
#  FeCu.cr <- c("chalcopyrite", "bornite")
#  Cu.cr <- c("copper", "cuprite", "tenorite", "chalcocite", "covellite")
#  FeCu.abbrv <- c("Ccp", "Bn", "Cu", "Cpr", "Tnr", "Cct", "Cv")
#  # Aqueous species
#  iFe.aq <- retrieve("Fe", c("S", "O", "H", "Cl"), "aq")
#  Fe.aq <- info(iFe.aq)$name
#  iCu.aq <- retrieve("Cu", c("S", "O", "H", "Cl"), "aq")
#  Cu.aq <- info(iCu.aq)$name
#  # Expressions for making the legend
#  TPexpr <- describe.property(c("T", "P"), c(T, "Psat"))
#  Sexpr <- as.expression(bquote(sum(S) == .(10^logmS)*m))
#  NaClexpr <- as.expression(bquote(NaCl == .(m_NaCl)*m))
#  aqexpr <- as.expression(bquote("("*aq*")"[italic(i)] == 10^.(logm_aq)*m))
#  
#  # Setup basis species
#  basis(c("Cu+", "pyrite", "H2S", "oxygen", "H2O", "H+", "Cl-"))
#  basis("H2S", logmS)
#  nacl <- NaCl(m_tot = m_NaCl, T = T, P = "Psat")
#  basis("Cl-", log10(nacl$m_Cl))
#  # Fe-bearing minerals
#  species(Fe.cr)
#  # Add aqueous species 20210220
#  species(iFe.aq, logm_aq, add = TRUE)
#  mFe <- mosaic(S.aq, pH = pH, O2 = O2, T = T, IS = nacl$IS)
#  # Start plot with just the fields for transparency effect
#  dFe <- diagram(mFe$A.species, lwd = 0, names = FALSE)
#  
#  # Cu-bearing minerals
#  species(c(FeCu.cr, Cu.cr))
#  # Add aqueous species 20210220
#  species(iCu.aq, logm_aq, add = TRUE)
#  
#  ## Mosaic with all Fe species as basis species
#  #mFeCu <- mosaic(list(S.aq, c(Fe.cr, Fe.aq)), pH = pH, O2 = O2, T = T, stable = list(NULL, dFe$predominant))
#  # Use only predominant Fe species as basis species (to speed up calculation) 20210224
#  predom <- dFe$predominant
#  ipredom <- sort(unique(as.numeric(predom)))
#  for(i in seq_along(ipredom)) predom[dFe$predominant == ipredom[i]] <- i
#  Fe.predom <- c(Fe.cr, Fe.aq)[ipredom]
#  # Use loga_aq argument to control the activity of aqueous species in mosaic calculation 20220722
#  # c(NA, logm_aq) means to use:
#  #   basis()'s value for logact of aqueous S species
#  #   logm_aq for logact of aqueous Fe species
#  mFeCu <- mosaic(list(S.aq, Fe.predom), pH = pH, O2 = O2, T = T, stable = list(NULL, predom), IS = nacl$IS, loga_aq = c(NA, logm_aq))
#  
#  # Adjust labels
#  bold <- c(rep(TRUE, length(FeCu.abbrv)), rep(FALSE, length(Cu.aq)))
#  names <- c(FeCu.abbrv, Cu.aq)
#  srt <- dx <- dy <- rep(0, length(names))
#  cex <- rep(1, length(names))
#  dx[names == "Cu"] <- -1.5
#  dx[names == "Bn"] <- 1.4
#  dx[names == "CuHS"] <- 1
#  dx[names == "Cu+"] <- -0.5
#  dy[names == "Cu"] <- 3
#  dx[names == "Cct"] <- -2
#  dy[names == "Cct"] <- 4
#  dy[names == "CuHS"] <- 1
#  dy[names == "Bn"] <- -0.9
#  dx[names == "CuCl2-"] <- -1
#  dy[names == "CuCl2-"] <- 2
#  cex[names == "Bn"] <- 0.8
#  srt[names == "Bn"] <- 85
#  # Highlight Ccp field
#  col.names <- col <- rep(2, nrow(mFeCu$A.species$species))
#  col[1] <- "#FF8C00"
#  col.names[1] <- "#FF8C00"
#  lwd <- rep(1, nrow(mFeCu$A.species$species))
#  lwd[1] <- 2
#  diagram(mFeCu$A.species, add = TRUE, lwd = lwd, col = col, col.names = col.names, names = names, bold = bold, dx = dx, dy = dy, cex.names = cex, srt = srt)
#  # Add second Cu label
#  text(12.3, -47, "Cu", col = 2, font = 2)
#  
#  # Plot the Fe-system lines and names "on top" so they are not covered by fill colors
#  diagram(mFe$A.bases, add = TRUE, lty = 2, col = 4, names = FALSE, fill = NA)
#  bold <- c(rep(TRUE, length(Fe.abbrv)), rep(FALSE, length(Fe.aq)))
#  names <- c(Fe.abbrv, Fe.aq)
#  srt <- dx <- dy <- rep(0, length(names))
#  cex <- rep(1, length(names))
#  dy[names == "Hem"] <- 0.5
#  dy[names == "Mag"] <- -0.8
#  dx[names == "Hem"] <- -0.5
#  dx[names == "Mag"] <- 0.25
#  dx[names == "Fe+2"] <- 0.5
#  dx[names == "FeO2-"] <- 1
#  dy[names == "FeO2-"] <- -3
#  dx[names == "HFeO2-"] <- -0.5
#  dy[names == "HFeO2-"] <- 1
#  srt[names == "Mag"] <- 83
#  srt[names == "FeSO4"] <- 90
#  diagram(mFe$A.species, add = TRUE, lwd = 2, names = names, bold = bold, dx = dx, dy = dy, cex.names = cex, srt = srt, fill = NA)
#  legend("topright", c(TPexpr, Sexpr, NaClexpr, aqexpr), bty = "n")
#  title("Cu-Fe-S-O-H-Cl (minerals and aqueous species)", font.main = 1)
#  
#  # Restore default OBIGT database
#  OBIGT()

## ----stack2, echo = FALSE, results = "hide", message = FALSE, fig.width = 6, fig.height = 5, out.width = "75%", fig.align = "center", pngquant = pngquant----
# Define system
pH <- c(0, 14, res2)
O2 <- c(-48, -33, res2)
T <- 200
logmS <- -2
m_NaCl <- 0.1
logm_aq <- -6 # for both Fe- and Cu-bearing aq species
# Basis species
S.aq <- c("H2S", "HS-", "HSO4-", "SO4-2")
# Minerals
Fe.cr <- c("pyrite", "pyrrhotite", "magnetite", "hematite")
Fe.abbrv <- c("Py", "Po", "Mag", "Hem")
FeCu.cr <- c("chalcopyrite", "bornite")
Cu.cr <- c("copper", "cuprite", "tenorite", "chalcocite", "covellite")
FeCu.abbrv <- c("Ccp", "Bn", "Cu", "Cpr", "Tnr", "Cct", "Cv")
# Aqueous species
iFe.aq <- retrieve("Fe", c("S", "O", "H", "Cl"), "aq")
Fe.aq <- info(iFe.aq)$name
iCu.aq <- retrieve("Cu", c("S", "O", "H", "Cl"), "aq")
Cu.aq <- info(iCu.aq)$name
# Expressions for making the legend
TPexpr <- describe.property(c("T", "P"), c(T, "Psat"))
Sexpr <- as.expression(bquote(sum(S) == .(10^logmS)*m))
NaClexpr <- as.expression(bquote(NaCl == .(m_NaCl)*m))
aqexpr <- as.expression(bquote("("*aq*")"[italic(i)] == 10^.(logm_aq)*m))

# Setup basis species
basis(c("Cu+", "pyrite", "H2S", "oxygen", "H2O", "H+", "Cl-"))
basis("H2S", logmS)
nacl <- NaCl(m_tot = m_NaCl, T = T, P = "Psat")
basis("Cl-", log10(nacl$m_Cl))
# Fe-bearing minerals
species(Fe.cr)
# Add aqueous species 20210220
species(iFe.aq, logm_aq, add = TRUE)
mFe <- mosaic(S.aq, pH = pH, O2 = O2, T = T, IS = nacl$IS)
# Start plot with just the fields for transparency effect
dFe <- diagram(mFe$A.species, lwd = 0, names = FALSE)

# Cu-bearing minerals
species(c(FeCu.cr, Cu.cr))
# Add aqueous species 20210220
species(iCu.aq, logm_aq, add = TRUE)

## Mosaic with all Fe species as basis species
#mFeCu <- mosaic(list(S.aq, c(Fe.cr, Fe.aq)), pH = pH, O2 = O2, T = T, stable = list(NULL, dFe$predominant))
# Use only predominant Fe species as basis species (to speed up calculation) 20210224
predom <- dFe$predominant
ipredom <- sort(unique(as.numeric(predom)))
for(i in seq_along(ipredom)) predom[dFe$predominant == ipredom[i]] <- i
Fe.predom <- c(Fe.cr, Fe.aq)[ipredom]
# Use loga_aq argument to control the activity of aqueous species in mosaic calculation 20220722
# c(NA, logm_aq) means to use:
#   basis()'s value for logact of aqueous S species
#   logm_aq for logact of aqueous Fe species
mFeCu <- mosaic(list(S.aq, Fe.predom), pH = pH, O2 = O2, T = T, stable = list(NULL, predom), IS = nacl$IS, loga_aq = c(NA, logm_aq))

# Adjust labels
bold <- c(rep(TRUE, length(FeCu.abbrv)), rep(FALSE, length(Cu.aq)))
names <- c(FeCu.abbrv, Cu.aq)
srt <- dx <- dy <- rep(0, length(names))
cex <- rep(1, length(names))
dx[names == "Cu"] <- -1.5
dx[names == "Bn"] <- 1.4
dx[names == "CuHS"] <- 1
dx[names == "Cu+"] <- -0.5
dy[names == "Cu"] <- 3
dx[names == "Cct"] <- -2
dy[names == "Cct"] <- 4
dy[names == "CuHS"] <- 1
dy[names == "Bn"] <- -0.9
dx[names == "CuCl2-"] <- -1
dy[names == "CuCl2-"] <- 2
cex[names == "Bn"] <- 0.8
srt[names == "Bn"] <- 85
# Highlight Ccp field
col.names <- col <- rep(2, nrow(mFeCu$A.species$species))
col[1] <- "#FF8C00"
col.names[1] <- "#FF8C00"
lwd <- rep(1, nrow(mFeCu$A.species$species))
lwd[1] <- 2
diagram(mFeCu$A.species, add = TRUE, lwd = lwd, col = col, col.names = col.names, names = names, bold = bold, dx = dx, dy = dy, cex.names = cex, srt = srt)
# Add second Cu label
text(12.3, -47, "Cu", col = 2, font = 2)

# Plot the Fe-system lines and names "on top" so they are not covered by fill colors
diagram(mFe$A.bases, add = TRUE, lty = 2, col = 4, names = FALSE, fill = NA)
bold <- c(rep(TRUE, length(Fe.abbrv)), rep(FALSE, length(Fe.aq)))
names <- c(Fe.abbrv, Fe.aq)
srt <- dx <- dy <- rep(0, length(names))
cex <- rep(1, length(names))
dy[names == "Hem"] <- 0.5
dy[names == "Mag"] <- -0.8
dx[names == "Hem"] <- -0.5
dx[names == "Mag"] <- 0.25
dx[names == "Fe+2"] <- 0.5
dx[names == "FeO2-"] <- 1
dy[names == "FeO2-"] <- -3
dx[names == "HFeO2-"] <- -0.5
dy[names == "HFeO2-"] <- 1
srt[names == "Mag"] <- 83
srt[names == "FeSO4"] <- 90
diagram(mFe$A.species, add = TRUE, lwd = 2, names = names, bold = bold, dx = dx, dy = dy, cex.names = cex, srt = srt, fill = NA)
legend("topright", c(TPexpr, Sexpr, NaClexpr, aqexpr), bty = "n")
title("Cu-Fe-S-O-H-Cl (minerals and aqueous species)", font.main = 1)

# Restore default OBIGT database
OBIGT()

## ----stack2.refs, message = FALSE---------------------------------------------
minerals <- list(Fe.cr = Fe.cr, Cu.cr = Cu.cr, FeCu.cr = FeCu.cr)
aqueous <- list(S.aq = S.aq, Fe.aq = Fe.aq, Cu.aq = Cu.aq)
allspecies <- c(minerals, aqueous)
iall <- lapply(allspecies, info)
allkeys <- lapply(iall, function(x) thermo.refs(x)$key)
allkeys

## ----stack2.cite, results = "asis"--------------------------------------------
allyears <- lapply(iall, function(x) thermo.refs(x)$year)
o <- order(unlist(allyears))
keys <- gsub("\\..*", "", unique(unlist(allkeys)[o]))
cat(paste(paste0("@", keys), collapse = "; "))

## ----mixing2, eval = FALSE----------------------------------------------------
#  par(mfrow = c(2, 2))
#  
#  logaH2S <- -2
#  T <- 200
#  pH <- c(0, 14)
#  O2 <- c(-60, -25)
#  basis(c("Cu+", "Fe+2", "H2S", "oxygen", "H2O", "H+"))
#  basis("H2S", logaH2S)
#  
#  S.aq <- c("H2S", "HS-", "HSO4-", "SO4-2")
#  Fe.cr <- c("pyrite", "pyrrhotite", "magnetite", "hematite")
#  Cu.cr <- c("copper", "cuprite", "tenorite", "chalcocite", "covellite")
#  FeCu.cr <- c("chalcopyrite", "bornite")
#  
#  species(Fe.cr)
#  mFe <- mosaic(S.aq, pH = pH, O2 = O2, T = T)
#  diagram(mFe$A.bases, lty = 2, col = 4, col.names = 4, italic = TRUE, dx = c(0, 1, 0, 0), dy = c(-1, 0, -2, 0))
#  names <- info(info(Fe.cr))$abbrv
#  dFe <- diagram(mFe$A.species, add = TRUE, names = names)
#  
#  species(Cu.cr)
#  mCu <- mosaic(S.aq, pH = pH, O2 = O2, T = T)
#  names <- info(info(Cu.cr))$abbrv
#  col.names <- rep(2, length(names))
#  col.names[1] <- 0
#  dCu <- diagram(mCu$A.species, add = TRUE, col = 2, col.names = col.names, names = names)
#  text(12, -55, "Cu", col = 2)
#  legend("topright", legend = lTP(T, "Psat"), bty = "n")
#  title(paste("Fe-S-O-H and Cu-S-O-H; Total S =", 10^logaH2S, "m"))
#  
#  species(FeCu.cr)
#  mFeCu <- mosaic(S.aq, pH = pH, O2 = O2, T = T)
#  names <- info(info(FeCu.cr))$abbrv
#  dFeCu <- diagram(mFeCu$A.species, plot.it = FALSE, names = names)
#  
#  fill <- function(a) {
#    ifelse(grepl("Ccp", a$species$name), "#FF8C0088",
#      ifelse(grepl("Bn", a$species$name), "#DC143C88", NA)
#  )}
#  srt <- function(a) ifelse(a$species$name %in% c("Mag+Bn", "Mag+Ccp"), 80, 0)
#  cex.names <- function(a) ifelse(a$species$name %in% c("Mag+Bn", "Mag+Ccp", "Mag+Cct"), 0.8, 1)
#  dx <- function(a) sapply(a$species$name, switch, "Mag+Bn" = 0.15, "Mag+Cct" = 0.5, 0)
#  dy <- function(a) sapply(a$species$name, switch, "Py+Bn" = -1, "Po+Bn" = -0.8, "Po+Cu" = -0.8, "Mag+Ccp" = -1, 0)
#  
#  a11 <- mix(dFe, dCu, dFeCu)
#  diagram(a11, fill = fill(a11), srt = srt(a11), min.area = 0.01, dx = dx(a11), dy = dy(a11), cex.names = cex.names(a11))
#  title("Fe:Cu = 1:1")
#  
#  a21 <- mix(dFe, dCu, dFeCu, c(2, 1))
#  diagram(a21, fill = fill(a21), srt = srt(a21), min.area = 0.01, dx = dx(a21), dy = dy(a21), cex.names = cex.names(a21))
#  title("Fe:Cu = 2:1")
#  
#  a12 <- mix(dFe, dCu, dFeCu, c(1, 2))
#  diagram(a12, fill = fill(a12), srt = srt(a12), min.area = 0.01, dx = dx(a12), dy = dy(a12), cex.names = cex.names(a12))
#  title("Fe:Cu = 1:2")

## ----mixing2, echo = FALSE, results = "hide", message = FALSE, fig.width = 8, fig.height = 6.5, out.width = "100%", pngquant = pngquant----
par(mfrow = c(2, 2))

logaH2S <- -2
T <- 200
pH <- c(0, 14)
O2 <- c(-60, -25)
basis(c("Cu+", "Fe+2", "H2S", "oxygen", "H2O", "H+"))
basis("H2S", logaH2S)

S.aq <- c("H2S", "HS-", "HSO4-", "SO4-2")
Fe.cr <- c("pyrite", "pyrrhotite", "magnetite", "hematite")
Cu.cr <- c("copper", "cuprite", "tenorite", "chalcocite", "covellite")
FeCu.cr <- c("chalcopyrite", "bornite")

species(Fe.cr)
mFe <- mosaic(S.aq, pH = pH, O2 = O2, T = T)
diagram(mFe$A.bases, lty = 2, col = 4, col.names = 4, italic = TRUE, dx = c(0, 1, 0, 0), dy = c(-1, 0, -2, 0))
names <- info(info(Fe.cr))$abbrv
dFe <- diagram(mFe$A.species, add = TRUE, names = names)

species(Cu.cr)
mCu <- mosaic(S.aq, pH = pH, O2 = O2, T = T)
names <- info(info(Cu.cr))$abbrv
col.names <- rep(2, length(names))
col.names[1] <- 0
dCu <- diagram(mCu$A.species, add = TRUE, col = 2, col.names = col.names, names = names)
text(12, -55, "Cu", col = 2)
legend("topright", legend = lTP(T, "Psat"), bty = "n")
title(paste("Fe-S-O-H and Cu-S-O-H; Total S =", 10^logaH2S, "m"))

species(FeCu.cr)
mFeCu <- mosaic(S.aq, pH = pH, O2 = O2, T = T)
names <- info(info(FeCu.cr))$abbrv
dFeCu <- diagram(mFeCu$A.species, plot.it = FALSE, names = names)

fill <- function(a) {
  ifelse(grepl("Ccp", a$species$name), "#FF8C0088",
    ifelse(grepl("Bn", a$species$name), "#DC143C88", NA)
)}
srt <- function(a) ifelse(a$species$name %in% c("Mag+Bn", "Mag+Ccp"), 80, 0)
cex.names <- function(a) ifelse(a$species$name %in% c("Mag+Bn", "Mag+Ccp", "Mag+Cct"), 0.8, 1)
dx <- function(a) sapply(a$species$name, switch, "Mag+Bn" = 0.15, "Mag+Cct" = 0.5, 0)
dy <- function(a) sapply(a$species$name, switch, "Py+Bn" = -1, "Po+Bn" = -0.8, "Po+Cu" = -0.8, "Mag+Ccp" = -1, 0)

a11 <- mix(dFe, dCu, dFeCu)
diagram(a11, fill = fill(a11), srt = srt(a11), min.area = 0.01, dx = dx(a11), dy = dy(a11), cex.names = cex.names(a11))
title("Fe:Cu = 1:1")

a21 <- mix(dFe, dCu, dFeCu, c(2, 1))
diagram(a21, fill = fill(a21), srt = srt(a21), min.area = 0.01, dx = dx(a21), dy = dy(a21), cex.names = cex.names(a21))
title("Fe:Cu = 2:1")

a12 <- mix(dFe, dCu, dFeCu, c(1, 2))
diagram(a12, fill = fill(a12), srt = srt(a12), min.area = 0.01, dx = dx(a12), dy = dy(a12), cex.names = cex.names(a12))
title("Fe:Cu = 1:2")

## ----stack3, eval = FALSE-----------------------------------------------------
#  T <- 125
#  layout(matrix(c(1, 2, 3, 3), nrow = 2), widths = c(1, 1.5))
#  
#  # Fe-S-O-H diagram
#  basis(c("copper", "hematite", "S2", "oxygen", "H2O", "H+", "Cl-"))
#  bFe <- species(c("hematite", "magnetite", "pyrite"))$name
#  aFe <- affinity(S2 = c(-34, -10, res1), O2 = c(-55, -40, res1), T = T)
#  # Order species by a function of composition to make colors
#  oFe <- order(aFe$species$S2 - aFe$species$O2)
#  fill <- terrain.colors(length(oFe), alpha = 0.2)[oFe]
#  abbrv <- info(aFe$species$ispecies)$abbrv
#  dFe <- diagram(aFe, names = abbrv, fill = fill)
#  title("Fe-S-O-H")
#  
#  # Cu-Fe-S-O-H diagram based on reactions with the
#  # stable Fe-bearing minerals (mosaic stack)
#  bCu <- species(c("copper", "chalcocite", "covellite", "chalcopyrite", "bornite"))$name
#  mCu <- mosaic(bFe, S2 = c(-34, -10, res1), O2 = c(-55, -40, res1),
#                T = T, stable = list(dFe$predominant))
#  oCu <- order(mCu$A.species$species$S2 - mCu$A.species$species$O2)
#  fill <- terrain.colors(length(oCu), alpha = 0.2)[oCu]
#  abbrv <- info(mCu$A.species$species$ispecies)$abbrv
#  dCu <- diagram(mCu$A.species, names = abbrv, fill = fill, dx = c(0, 0, 0, 0, 1.8))
#  title("Cu-Fe-S-O-H")
#  
#  # Mash the diagrams and adjust labels
#  aFeCu <- mash(dFe, dCu)
#  names <- aFeCu$species$name
#  dx <- dy <- srt <- rep(0, length(names))
#  cex <- rep(1, length(names))
#  cex[names %in% c("Hem+Ccp", "Hem+Cv")] <- 0.8
#  srt[names %in% c("Mag+Cu", "Hem+Cu")] <- 90
#  srt[names %in% c("Mag+Bn", "Hem+Bn")] <- 63
#  srt[names %in% c("Mag+Ccp", "Hem+Ccp")] <- 68
#  srt[names %in% c("Py+Bn", "Py+Cv")] <- 90
#  dx[names == "Hem+Ccp"] <- -0.4
#  dy[names == "Hem+Ccp"] <- -0.5
#  oFeCu <- order(aFeCu$species$S2 - aFeCu$species$O2)
#  fill <- terrain.colors(length(oFeCu), alpha = 0.2)[oFeCu]
#  diagram(aFeCu, cex.names = cex, srt = srt, dx = dx, dy = dy, fill = fill)
#  legend("topleft", legend = lTP(T, "Psat"), bg = "white")
#  title("Cu-Fe-S-O-H")

## ----stack3, echo = FALSE, results = "hide", message = FALSE, fig.width = 6, fig.height = 4, out.width = "100%", pngquant = pngquant----
T <- 125
layout(matrix(c(1, 2, 3, 3), nrow = 2), widths = c(1, 1.5))

# Fe-S-O-H diagram
basis(c("copper", "hematite", "S2", "oxygen", "H2O", "H+", "Cl-"))
bFe <- species(c("hematite", "magnetite", "pyrite"))$name
aFe <- affinity(S2 = c(-34, -10, res1), O2 = c(-55, -40, res1), T = T)
# Order species by a function of composition to make colors
oFe <- order(aFe$species$S2 - aFe$species$O2)
fill <- terrain.colors(length(oFe), alpha = 0.2)[oFe]
abbrv <- info(aFe$species$ispecies)$abbrv
dFe <- diagram(aFe, names = abbrv, fill = fill)
title("Fe-S-O-H")

# Cu-Fe-S-O-H diagram based on reactions with the
# stable Fe-bearing minerals (mosaic stack)
bCu <- species(c("copper", "chalcocite", "covellite", "chalcopyrite", "bornite"))$name
mCu <- mosaic(bFe, S2 = c(-34, -10, res1), O2 = c(-55, -40, res1),
              T = T, stable = list(dFe$predominant))
oCu <- order(mCu$A.species$species$S2 - mCu$A.species$species$O2)
fill <- terrain.colors(length(oCu), alpha = 0.2)[oCu]
abbrv <- info(mCu$A.species$species$ispecies)$abbrv
dCu <- diagram(mCu$A.species, names = abbrv, fill = fill, dx = c(0, 0, 0, 0, 1.8))
title("Cu-Fe-S-O-H")

# Mash the diagrams and adjust labels
aFeCu <- mash(dFe, dCu)
names <- aFeCu$species$name
dx <- dy <- srt <- rep(0, length(names))
cex <- rep(1, length(names))
cex[names %in% c("Hem+Ccp", "Hem+Cv")] <- 0.8
srt[names %in% c("Mag+Cu", "Hem+Cu")] <- 90
srt[names %in% c("Mag+Bn", "Hem+Bn")] <- 63
srt[names %in% c("Mag+Ccp", "Hem+Ccp")] <- 68
srt[names %in% c("Py+Bn", "Py+Cv")] <- 90
dx[names == "Hem+Ccp"] <- -0.4
dy[names == "Hem+Ccp"] <- -0.5
oFeCu <- order(aFeCu$species$S2 - aFeCu$species$O2)
fill <- terrain.colors(length(oFeCu), alpha = 0.2)[oFeCu]
diagram(aFeCu, cex.names = cex, srt = srt, dx = dx, dy = dy, fill = fill)
legend("topleft", legend = lTP(T, "Psat"), bg = "white")
title("Cu-Fe-S-O-H")

## ----solubility, eval = FALSE-------------------------------------------------
#  par(mfrow = c(1, 3))
#  basis("pH", 6)
#  # Estimate the molality of Cl for ca. 80,000 mg/kg solution (Table 2 of Sverjensky, 1987)
#  m_tot <- 80000 / mass("Cl") / 1000
#  calc <- NaCl(T = T, m_tot = m_tot)
#  # Use log molality here, not log activity, because
#  # activity coefficients are calculated by setting IS below
#  basis("Cl-", log10(calc$m_Cl))
#  # Initial setup: dissolve a single mineral to form aqueous Cu complexes
#  species("copper")
#  iaq <- retrieve("Cu", c("O", "H", "Cl", "S"), "aq")
#  
#  # Function to calculate solubility of Cu for stable assemblages of Fe and Cu minerals
#  # (i.e. equilibrium is imposed with all of these minerals, not only Cu(s))
#  mfun <- function() {
#    s <- solubility(iaq, bases = list(bFe, bCu), S2 = c(-34, -10, res1), O2 = c(-55, -40, res1),
#      T = T, IS = calc$IS, stable = list(dFe$predominant, dCu$predominant))
#    s <- convert(s, "ppm")
#    diagram(aFeCu, names = NA, col = "gray", fill = fill)
#    diagram(s, type = "loga.balance", levels = 10^(-3:3), add = TRUE)
#    diagram(s, type = "loga.balance", levels = 35, add = TRUE, lwd = 3, col = 6, contour.method = NA)
#  }
#  # DIAGRAM 1
#  mfun()
#  title("Cu (ppm)")
#  
#  # Calculate logK for CuCl2- dissociation at 125 °C
#  logK <- subcrt(c("CuCl2-", "Cu+", "Cl-"), c(-1, 1, 2), T = T)$out$logK
#  # Sverjensky (1987) used Helgeson (1969) value, which is ca. -5.2
#  dlogK <- logK - -5.2
#  # Calculate the difference in ΔG° corresponding to this logK difference
#  dG_J <- convert(dlogK, "G", T = convert(T, "K"))
#  # We should use calories here because the database values are in calories 20220604
#  stopifnot(info(info("CuCl2-"))$E_units == "cal")
#  dG_cal <- convert(dG_J, "cal")
#  # Apply this difference to the CuCl2- entry in OBIGT
#  newG <- info(info("CuCl2-"))$G + dG_cal
#  mod.OBIGT("CuCl2-", G = newG)
#  
#  # Do the same thing for CuCl3-2
#  logK <- subcrt(c("CuCl3-2", "Cu+", "Cl-"), c(-1, 1, 3), T = T)$out$logK
#  dlogK <- logK - -5.6
#  dG_J <- convert(dlogK, "G", T = convert(T, "K"))
#  stopifnot(info(info("CuCl3-2"))$E_units == "cal")
#  dG_cal <- convert(dG_J, "cal")
#  newG <- info(info("CuCl3-2"))$G + dG_cal
#  mod.OBIGT("CuCl3-2", G = newG)
#  
#  # DIAGRAM 2
#  mfun()
#  title("Cu (ppm)", line = 1.7)
#  CuCl2 <- expr.species("CuCl2-")
#  CuCl3 <- expr.species("CuCl3-2")
#  title(bquote("Helgeson (1969)"~.(CuCl2)~and~.(CuCl3)), line = 0.9)
#  
#  # Set up system to dissolve S2(gas)
#  basis(c("S2", "copper", "hematite", "oxygen", "H2O", "H+", "Cl-"))
#  basis("pH", 6)
#  species("S2")
#  # Calculate concentration of SO4-2
#  iaq <- info("SO4-2")
#  s <- solubility(iaq, S2 = c(-34, -10, res1), O2 = c(-55, -40, res1), T = T, IS = calc$IS, in.terms.of = "SO4-2")
#  s <- convert(s, "ppm")
#  # DIAGRAM 3
#  diagram(aFeCu, names = NA, col = "gray", fill = fill)
#  diagram(s, type = "loga.balance", levels = 10^(-3:3), add = TRUE)
#  diagram(s, type = "loga.balance", levels = 35, add = TRUE, lwd = 3, col = 6, contour.method = NA)
#  title(bquote(bold(.(expr.species("SO4-2"))~"(ppm)")))

## ----solubility, echo = FALSE, results = "hide", message = FALSE, fig.width = 7, fig.height = 3, out.width = "100%", fig.align = "center", pngquant = pngquant----
par(mfrow = c(1, 3))
basis("pH", 6)
# Estimate the molality of Cl for ca. 80,000 mg/kg solution (Table 2 of Sverjensky, 1987)
m_tot <- 80000 / mass("Cl") / 1000
calc <- NaCl(T = T, m_tot = m_tot)
# Use log molality here, not log activity, because
# activity coefficients are calculated by setting IS below
basis("Cl-", log10(calc$m_Cl))
# Initial setup: dissolve a single mineral to form aqueous Cu complexes
species("copper")
iaq <- retrieve("Cu", c("O", "H", "Cl", "S"), "aq")

# Function to calculate solubility of Cu for stable assemblages of Fe and Cu minerals
# (i.e. equilibrium is imposed with all of these minerals, not only Cu(s))
mfun <- function() {
  s <- solubility(iaq, bases = list(bFe, bCu), S2 = c(-34, -10, res1), O2 = c(-55, -40, res1),
    T = T, IS = calc$IS, stable = list(dFe$predominant, dCu$predominant))
  s <- convert(s, "ppm")
  diagram(aFeCu, names = NA, col = "gray", fill = fill)
  diagram(s, type = "loga.balance", levels = 10^(-3:3), add = TRUE)
  diagram(s, type = "loga.balance", levels = 35, add = TRUE, lwd = 3, col = 6, contour.method = NA)
}
# DIAGRAM 1
mfun()
title("Cu (ppm)")

# Calculate logK for CuCl2- dissociation at 125 °C
logK <- subcrt(c("CuCl2-", "Cu+", "Cl-"), c(-1, 1, 2), T = T)$out$logK
# Sverjensky (1987) used Helgeson (1969) value, which is ca. -5.2
dlogK <- logK - -5.2
# Calculate the difference in ΔG° corresponding to this logK difference
dG_J <- convert(dlogK, "G", T = convert(T, "K"))
# We should use calories here because the database values are in calories 20220604
stopifnot(info(info("CuCl2-"))$E_units == "cal")
dG_cal <- convert(dG_J, "cal")
# Apply this difference to the CuCl2- entry in OBIGT
newG <- info(info("CuCl2-"))$G + dG_cal
mod.OBIGT("CuCl2-", G = newG)

# Do the same thing for CuCl3-2
logK <- subcrt(c("CuCl3-2", "Cu+", "Cl-"), c(-1, 1, 3), T = T)$out$logK
dlogK <- logK - -5.6
dG_J <- convert(dlogK, "G", T = convert(T, "K"))
stopifnot(info(info("CuCl3-2"))$E_units == "cal")
dG_cal <- convert(dG_J, "cal")
newG <- info(info("CuCl3-2"))$G + dG_cal
mod.OBIGT("CuCl3-2", G = newG)

# DIAGRAM 2
mfun()
title("Cu (ppm)", line = 1.7)
CuCl2 <- expr.species("CuCl2-")
CuCl3 <- expr.species("CuCl3-2")
title(bquote("Helgeson (1969)"~.(CuCl2)~and~.(CuCl3)), line = 0.9)

# Set up system to dissolve S2(gas)
basis(c("S2", "copper", "hematite", "oxygen", "H2O", "H+", "Cl-"))
basis("pH", 6)
species("S2")
# Calculate concentration of SO4-2
iaq <- info("SO4-2")
s <- solubility(iaq, S2 = c(-34, -10, res1), O2 = c(-55, -40, res1), T = T, IS = calc$IS, in.terms.of = "SO4-2")
s <- convert(s, "ppm")
# DIAGRAM 3
diagram(aFeCu, names = NA, col = "gray", fill = fill)
diagram(s, type = "loga.balance", levels = 10^(-3:3), add = TRUE)
diagram(s, type = "loga.balance", levels = 35, add = TRUE, lwd = 3, col = 6, contour.method = NA)
title(bquote(bold(.(expr.species("SO4-2"))~"(ppm)")))

## ----NaCl---------------------------------------------------------------------
# Ionic strength
calc$IS
# Logarithm of activity of Cl-
log10(calc$m_Cl * calc$gam_Cl)

## ----logK, message = FALSE----------------------------------------------------
# logK values interpolated from Table 5 of Helgeson (1969)
subcrt(c("CuCl2-", "Cu+", "Cl-"), c(-1, 1, 2), T = T)$out$logK
subcrt(c("CuCl3-2", "Cu+", "Cl-"), c(-1, 1, 3), T = T)$out$logK
# Default OBIGT database
reset()
subcrt(c("CuCl2-", "Cu+", "Cl-"), c(-1, 1, 2), T = T)$out$logK
subcrt(c("CuCl3-2", "Cu+", "Cl-"), c(-1, 1, 3), T = T)$out$logK

## ----iCu.aq-------------------------------------------------------------------
names(iCu.aq)

## ----rebalance, eval = FALSE--------------------------------------------------
#  mat <- matrix(c(1, 1, 2, 2, 3, 3, 4, 4, 4, 5, 5, 5), byrow = TRUE, nrow = 2)
#  layout(mat)
#  par(font.main = 1)
#  
#  basis(c("Fe+2", "Cu+", "hydrogen sulfide", "oxygen", "H2O", "H+"))
#  xlab <- ratlab("Fe+2", "Cu+")
#  
#  ### PRIMARY balancing
#  
#  # Only Fe-bearing minerals
#  species(c("pyrite", "pyrrhotite", "magnetite", "hematite"))
#  aFe <- affinity("Fe+2" = c(0, 12), O2 = c(-40, -16), T = 400, P = 2000)
#  names <- info(aFe$species$ispecies)$abbrv
#  dFe <- diagram(aFe, xlab = xlab, names = names)
#  title(bquote("Only Fe; 1° balance:" ~ .(expr.species(dFe$balance))))
#  label.plot("A")
#  
#  # Only Cu-bearing minerals
#  species(c("covellite", "chalcocite", "tenorite", "cuprite"))
#  aCu <- affinity(aFe)  # argument recall
#  names <- info(aCu$species$ispecies)$abbrv
#  dCu <- diagram(aCu, xlab = xlab, names = names)
#  title(bquote("Only Cu; 1° balance:" ~ .(expr.species(dCu$balance))))
#  label.plot("B")
#  
#  # Only Fe- AND Cu-bearing minerals
#  species(c("chalcopyrite", "bornite"))
#  aFeCu <- affinity(aFe)
#  names <- info(aFeCu$species$ispecies)$abbrv
#  dFeCu <- diagram(aFeCu, xlab = xlab, balance = "H+", names = names)
#  title(bquote("Only Fe+Cu; 1° balance:" ~ .(expr.species(dFeCu$balance))))
#  label.plot("C")
#  
#  ### SECONDARY balancing
#  
#  # Fe- or Cu-bearing minerals
#  ad1 <- rebalance(dFe, dCu, balance = "H+")
#  names <- info(ad1$species$ispecies)$abbrv
#  d1 <- diagram(ad1, xlab = xlab, balance = 1, names = names)
#  title(bquote("Only Fe or Cu; 2° balance:" ~ .(expr.species("H+"))))
#  label.plot("D")
#  
#  # All minerals
#  d1$values <- c(dFe$values, dCu$values)
#  ad2 <- rebalance(d1, dFeCu, balance = "H+")
#  names <- info(ad2$species$ispecies)$abbrv
#  diagram(ad2, xlab = xlab, balance = 1, names = names)
#  title(bquote("Fe and/or Cu; 2° balance:" ~ .(expr.species("H+"))))
#  label.plot("E")
#  
#  db <- describe.basis(3)
#  leg <- lex(lTP(400, 2000), db)
#  legend("bottomleft", legend = leg, bty = "n")

## ----rebalance, eval = FALSE, echo = 5:6--------------------------------------
#  mat <- matrix(c(1, 1, 2, 2, 3, 3, 4, 4, 4, 5, 5, 5), byrow = TRUE, nrow = 2)
#  layout(mat)
#  par(font.main = 1)
#  
#  basis(c("Fe+2", "Cu+", "hydrogen sulfide", "oxygen", "H2O", "H+"))
#  xlab <- ratlab("Fe+2", "Cu+")
#  
#  ### PRIMARY balancing
#  
#  # Only Fe-bearing minerals
#  species(c("pyrite", "pyrrhotite", "magnetite", "hematite"))
#  aFe <- affinity("Fe+2" = c(0, 12), O2 = c(-40, -16), T = 400, P = 2000)
#  names <- info(aFe$species$ispecies)$abbrv
#  dFe <- diagram(aFe, xlab = xlab, names = names)
#  title(bquote("Only Fe; 1° balance:" ~ .(expr.species(dFe$balance))))
#  label.plot("A")
#  
#  # Only Cu-bearing minerals
#  species(c("covellite", "chalcocite", "tenorite", "cuprite"))
#  aCu <- affinity(aFe)  # argument recall
#  names <- info(aCu$species$ispecies)$abbrv
#  dCu <- diagram(aCu, xlab = xlab, names = names)
#  title(bquote("Only Cu; 1° balance:" ~ .(expr.species(dCu$balance))))
#  label.plot("B")
#  
#  # Only Fe- AND Cu-bearing minerals
#  species(c("chalcopyrite", "bornite"))
#  aFeCu <- affinity(aFe)
#  names <- info(aFeCu$species$ispecies)$abbrv
#  dFeCu <- diagram(aFeCu, xlab = xlab, balance = "H+", names = names)
#  title(bquote("Only Fe+Cu; 1° balance:" ~ .(expr.species(dFeCu$balance))))
#  label.plot("C")
#  
#  ### SECONDARY balancing
#  
#  # Fe- or Cu-bearing minerals
#  ad1 <- rebalance(dFe, dCu, balance = "H+")
#  names <- info(ad1$species$ispecies)$abbrv
#  d1 <- diagram(ad1, xlab = xlab, balance = 1, names = names)
#  title(bquote("Only Fe or Cu; 2° balance:" ~ .(expr.species("H+"))))
#  label.plot("D")
#  
#  # All minerals
#  d1$values <- c(dFe$values, dCu$values)
#  ad2 <- rebalance(d1, dFeCu, balance = "H+")
#  names <- info(ad2$species$ispecies)$abbrv
#  diagram(ad2, xlab = xlab, balance = 1, names = names)
#  title(bquote("Fe and/or Cu; 2° balance:" ~ .(expr.species("H+"))))
#  label.plot("E")
#  
#  db <- describe.basis(3)
#  leg <- lex(lTP(400, 2000), db)
#  legend("bottomleft", legend = leg, bty = "n")

## ----rebalance, eval = FALSE, echo = 10:33------------------------------------
#  mat <- matrix(c(1, 1, 2, 2, 3, 3, 4, 4, 4, 5, 5, 5), byrow = TRUE, nrow = 2)
#  layout(mat)
#  par(font.main = 1)
#  
#  basis(c("Fe+2", "Cu+", "hydrogen sulfide", "oxygen", "H2O", "H+"))
#  xlab <- ratlab("Fe+2", "Cu+")
#  
#  ### PRIMARY balancing
#  
#  # Only Fe-bearing minerals
#  species(c("pyrite", "pyrrhotite", "magnetite", "hematite"))
#  aFe <- affinity("Fe+2" = c(0, 12), O2 = c(-40, -16), T = 400, P = 2000)
#  names <- info(aFe$species$ispecies)$abbrv
#  dFe <- diagram(aFe, xlab = xlab, names = names)
#  title(bquote("Only Fe; 1° balance:" ~ .(expr.species(dFe$balance))))
#  label.plot("A")
#  
#  # Only Cu-bearing minerals
#  species(c("covellite", "chalcocite", "tenorite", "cuprite"))
#  aCu <- affinity(aFe)  # argument recall
#  names <- info(aCu$species$ispecies)$abbrv
#  dCu <- diagram(aCu, xlab = xlab, names = names)
#  title(bquote("Only Cu; 1° balance:" ~ .(expr.species(dCu$balance))))
#  label.plot("B")
#  
#  # Only Fe- AND Cu-bearing minerals
#  species(c("chalcopyrite", "bornite"))
#  aFeCu <- affinity(aFe)
#  names <- info(aFeCu$species$ispecies)$abbrv
#  dFeCu <- diagram(aFeCu, xlab = xlab, balance = "H+", names = names)
#  title(bquote("Only Fe+Cu; 1° balance:" ~ .(expr.species(dFeCu$balance))))
#  label.plot("C")
#  
#  ### SECONDARY balancing
#  
#  # Fe- or Cu-bearing minerals
#  ad1 <- rebalance(dFe, dCu, balance = "H+")
#  names <- info(ad1$species$ispecies)$abbrv
#  d1 <- diagram(ad1, xlab = xlab, balance = 1, names = names)
#  title(bquote("Only Fe or Cu; 2° balance:" ~ .(expr.species("H+"))))
#  label.plot("D")
#  
#  # All minerals
#  d1$values <- c(dFe$values, dCu$values)
#  ad2 <- rebalance(d1, dFeCu, balance = "H+")
#  names <- info(ad2$species$ispecies)$abbrv
#  diagram(ad2, xlab = xlab, balance = 1, names = names)
#  title(bquote("Fe and/or Cu; 2° balance:" ~ .(expr.species("H+"))))
#  label.plot("E")
#  
#  db <- describe.basis(3)
#  leg <- lex(lTP(400, 2000), db)
#  legend("bottomleft", legend = leg, bty = "n")

## ----rebalance, eval = FALSE, echo = 36:49------------------------------------
#  mat <- matrix(c(1, 1, 2, 2, 3, 3, 4, 4, 4, 5, 5, 5), byrow = TRUE, nrow = 2)
#  layout(mat)
#  par(font.main = 1)
#  
#  basis(c("Fe+2", "Cu+", "hydrogen sulfide", "oxygen", "H2O", "H+"))
#  xlab <- ratlab("Fe+2", "Cu+")
#  
#  ### PRIMARY balancing
#  
#  # Only Fe-bearing minerals
#  species(c("pyrite", "pyrrhotite", "magnetite", "hematite"))
#  aFe <- affinity("Fe+2" = c(0, 12), O2 = c(-40, -16), T = 400, P = 2000)
#  names <- info(aFe$species$ispecies)$abbrv
#  dFe <- diagram(aFe, xlab = xlab, names = names)
#  title(bquote("Only Fe; 1° balance:" ~ .(expr.species(dFe$balance))))
#  label.plot("A")
#  
#  # Only Cu-bearing minerals
#  species(c("covellite", "chalcocite", "tenorite", "cuprite"))
#  aCu <- affinity(aFe)  # argument recall
#  names <- info(aCu$species$ispecies)$abbrv
#  dCu <- diagram(aCu, xlab = xlab, names = names)
#  title(bquote("Only Cu; 1° balance:" ~ .(expr.species(dCu$balance))))
#  label.plot("B")
#  
#  # Only Fe- AND Cu-bearing minerals
#  species(c("chalcopyrite", "bornite"))
#  aFeCu <- affinity(aFe)
#  names <- info(aFeCu$species$ispecies)$abbrv
#  dFeCu <- diagram(aFeCu, xlab = xlab, balance = "H+", names = names)
#  title(bquote("Only Fe+Cu; 1° balance:" ~ .(expr.species(dFeCu$balance))))
#  label.plot("C")
#  
#  ### SECONDARY balancing
#  
#  # Fe- or Cu-bearing minerals
#  ad1 <- rebalance(dFe, dCu, balance = "H+")
#  names <- info(ad1$species$ispecies)$abbrv
#  d1 <- diagram(ad1, xlab = xlab, balance = 1, names = names)
#  title(bquote("Only Fe or Cu; 2° balance:" ~ .(expr.species("H+"))))
#  label.plot("D")
#  
#  # All minerals
#  d1$values <- c(dFe$values, dCu$values)
#  ad2 <- rebalance(d1, dFeCu, balance = "H+")
#  names <- info(ad2$species$ispecies)$abbrv
#  diagram(ad2, xlab = xlab, balance = 1, names = names)
#  title(bquote("Fe and/or Cu; 2° balance:" ~ .(expr.species("H+"))))
#  label.plot("E")
#  
#  db <- describe.basis(3)
#  leg <- lex(lTP(400, 2000), db)
#  legend("bottomleft", legend = leg, bty = "n")

## ----rebalance, echo = FALSE, results = "hide", message = FALSE, fig.width = 6.5, fig.height = 5, out.width = "100%", fig.align = "center", pngquant = pngquant----
mat <- matrix(c(1, 1, 2, 2, 3, 3, 4, 4, 4, 5, 5, 5), byrow = TRUE, nrow = 2)
layout(mat)
par(font.main = 1)

basis(c("Fe+2", "Cu+", "hydrogen sulfide", "oxygen", "H2O", "H+"))
xlab <- ratlab("Fe+2", "Cu+")

### PRIMARY balancing

# Only Fe-bearing minerals
species(c("pyrite", "pyrrhotite", "magnetite", "hematite"))
aFe <- affinity("Fe+2" = c(0, 12), O2 = c(-40, -16), T = 400, P = 2000)
names <- info(aFe$species$ispecies)$abbrv
dFe <- diagram(aFe, xlab = xlab, names = names)
title(bquote("Only Fe; 1° balance:" ~ .(expr.species(dFe$balance))))
label.plot("A")

# Only Cu-bearing minerals
species(c("covellite", "chalcocite", "tenorite", "cuprite"))
aCu <- affinity(aFe)  # argument recall
names <- info(aCu$species$ispecies)$abbrv
dCu <- diagram(aCu, xlab = xlab, names = names)
title(bquote("Only Cu; 1° balance:" ~ .(expr.species(dCu$balance))))
label.plot("B")

# Only Fe- AND Cu-bearing minerals
species(c("chalcopyrite", "bornite"))
aFeCu <- affinity(aFe)
names <- info(aFeCu$species$ispecies)$abbrv
dFeCu <- diagram(aFeCu, xlab = xlab, balance = "H+", names = names)
title(bquote("Only Fe+Cu; 1° balance:" ~ .(expr.species(dFeCu$balance))))
label.plot("C")

### SECONDARY balancing

# Fe- or Cu-bearing minerals
ad1 <- rebalance(dFe, dCu, balance = "H+")
names <- info(ad1$species$ispecies)$abbrv
d1 <- diagram(ad1, xlab = xlab, balance = 1, names = names)
title(bquote("Only Fe or Cu; 2° balance:" ~ .(expr.species("H+"))))
label.plot("D")

# All minerals
d1$values <- c(dFe$values, dCu$values)
ad2 <- rebalance(d1, dFeCu, balance = "H+")
names <- info(ad2$species$ispecies)$abbrv
diagram(ad2, xlab = xlab, balance = 1, names = names)
title(bquote("Fe and/or Cu; 2° balance:" ~ .(expr.species("H+"))))
label.plot("E")

db <- describe.basis(3)
leg <- lex(lTP(400, 2000), db)
legend("bottomleft", legend = leg, bty = "n")

## ----non-metal, results = "hide", message = FALSE, fig.width = 6, fig.height = 5, out.width = "80%", fig.align = "center", pngquant = pngquant----
basis(c("Fe+2", "Cu+", "hydrogen sulfide", "oxygen", "H2O", "H+"))
basis("H2S", 2)
species(c("pyrite", "magnetite", "hematite", "covellite", "tenorite",
          "chalcopyrite", "bornite"))
a <- affinity("Cu+" = c(-8, 2, 500), "Fe+2" = c(-4, 12, 500), T = 400, P = 2000)
names <- info(a$species$ispecies)$abbrv
d <- diagram(a, xlab = ratlab("Cu+"), ylab = ratlab("Fe+2"), balance = "O2", names = names)
title(bquote("Cu-Fe-S-O-H; 1° balance:" ~ .(expr.species(d$balance))))
# Add saturation lines
species(c("pyrrhotite", "ferrous-oxide", "chalcocite", "cuprite"))
asat <- affinity(a)  # argument recall
names <- asat$species$name
names[2] <- "ferrous oxide"
diagram(asat, type = "saturation", add = TRUE, lty = 2, col = 4, names = names)
legend("topleft", legend = lTP(400, 2000), bty = "n")

## ----mosaic-combo, eval = FALSE-----------------------------------------------
#  ## METHOD 1: Keff equation (Robinson et al., 2021)
#  
#  # A function to calculate Keff for any combination of T and pH
#  Keff <- function(T = 25, pH = 7) {
#    # Make T and pH the same length
#    len <- max(length(T), length(pH))
#    T <- rep(T, length.out = len)
#    pH <- rep(pH, length.out = len)
#    # Calculate activity of H+
#    aH <- 10^(-pH)
#    # Calculate logKs
#    logK1 <- subcrt(c("acetic acid", "NH3", "acetamide", "H2O"), c(-1, -1, 1, 1), T = T)$out$logK
#    logK2 <- subcrt(c("acetic acid", "acetate", "H+"), c(-1, 1, 1), T = T)$out$logK
#    logK3 <- subcrt(c("NH3", "H+", "NH4+"), c(-1, -1, 1), T = T)$out$logK
#    # Calculate Ks
#    K1 <- 10^logK1
#    K2 <- 10^logK2
#    K3 <- 10^logK3
#    # Calculate Keff (Eq. 7)
#    Keff <- K1 * (1 + K2 / aH) ^ -1 * (1 + K3 * aH) ^ -1
#    Keff
#  }
#  
#  # Calculate logKeff as a function of pH at 100 °C
#  res <- 128
#  pH <- seq(0, 14, length.out = res)
#  T <- 100
#  logKeff <- log10(Keff(pH = pH, T = T))
#  
#  # Calculate activity of acetamide for
#  # acetic acid + acetate = 0.01 m
#  # ammonia + ammonium = 0.001 m
#  logAc <- log10(0.01)
#  logAm <- log10(0.001)
#  logAcAm <- logKeff + logAc + logAm
#  
#  ## METHOD 2: Mosaic combo
#  
#  # Define total activities
#  a_N <- 0.001
#  # This is 2 * 0.01 because acetic acid has 2 carbons
#  a_C <- 2 * 0.01
#  loga_N <- log10(a_N)
#  loga_C <- log10(a_C)
#  # Setup basis species
#  basis(c("CO2", "NH3", "O2", "H2O", "H+"))
#  basis("NH3", loga_N)
#  # Load all C-bearing species (including acetamide)
#  species(c("acetamide", "acetic acid", "acetate"))
#  # Calculate distribution of C-bearing species accounting for ammonia/ammonium speciation
#  m <- mosaic(c("NH3", "NH4+"), pH = c(0, 14, res), T = T)
#  e <- equilibrate(m, loga.balance = loga_C)
#  
#  # Plot and label diagram
#  # Start with empty diagram
#  diagram(e, ylim = c(-8, -0), lty = 0, names = FALSE)
#  # Add pH = 6 line
#  abline(v = 6, col = "gray60", lty = 5)
#  # Add line for acetamide activity calculated with Keff
#  lines(pH, logAcAm, col = 4, lwd = 6, lty = 2)
#  # Add lines from CHNOSZ calculations
#  diagram(e, add = TRUE, lty = c(2, 3, 1, 2, 3), lwd = c(1, 1, 2, 1, 1), dx = c(-0.2, 0.2, -2.5, 0, 0), dy = c(0.1, 0.1, -1, 0.1, 0.1), srt = c(0, 0, 52, 0, 0))
#  tN <- paste("Total N in basis species =", format(a_N, scientific = FALSE), "m")
#  tC <- paste("Total C in formed species =", format(a_C, scientific = FALSE), "m")
#  title(main = paste(tN, tC, sep = "\n"), font.main = 1)
#  legend("topright", legend = lTP(T, "Psat"), bty = "n")
#  # Check that we got equal values
#  stopifnot(all.equal(as.numeric(e$loga.equil[[3]]), logAcAm, tol = 1e-3, scale = 1))

## ----mosaic-combo, echo = FALSE, results = "hide", message = FALSE, fig.width = 6, fig.height = 5, out.width = "80%", fig.align = "center", pngquant = pngquant----
## METHOD 1: Keff equation (Robinson et al., 2021)

# A function to calculate Keff for any combination of T and pH
Keff <- function(T = 25, pH = 7) {
  # Make T and pH the same length
  len <- max(length(T), length(pH))
  T <- rep(T, length.out = len)
  pH <- rep(pH, length.out = len)
  # Calculate activity of H+
  aH <- 10^(-pH)
  # Calculate logKs
  logK1 <- subcrt(c("acetic acid", "NH3", "acetamide", "H2O"), c(-1, -1, 1, 1), T = T)$out$logK
  logK2 <- subcrt(c("acetic acid", "acetate", "H+"), c(-1, 1, 1), T = T)$out$logK
  logK3 <- subcrt(c("NH3", "H+", "NH4+"), c(-1, -1, 1), T = T)$out$logK
  # Calculate Ks
  K1 <- 10^logK1
  K2 <- 10^logK2
  K3 <- 10^logK3
  # Calculate Keff (Eq. 7)
  Keff <- K1 * (1 + K2 / aH) ^ -1 * (1 + K3 * aH) ^ -1
  Keff
}

# Calculate logKeff as a function of pH at 100 °C
res <- 128
pH <- seq(0, 14, length.out = res)
T <- 100
logKeff <- log10(Keff(pH = pH, T = T))

# Calculate activity of acetamide for
# acetic acid + acetate = 0.01 m
# ammonia + ammonium = 0.001 m
logAc <- log10(0.01)
logAm <- log10(0.001)
logAcAm <- logKeff + logAc + logAm

## METHOD 2: Mosaic combo

# Define total activities
a_N <- 0.001
# This is 2 * 0.01 because acetic acid has 2 carbons
a_C <- 2 * 0.01
loga_N <- log10(a_N)
loga_C <- log10(a_C)
# Setup basis species
basis(c("CO2", "NH3", "O2", "H2O", "H+"))
basis("NH3", loga_N)
# Load all C-bearing species (including acetamide)
species(c("acetamide", "acetic acid", "acetate"))
# Calculate distribution of C-bearing species accounting for ammonia/ammonium speciation
m <- mosaic(c("NH3", "NH4+"), pH = c(0, 14, res), T = T)
e <- equilibrate(m, loga.balance = loga_C)

# Plot and label diagram
# Start with empty diagram
diagram(e, ylim = c(-8, -0), lty = 0, names = FALSE)
# Add pH = 6 line
abline(v = 6, col = "gray60", lty = 5)
# Add line for acetamide activity calculated with Keff
lines(pH, logAcAm, col = 4, lwd = 6, lty = 2)
# Add lines from CHNOSZ calculations
diagram(e, add = TRUE, lty = c(2, 3, 1, 2, 3), lwd = c(1, 1, 2, 1, 1), dx = c(-0.2, 0.2, -2.5, 0, 0), dy = c(0.1, 0.1, -1, 0.1, 0.1), srt = c(0, 0, 52, 0, 0))
tN <- paste("Total N in basis species =", format(a_N, scientific = FALSE), "m")
tC <- paste("Total C in formed species =", format(a_C, scientific = FALSE), "m")
title(main = paste(tN, tC, sep = "\n"), font.main = 1)
legend("topright", legend = lTP(T, "Psat"), bty = "n")
# Check that we got equal values
stopifnot(all.equal(as.numeric(e$loga.equil[[3]]), logAcAm, tol = 1e-3, scale = 1))

