## ----setup, include=FALSE-----------------------------------------------------
options(width = 80)
## use pngquant to optimize PNG images
library(knitr)
knit_hooks$set(pngquant = hook_pngquant)
pngquant <- "--speed=1 --quality=0-25"
if (!nzchar(Sys.which("pngquant"))) pngquant <- NULL
## logK with a thin space 20200627
logK <- "log&thinsp;<i>K</i>"

## ----CHNOSZ_reset, include=FALSE----------------------------------------------
library(CHNOSZ)
reset()

## ----mash, echo = 1:8, eval = FALSE-------------------------------------------
#  par(mfrow = c(1, 2))
#  basis("CHNOS+")
#  species(c("CH4", "CO2", "HCO3-", "CO3-2"))
#  aC <- affinity(pH = c(0, 14), O2 = c(-75, -60))
#  dC <- diagram(aC)
#  species(c("H2S", "HS-", "HSO4-", "SO4-2"))
#  aS <- affinity(aC)  # argument recall
#  dS <- diagram(aS, add = TRUE, col = 4, col.names = 4)
#  aCS <- mash(dC, dS)
#  diagram(aCS)
#  legend("topright", legend = lTP(25, 1), bty = "n")

## ----mash, echo = 9:11,  results = "hide", message = FALSE, fig.width = 10, fig.height = 5, out.width = "100%"----
par(mfrow = c(1, 2))
basis("CHNOS+")
species(c("CH4", "CO2", "HCO3-", "CO3-2"))
aC <- affinity(pH = c(0, 14), O2 = c(-75, -60))
dC <- diagram(aC)
species(c("H2S", "HS-", "HSO4-", "SO4-2"))
aS <- affinity(aC)  # argument recall
dS <- diagram(aS, add = TRUE, col = 4, col.names = 4)
aCS <- mash(dC, dS)
diagram(aCS)
legend("topright", legend = lTP(25, 1), bty = "n")

## ----materials, message = FALSE, results = "hide"-----------------------------
## Formation energies (eV / atom) for solids from Materials API, e.g.
# from pymatgen import MPRester
# m = MPRester("USER_API_KEY")
# m.query(criteria={"task_id": "mp-1279742"}, properties=["formation_energy_per_atom"])
# mp-13, mp-1279742, mp-19306, mp-19770
Fe.cr <- c(Fe = 0, FeO = -1.72803, Fe3O4 = -1.85868, Fe2O3 = -1.90736)  # 20201109
# mp-146, mp-18937, mp-1275946, mp-19094, mp-756395, mp-25279
V.cr <- c(V = 0, V2O3 = -2.52849, V3O5 = -2.52574, VO2 = -2.48496, V3O7 = -2.32836, V2O5 = -2.29524)  # 20201109

# Convert formation energies from eV / atom to eV / molecule
natom.Fe <- sapply(makeup(names(Fe.cr)), sum)
Fe.cr <- Fe.cr * natom.Fe
natom.V <- sapply(makeup(names(V.cr)), sum)
V.cr <- V.cr * natom.V

# Convert formation energies from eV / molecule to J / mol
eVtoJ <- function(eV) eV * 1.602176634e-19 * 6.02214076e23
Fe.cr <- eVtoJ(Fe.cr)
V.cr <- eVtoJ(V.cr)

# Gibbs energies of formation (J / mol) for aqueous species from Wagman et al.
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
# Be sure to set E_units to Joules; the default in OBIGT is calories!
# This function modifies OBIGT and returns the species indices of the affected species
modfun <- function(x, state) sapply(seq_along(x), function(i) {
  mod.OBIGT(names(x)[i], formula = names(x)[i], state = state, E_units = "J", G = x[i])
})
iFe.cr <- modfun(Fe.cr, "cr")
iFe.aq <- modfun(Fe.aq, "aq")
iV.cr <- modfun(V.cr, "cr")
iV.aq <- modfun(V.aq, "aq")

# Formation energies (eV / atom) for bimetallic solids from Materials API
# mp-1335, mp-1079399, mp-866134, mp-558525, mp-504509 (triclinic FeVO4)
FeV.cr <- c(FeV = -0.12928, FeV3 = -0.17128, Fe3V = -0.12854, Fe2V4O13 = -2.23833, FeVO4 = -1.75611)  # 20201109
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
#  aFe <- affinity(pH = c(4, 10), Eh = c(-1.5, 0))
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
#  diagram(a11, min.area = 0.01)
#  title("Fe:V = 1:1")
#  label.figure(lTP(25, 1), xfrac = 0.12)
#  # 1:3 mixture
#  a13 <- mix(dFe, dV, dFeV, c(1, 3))
#  diagram(a13, min.area = 0.01)
#  title("Fe:V = 1:3")
#  # 1:5 mixture
#  a15 <- mix(dFe, dV, dFeV, c(1, 5))
#  diagram(a15, min.area = 0.01)
#  title("Fe:V = 1:5")

## ----mixing1, echo = 16:37, message = FALSE, results = "hide", fig.width = 9, fig.height = 3, out.width = "100%", out.extra='class="full-width"', pngquant = pngquant----
par(mfrow = c(1, 3))
loga.Fe <- -5
loga.V <- -5
# Fe-O-H diagram
basis(c("VO+2", "Fe+2", "H2O", "e-", "H+"))
species(c(iFe.aq, iFe.cr))
species(1:length(iFe.aq), loga.Fe)
aFe <- affinity(pH = c(4, 10), Eh = c(-1.5, 0))
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
diagram(a11, min.area = 0.01)
title("Fe:V = 1:1")
label.figure(lTP(25, 1), xfrac = 0.12)
# 1:3 mixture
a13 <- mix(dFe, dV, dFeV, c(1, 3))
diagram(a13, min.area = 0.01)
title("Fe:V = 1:3")
# 1:5 mixture
a15 <- mix(dFe, dV, dFeV, c(1, 5))
diagram(a15, min.area = 0.01)
title("Fe:V = 1:5")

## ----FeVO4, eval = FALSE, echo = 1:20-----------------------------------------
#  par(mfrow = c(1, 2))
#  # Fe-bearing species
#  basis(c("VO+2", "Fe+2", "H2O", "e-", "H+"))
#  species(c(iFe.aq, iFe.cr))$name
#  species(1:length(iFe.aq), loga.Fe)
#  aFe <- affinity(pH = c(0, 14, 400), Eh = c(-1.5, 2, 400))
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
#  d11 <- diagram(a11, min.area = 0.01)
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
#  image(aFeVO4_vs_stable, col = rev(topo.colors(100, 0.7)), add = TRUE)
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

## ----FeVO4, echo = 22:34, message = FALSE, results = "hide", fig.width = 10, fig.height = 5, out.width = "100%", pngquant = pngquant----
par(mfrow = c(1, 2))
# Fe-bearing species
basis(c("VO+2", "Fe+2", "H2O", "e-", "H+"))
species(c(iFe.aq, iFe.cr))$name
species(1:length(iFe.aq), loga.Fe)
aFe <- affinity(pH = c(0, 14, 400), Eh = c(-1.5, 2, 400))
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
d11 <- diagram(a11, min.area = 0.01)
water.lines(d11, col = "orangered")

# Calculate affinity of FeVO4
species("FeVO4")
aFeVO4 <- affinity(aFe)  # argument recall
# Calculate difference from stable species
aFeVO4_vs_stable <- aFeVO4$values[[1]] - d11$predominant.values
# Overlay lines from diagram on color map
diagram(a11, fill = NA, names = FALSE, limit.water = FALSE)
opar <- par(usr = c(0, 1, 0, 1))
image(aFeVO4_vs_stable, col = rev(topo.colors(100, 0.7)), add = TRUE)
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

## ----Gpbx_min, echo = 2:8, message = FALSE, fig.keep = "none"-----------------
plot(1:10) # so we can run "points" in this chunk
imax <- arrayInd(which.max(aFeVO4_vs_stable), dim(aFeVO4_vs_stable))
pH <- d11$vals$pH[imax[1]]
Eh <- d11$vals$Eh[imax[2]]
points(pH, Eh, pch = 10, cex = 2, lwd = 2, col = "gold")
stable <- d11$names[d11$predominant[imax]]
text(pH, Eh, stable, adj = c(0.3, 2), cex = 1.2, col = "gold")
range(aFeVO4_vs_stable[d11$predominant == d11$predominant[imax]])

## ----hull, echo = 1:7, message = FALSE----------------------------------------
b <- basis(c("Fe2O3", "Fe2V4O13", "O2"))
cal_mol <- subcrt("FeVO4", 1, T = 25)$out$G
convert(cal_mol, "logK")
J_mol <- convert(cal_mol, "J")
eV_mol <- J_mol / 1.602176634e-19
eV_atom <- eV_mol / 6.02214076e23 / 6
round(eV_atom, 3)
stopifnot(round(eV_atom, 3) == 0.413)

## ----reset, message = FALSE---------------------------------------------------
reset()

## ----stack1_1, results = "hide", message = FALSE------------------------------
logaH2S <- -2
T <- 200
pH <- c(0, 14, 500)
O2 <- c(-45, -30, 500)
basis(c("Cu+", "pyrite", "H2S", "oxygen", "H2O", "H+"))
basis("H2S", logaH2S)
S.aq <- c("H2S", "HS-", "HSO4-", "SO4-2")
Fe.cr <- c("pyrite", "pyrrhotite", "magnetite", "hematite")

## ----stack1_2, eval = FALSE, echo = 1:4---------------------------------------
#  species(Fe.cr)
#  mFe <- mosaic(S.aq, pH = pH, O2 = O2, T = T)
#  diagram(mFe$A.bases, lty = 2, col = 4, col.names = 4, italic = TRUE)
#  dFe <- diagram(mFe$A.species, add = TRUE, lwd = 2)
#  FeCu.cr <- c("chalcopyrite", "bornite")
#  Cu.cr <- c("copper", "cuprite", "tenorite", "chalcocite", "covellite")
#  species(c(FeCu.cr, Cu.cr))
#  mFeCu <- mosaic(list(S.aq, Fe.cr), pH = pH, O2 = O2,
#                T = T, stable = list(NULL, dFe$predominant))
#  diagram(mFeCu$A.species, add = TRUE, col = 2, col.names = 2, bold = TRUE)
#  col <- c("#FF8C00", rep(NA, 6))
#  diagram(mFeCu$A.species, add = TRUE, col = col, lwd = 2, col.names = col, bold = TRUE)
#  TP <- describe.property(c("T", "P"), c(T, "Psat"))
#  legend("topright", TP, bty = "n")
#  title(paste("Cu-Fe-S-O-H; Total S =", 10^logaH2S, "m"))

## ----stack1_2, echo=5:15, results = "hide", message = FALSE, fig.width = 6, fig.height = 5, out.width = "75%", fig.align = "center", pngquant = pngquant----
species(Fe.cr)
mFe <- mosaic(S.aq, pH = pH, O2 = O2, T = T)
diagram(mFe$A.bases, lty = 2, col = 4, col.names = 4, italic = TRUE)
dFe <- diagram(mFe$A.species, add = TRUE, lwd = 2)
FeCu.cr <- c("chalcopyrite", "bornite")
Cu.cr <- c("copper", "cuprite", "tenorite", "chalcocite", "covellite")
species(c(FeCu.cr, Cu.cr))
mFeCu <- mosaic(list(S.aq, Fe.cr), pH = pH, O2 = O2,
              T = T, stable = list(NULL, dFe$predominant))
diagram(mFeCu$A.species, add = TRUE, col = 2, col.names = 2, bold = TRUE)
col <- c("#FF8C00", rep(NA, 6))
diagram(mFeCu$A.species, add = TRUE, col = col, lwd = 2, col.names = col, bold = TRUE)
TP <- describe.property(c("T", "P"), c(T, "Psat"))
legend("topright", TP, bty = "n")
title(paste("Cu-Fe-S-O-H; Total S =", 10^logaH2S, "m"))

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
#  diagram(mFe$A.bases, lty = 2, col = 4, col.names = 4, italic = TRUE)
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
#  srt <- function(a) ifelse(a$species$name %in% c("Mt+Bn", "Mt+Cct", "Mt+Ccp"), 80, 0)
#  
#  a11 <- mix(dFe, dCu, dFeCu)
#  diagram(a11, fill = fill(a11), srt = srt(a11), min.area = 0.01)
#  title("Fe:Cu = 1:1")
#  
#  a21 <- mix(dFe, dCu, dFeCu, c(2, 1))
#  diagram(a21, fill = fill(a21), srt = srt(a21), min.area = 0.01)
#  title("Fe:Cu = 2:1")
#  
#  a12 <- mix(dFe, dCu, dFeCu, c(1, 2))
#  diagram(a12, fill = fill(a12), srt = srt(a12), min.area = 0.01)
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
diagram(mFe$A.bases, lty = 2, col = 4, col.names = 4, italic = TRUE)
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
srt <- function(a) ifelse(a$species$name %in% c("Mt+Bn", "Mt+Cct", "Mt+Ccp"), 80, 0)

a11 <- mix(dFe, dCu, dFeCu)
diagram(a11, fill = fill(a11), srt = srt(a11), min.area = 0.01)
title("Fe:Cu = 1:1")

a21 <- mix(dFe, dCu, dFeCu, c(2, 1))
diagram(a21, fill = fill(a21), srt = srt(a21), min.area = 0.01)
title("Fe:Cu = 2:1")

a12 <- mix(dFe, dCu, dFeCu, c(1, 2))
diagram(a12, fill = fill(a12), srt = srt(a12), min.area = 0.01)
title("Fe:Cu = 1:2")

## ----stack2, eval = FALSE-----------------------------------------------------
#  T <- 125
#  layout(matrix(c(1, 2, 3, 3), nrow = 2), widths = c(1, 1.5))
#  
#  # Fe-S-O-H diagram
#  basis(c("copper", "hematite", "S2", "oxygen", "H2O", "H+", "Cl-"))
#  bFe <- species(c("hematite", "magnetite", "pyrite"))$name
#  aFe <- affinity(S2 = c(-34, -10), O2 = c(-55, -40), T = T)
#  # Order species by a function of composition to make colors
#  oFe <- order(aFe$species$S2 - aFe$species$O2)
#  fill <- terrain.colors(length(oFe), alpha = 0.3)[oFe]
#  abbrv <- info(aFe$species$ispecies)$abbrv
#  dFe <- diagram(aFe, names = abbrv, fill = fill)
#  title("Fe-S-O-H")
#  
#  # Cu-Fe-S-O-H diagram based on reactions with the
#  # stable Fe-bearing minerals (mosaic stack)
#  bCu <- species(c("copper", "chalcocite", "covellite", "chalcopyrite", "bornite"))$name
#  mCu <- mosaic(bFe, S2 = c(-34, -10), O2 = c(-55, -40),
#                T = T, stable = list(dFe$predominant))
#  oCu <- order(mCu$A.species$species$S2 - mCu$A.species$species$O2)
#  fill <- terrain.colors(length(oCu), alpha = 0.3)[oCu]
#  abbrv <- info(mCu$A.species$species$ispecies)$abbrv
#  dCu <- diagram(mCu$A.species, names = abbrv, fill = fill)
#  title("Cu-Fe-S-O-H")
#  
#  # Mash the diagrams and adjust labels
#  aFeCu <- mash(dFe, dCu)
#  names <- aFeCu$species$name
#  srt <- rep(0, length(names))
#  srt[names %in% c("Mt+Cu", "Hm+Cu")] <- 90
#  srt[names %in% c("Mt+Bn", "Hm+Bn")] <- 63
#  srt[names %in% c("Mt+Ccp", "Hm+Ccp")] <- 66
#  srt[names %in% c("Py+Bn", "Py+Cv")] <- 90
#  cex <- rep(1, length(names))
#  cex[names %in% c("Hm+Ccp", "Hm+Cv")] <- 0.8
#  oFeCu <- order(aFeCu$species$S2 - aFeCu$species$O2)
#  fill <- terrain.colors(length(oFeCu), alpha = 0.3)[oFeCu]
#  diagram(aFeCu, cex.names = cex, srt = srt, fill = fill)
#  legend("topleft", legend = lTP(T, "Psat"), bg = "white")
#  title("Cu-Fe-S-O-H")

## ----stack2, echo = FALSE, results = "hide", message = FALSE, fig.width = 6, fig.height = 4, out.width = "100%", pngquant = pngquant----
T <- 125
layout(matrix(c(1, 2, 3, 3), nrow = 2), widths = c(1, 1.5))

# Fe-S-O-H diagram
basis(c("copper", "hematite", "S2", "oxygen", "H2O", "H+", "Cl-"))
bFe <- species(c("hematite", "magnetite", "pyrite"))$name
aFe <- affinity(S2 = c(-34, -10), O2 = c(-55, -40), T = T)
# Order species by a function of composition to make colors
oFe <- order(aFe$species$S2 - aFe$species$O2)
fill <- terrain.colors(length(oFe), alpha = 0.3)[oFe]
abbrv <- info(aFe$species$ispecies)$abbrv
dFe <- diagram(aFe, names = abbrv, fill = fill)
title("Fe-S-O-H")

# Cu-Fe-S-O-H diagram based on reactions with the
# stable Fe-bearing minerals (mosaic stack)
bCu <- species(c("copper", "chalcocite", "covellite", "chalcopyrite", "bornite"))$name
mCu <- mosaic(bFe, S2 = c(-34, -10), O2 = c(-55, -40),
              T = T, stable = list(dFe$predominant))
oCu <- order(mCu$A.species$species$S2 - mCu$A.species$species$O2)
fill <- terrain.colors(length(oCu), alpha = 0.3)[oCu]
abbrv <- info(mCu$A.species$species$ispecies)$abbrv
dCu <- diagram(mCu$A.species, names = abbrv, fill = fill)
title("Cu-Fe-S-O-H")

# Mash the diagrams and adjust labels
aFeCu <- mash(dFe, dCu)
names <- aFeCu$species$name
srt <- rep(0, length(names))
srt[names %in% c("Mt+Cu", "Hm+Cu")] <- 90
srt[names %in% c("Mt+Bn", "Hm+Bn")] <- 63
srt[names %in% c("Mt+Ccp", "Hm+Ccp")] <- 66
srt[names %in% c("Py+Bn", "Py+Cv")] <- 90
cex <- rep(1, length(names))
cex[names %in% c("Hm+Ccp", "Hm+Cv")] <- 0.8
oFeCu <- order(aFeCu$species$S2 - aFeCu$species$O2)
fill <- terrain.colors(length(oFeCu), alpha = 0.3)[oFeCu]
diagram(aFeCu, cex.names = cex, srt = srt, fill = fill)
legend("topleft", legend = lTP(T, "Psat"), bg = "white")
title("Cu-Fe-S-O-H")

## ----solubility, eval = FALSE-------------------------------------------------
#  # Set up plot and system with aqueous Cu species
#  par(mfrow = c(1, 3))
#  basis("pH", 6)
#  iCu.aq <- retrieve("Cu", c("O", "H", "Cl", "S"), "aq")
#  species(iCu.aq)
#  # Estimate the molality of Cl for ca. 80,000 mg/kg solution (Table 2 of Sverjensky, 1987)
#  m_tot <- 80000 / mass("Cl") / 1000
#  calc <- NaCl(T = T, m_tot = m_tot)
#  # Use log molality here, not log activity, because
#  # activity coefficients are calculated by setting IS below
#  basis("Cl-", log10(calc$m_Cl))
#  
#  # Calculate affinities for aqueous Cu species while changing both Fe and Cu minerals
#  mfun <- function() {
#    mFeCu <- mosaic(list(bFe, bCu), S2 = c(-34, -10), O2 = c(-55, -40),
#      T = T, IS = calc$IS, stable = list(dFe$predominant, dCu$predominant))
#    # Calculate concentration of Cu
#    s <- solubility(mFeCu$A.species)
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
#  dG <- convert(dlogK, "G", T = convert(T, "K"))
#  # Apply this difference to the CuCl2- entry in OBIGT
#  newG <- info(info("CuCl2-"))$G + dG
#  mod.OBIGT("CuCl2-", G = newG)
#  
#  # Do the same thing for CuCl3-2
#  logK <- subcrt(c("CuCl3-2", "Cu+", "Cl-"), c(-1, 1, 3), T = T)$out$logK
#  dlogK <- logK - -5.6
#  dG <- convert(dlogK, "G", T = convert(T, "K"))
#  newG <- info(info("CuCl3-2"))$G + dG
#  mod.OBIGT("CuCl3-2", G = newG)
#  
#  # DIAGRAM 2
#  mfun()
#  title("Cu (ppm)", line = 1.7)
#  CuCl2 <- expr.species("CuCl2-")
#  CuCl3 <- expr.species("CuCl3-2")
#  title(bquote("Helgeson (1969)"~.(CuCl2)~and~.(CuCl3)), line = 0.9)
#  
#  # Set up system with SO4-2 (to dissolve S2(gas))
#  species("SO4-2")
#  aSO4 <- affinity(S2 = c(-34, -10), O2 = c(-55, -40), T = T, IS = calc$IS)
#  # Calculate concentration of SO4-2
#  s <- solubility(aSO4, in.terms.of = "SO4-2")
#  s <- convert(s, "ppm")
#  # DIAGRAM 3
#  diagram(aFeCu, names = NA, col = "gray", fill = fill)
#  diagram(s, type = "loga.balance", levels = 10^(-3:3), add = TRUE)
#  diagram(s, type = "loga.balance", levels = 35, add = TRUE, lwd = 3, col = 6, contour.method = NA)
#  title(bquote(bold(.(expr.species("SO4-2"))~"(ppm)")))

## ----solubility, echo = FALSE, results = "hide", message = FALSE, fig.width = 7, fig.height = 3, out.width = "100%", fig.align = "center", pngquant = pngquant----
# Set up plot and system with aqueous Cu species
par(mfrow = c(1, 3))
basis("pH", 6)
iCu.aq <- retrieve("Cu", c("O", "H", "Cl", "S"), "aq")
species(iCu.aq)
# Estimate the molality of Cl for ca. 80,000 mg/kg solution (Table 2 of Sverjensky, 1987)
m_tot <- 80000 / mass("Cl") / 1000
calc <- NaCl(T = T, m_tot = m_tot)
# Use log molality here, not log activity, because
# activity coefficients are calculated by setting IS below
basis("Cl-", log10(calc$m_Cl))

# Calculate affinities for aqueous Cu species while changing both Fe and Cu minerals
mfun <- function() {
  mFeCu <- mosaic(list(bFe, bCu), S2 = c(-34, -10), O2 = c(-55, -40),
    T = T, IS = calc$IS, stable = list(dFe$predominant, dCu$predominant))
  # Calculate concentration of Cu
  s <- solubility(mFeCu$A.species)
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
dG <- convert(dlogK, "G", T = convert(T, "K"))
# Apply this difference to the CuCl2- entry in OBIGT
newG <- info(info("CuCl2-"))$G + dG
mod.OBIGT("CuCl2-", G = newG)

# Do the same thing for CuCl3-2
logK <- subcrt(c("CuCl3-2", "Cu+", "Cl-"), c(-1, 1, 3), T = T)$out$logK
dlogK <- logK - -5.6
dG <- convert(dlogK, "G", T = convert(T, "K"))
newG <- info(info("CuCl3-2"))$G + dG
mod.OBIGT("CuCl3-2", G = newG)

# DIAGRAM 2
mfun()
title("Cu (ppm)", line = 1.7)
CuCl2 <- expr.species("CuCl2-")
CuCl3 <- expr.species("CuCl3-2")
title(bquote("Helgeson (1969)"~.(CuCl2)~and~.(CuCl3)), line = 0.9)

# Set up system with SO4-2 (to dissolve S2(gas))
species("SO4-2")
aSO4 <- affinity(S2 = c(-34, -10), O2 = c(-55, -40), T = T, IS = calc$IS)
# Calculate concentration of SO4-2
s <- solubility(aSO4, in.terms.of = "SO4-2")
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
#  dFe <- diagram(aFe, xlab = xlab)
#  title(paste("Only Fe; 1° balance:", dFe$balance))
#  label.plot("A")
#  
#  # Only Cu-bearing minerals
#  species(c("covellite", "chalcocite", "tenorite", "cuprite"))
#  aCu <- affinity(aFe)  # argument recall
#  dCu <- diagram(aCu, xlab = xlab)
#  title(paste("Only Cu; 1° balance:", dCu$balance))
#  label.plot("B")
#  
#  # Only Fe- AND Cu-bearing minerals
#  species(c("chalcopyrite", "bornite"))
#  aFeCu <- affinity(aFe)
#  dFeCu <- diagram(aFeCu, xlab = xlab, balance = "H+")
#  title(paste("Only Fe+Cu; 1° balance:", dFeCu$balance))
#  label.plot("C")
#  
#  ### SECONDARY balancing
#  
#  # Fe- or Cu-bearing minerals
#  ad1 <- rebalance(dFe, dCu, balance = "H+")
#  d1 <- diagram(ad1, xlab = xlab, balance = 1)
#  title("Only Fe or Cu; 2° balance: H+")
#  label.plot("D")
#  
#  # All minerals
#  d1$values <- c(dFe$values, dCu$values)
#  ad2 <- rebalance(d1, dFeCu, balance = "H+")
#  abbrv <- info(ad2$species$ispecies)$abbrv
#  diagram(ad2, xlab = xlab, balance = 1, names = abbrv)
#  title("Fe and/or Cu; 2° balance: H+")
#  label.plot("E")
#  
#  db <- describe.basis(ibasis = 3)
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
#  dFe <- diagram(aFe, xlab = xlab)
#  title(paste("Only Fe; 1° balance:", dFe$balance))
#  label.plot("A")
#  
#  # Only Cu-bearing minerals
#  species(c("covellite", "chalcocite", "tenorite", "cuprite"))
#  aCu <- affinity(aFe)  # argument recall
#  dCu <- diagram(aCu, xlab = xlab)
#  title(paste("Only Cu; 1° balance:", dCu$balance))
#  label.plot("B")
#  
#  # Only Fe- AND Cu-bearing minerals
#  species(c("chalcopyrite", "bornite"))
#  aFeCu <- affinity(aFe)
#  dFeCu <- diagram(aFeCu, xlab = xlab, balance = "H+")
#  title(paste("Only Fe+Cu; 1° balance:", dFeCu$balance))
#  label.plot("C")
#  
#  ### SECONDARY balancing
#  
#  # Fe- or Cu-bearing minerals
#  ad1 <- rebalance(dFe, dCu, balance = "H+")
#  d1 <- diagram(ad1, xlab = xlab, balance = 1)
#  title("Only Fe or Cu; 2° balance: H+")
#  label.plot("D")
#  
#  # All minerals
#  d1$values <- c(dFe$values, dCu$values)
#  ad2 <- rebalance(d1, dFeCu, balance = "H+")
#  abbrv <- info(ad2$species$ispecies)$abbrv
#  diagram(ad2, xlab = xlab, balance = 1, names = abbrv)
#  title("Fe and/or Cu; 2° balance: H+")
#  label.plot("E")
#  
#  db <- describe.basis(ibasis = 3)
#  leg <- lex(lTP(400, 2000), db)
#  legend("bottomleft", legend = leg, bty = "n")

## ----rebalance, eval = FALSE, echo = 10:30------------------------------------
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
#  dFe <- diagram(aFe, xlab = xlab)
#  title(paste("Only Fe; 1° balance:", dFe$balance))
#  label.plot("A")
#  
#  # Only Cu-bearing minerals
#  species(c("covellite", "chalcocite", "tenorite", "cuprite"))
#  aCu <- affinity(aFe)  # argument recall
#  dCu <- diagram(aCu, xlab = xlab)
#  title(paste("Only Cu; 1° balance:", dCu$balance))
#  label.plot("B")
#  
#  # Only Fe- AND Cu-bearing minerals
#  species(c("chalcopyrite", "bornite"))
#  aFeCu <- affinity(aFe)
#  dFeCu <- diagram(aFeCu, xlab = xlab, balance = "H+")
#  title(paste("Only Fe+Cu; 1° balance:", dFeCu$balance))
#  label.plot("C")
#  
#  ### SECONDARY balancing
#  
#  # Fe- or Cu-bearing minerals
#  ad1 <- rebalance(dFe, dCu, balance = "H+")
#  d1 <- diagram(ad1, xlab = xlab, balance = 1)
#  title("Only Fe or Cu; 2° balance: H+")
#  label.plot("D")
#  
#  # All minerals
#  d1$values <- c(dFe$values, dCu$values)
#  ad2 <- rebalance(d1, dFeCu, balance = "H+")
#  abbrv <- info(ad2$species$ispecies)$abbrv
#  diagram(ad2, xlab = xlab, balance = 1, names = abbrv)
#  title("Fe and/or Cu; 2° balance: H+")
#  label.plot("E")
#  
#  db <- describe.basis(ibasis = 3)
#  leg <- lex(lTP(400, 2000), db)
#  legend("bottomleft", legend = leg, bty = "n")

## ----rebalance, eval = FALSE, echo = 33:43------------------------------------
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
#  dFe <- diagram(aFe, xlab = xlab)
#  title(paste("Only Fe; 1° balance:", dFe$balance))
#  label.plot("A")
#  
#  # Only Cu-bearing minerals
#  species(c("covellite", "chalcocite", "tenorite", "cuprite"))
#  aCu <- affinity(aFe)  # argument recall
#  dCu <- diagram(aCu, xlab = xlab)
#  title(paste("Only Cu; 1° balance:", dCu$balance))
#  label.plot("B")
#  
#  # Only Fe- AND Cu-bearing minerals
#  species(c("chalcopyrite", "bornite"))
#  aFeCu <- affinity(aFe)
#  dFeCu <- diagram(aFeCu, xlab = xlab, balance = "H+")
#  title(paste("Only Fe+Cu; 1° balance:", dFeCu$balance))
#  label.plot("C")
#  
#  ### SECONDARY balancing
#  
#  # Fe- or Cu-bearing minerals
#  ad1 <- rebalance(dFe, dCu, balance = "H+")
#  d1 <- diagram(ad1, xlab = xlab, balance = 1)
#  title("Only Fe or Cu; 2° balance: H+")
#  label.plot("D")
#  
#  # All minerals
#  d1$values <- c(dFe$values, dCu$values)
#  ad2 <- rebalance(d1, dFeCu, balance = "H+")
#  abbrv <- info(ad2$species$ispecies)$abbrv
#  diagram(ad2, xlab = xlab, balance = 1, names = abbrv)
#  title("Fe and/or Cu; 2° balance: H+")
#  label.plot("E")
#  
#  db <- describe.basis(ibasis = 3)
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
dFe <- diagram(aFe, xlab = xlab)
title(paste("Only Fe; 1° balance:", dFe$balance))
label.plot("A")

# Only Cu-bearing minerals
species(c("covellite", "chalcocite", "tenorite", "cuprite"))
aCu <- affinity(aFe)  # argument recall
dCu <- diagram(aCu, xlab = xlab)
title(paste("Only Cu; 1° balance:", dCu$balance))
label.plot("B")

# Only Fe- AND Cu-bearing minerals
species(c("chalcopyrite", "bornite"))
aFeCu <- affinity(aFe)
dFeCu <- diagram(aFeCu, xlab = xlab, balance = "H+")
title(paste("Only Fe+Cu; 1° balance:", dFeCu$balance))
label.plot("C")

### SECONDARY balancing

# Fe- or Cu-bearing minerals
ad1 <- rebalance(dFe, dCu, balance = "H+")
d1 <- diagram(ad1, xlab = xlab, balance = 1)
title("Only Fe or Cu; 2° balance: H+")
label.plot("D")

# All minerals
d1$values <- c(dFe$values, dCu$values)
ad2 <- rebalance(d1, dFeCu, balance = "H+")
abbrv <- info(ad2$species$ispecies)$abbrv
diagram(ad2, xlab = xlab, balance = 1, names = abbrv)
title("Fe and/or Cu; 2° balance: H+")
label.plot("E")

db <- describe.basis(ibasis = 3)
leg <- lex(lTP(400, 2000), db)
legend("bottomleft", legend = leg, bty = "n")

## ----non-metal, results = "hide", message = FALSE, fig.width = 6, fig.height = 5, out.width = "80%", fig.align = "center", pngquant = pngquant----
basis(c("Fe+2", "Cu+", "hydrogen sulfide", "oxygen", "H2O", "H+"))
basis("H2S", 2)
species(c("pyrite", "magnetite", "hematite", "covellite", "tenorite",
          "chalcopyrite", "bornite"))
a <- affinity("Cu+" = c(-8, 2, 500), "Fe+2" = c(-4, 12, 500), T = 400, P = 2000)
d <- diagram(a, xlab = ratlab("Cu+"), ylab = ratlab("Fe+2"), balance = "O2")
title(paste("Cu-Fe-S-O-H; 1° balance:", d$balance))
# Add saturation lines
species(c("pyrrhotite", "ferrous-oxide", "chalcocite", "cuprite"))
asat <- affinity(a)  # argument recall
diagram(asat, type = "saturation", add = TRUE, lty = 2, col = 4)
legend("topleft", legend = lTP(400, 2000), bty = "n")

## ----mosaic-combo, results = "hide", message = FALSE, fig.width = 6, fig.height = 5, out.width = "80%", fig.align = "center", pngquant = pngquant----
loga_N <- -4
loga_C <- -3
basis(c("CO2", "NH3", "O2", "H2O", "H+"))
basis("NH3", loga_N)
species(c("acetamide", "acetic acid", "acetate"))
m <- mosaic(c("NH3", "NH4+"), pH = c(0, 14))
e <- equilibrate(m, loga.balance = loga_C)
diagram(e, ylim = c(-10, -2))
tN <- paste("log(total N in basis species) =", loga_N)
tC <- paste("log(total C in formed species) =", loga_C)
title(main = paste(tN, tC, sep = "\n"), font.main = 1)
legend("topright", legend = lTP(25, 1), bty = "n")

