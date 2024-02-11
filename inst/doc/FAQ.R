## ----setup, include = FALSE---------------------------------------------------
library(CHNOSZ)
options(width = 80)
# Use pngquant to optimize PNG images
library(knitr)
knit_hooks$set(pngquant = hook_pngquant)
pngquant <- "--speed=1 --quality=0-25"
if (!nzchar(Sys.which("pngquant"))) pngquant <- NULL
# To make warnings appear in text box 20230619
# https://selbydavid.com/2017/06/18/rmarkdown-alerts/
knitr::knit_hooks$set(
   error = function(x, options) {
     paste('\n\n<div class="alert alert-danger">',
           gsub('##', '\n', gsub('^##\ Error', '**Error:**', x)),
           '</div>', sep = '\n')
   },
   warning = function(x, options) {
     paste('\n\n<div class="alert alert-warning">',
           gsub('##', '\n', gsub('^##\ Warning:', '**Warning:**', x)),
           '</div>', sep = '\n')
   },
   message = function(x, options) {
     paste('\n\n<div class="alert alert-info">',
           gsub('##', '\n', x),
           '</div>', sep = '\n')
   }
)

# Set dpi 20231129
knitr::opts_chunk$set(
  dpi = if(nzchar(Sys.getenv("CHNOSZ_BUILD_LARGE_VIGNETTES"))) 100 else 72
)

## ----echo = F, cache = F------------------------------------------------------
# Merge consecutive messages into a single div 20231114
knitr::knit_hooks$set(document = function(x){

  # Not sure why this is needed, but simply computing on 'x' doesn't work
  file <- tempfile()
  writeLines(x, file)
  x <- readLines(file)

  # Line numbers of the document with </div>
  enddiv <- which(x == "</div>")
  # Line numbers with <div class="alert alert-info">
  beginalert <- which(x == '<div class="alert alert-info">')
  # Find </div> followed <div class="alert alert-info"> (skip empty lines)
  removediv <- (enddiv + 2) %in% beginalert
  if(any(removediv)) {
    # Lines to remove
    rmlines1 <- enddiv[removediv]
    rmlines2 <- enddiv[removediv] + 1
    rmlines3 <- enddiv[removediv] + 2
    # Do the removing
    x <- x[-c(rmlines1, rmlines2, rmlines3)]
  }

  x

})

## ----HTML, include = FALSE----------------------------------------------------
NOTE <- '<span style="background-color: yellow;">NOTE</span>'
# CHNOSZ functions
equilibrate_ <- '<code>equilibrate()</code>'
info_ <- '<code>info()</code>'
thermo.refs_ <- '<code>thermo.refs()</code>'
# Math stuff
logK <- "log&thinsp;<i>K</i>"
Hplus <- "H<sup>+</sup>"
HCO2_ <- "HCO<sub>2</sub><sup>−</sup>"
HCO3_ <- "HCO<sub>3</sub><sup>−</sup>"
O2 <- "O<sub>2</sub>"
S2 <- "S<sub>2</sub>"
log <- "log&thinsp;"
aHCO2_ <- "<i>a</i><sub>HCO<sub>2</sub><sup>−</sup></sub>"
aHCO3_ <- "<i>a</i><sub>HCO<sub>3</sub><sup>−</sup></sub>"
logfO2 <- "log&thinsp;<i>f</i><sub>O<sub>2</sub></sub>"
Ctot <- "C<sub>tot</sub>"
C3H5O2_ <- "C<sub>3</sub>H<sub>5</sub>O<sub>2</sub><sup>−</sup>"
a3HCO3_ <- "<i>a</i><sup>3</sup><sub>HCO<sub>3</sub><sup>−</sup></sub>"
aC3H5O2_ <- "<i>a</i><sub>C<sub>3</sub>H<sub>5</sub>O<sub>2</sub><sup>−</sup></sub>"
a2HCO3_ <- "<i>a</i><sup>2</sup><sub>HCO<sub>3</sub><sup>−</sup></sub>"
logCtot <- "log&thinsp;C<sub>tot</sub>"
CO2 <- "CO<sub>2</sub>"
H2O <- "H<sub>2</sub>O"
S3minus <- "S<sub>3</sub><sup>-</sup>"
H2S <- "H<sub>2</sub>S"
SO4__ <- "SO<sub>4</sub><sup>-2</sup>"
Kplus <- "K<sup>+</sup>"
Naplus <- "Na<sup>+</sup>"
Clminus <- "Cl<sup>-</sup>"
H2 <- "H<sub>2</sub>"

## ----alanine_refs, message = FALSE--------------------------------------------
info(info("alanine"))[c("ref1", "ref2")]
thermo.refs(info("alanine"))

## ----PPM_refs, message = FALSE------------------------------------------------
basis(c("pyrite", "pyrrhotite", "oxygen"))
sres <- subcrt("magnetite", 1)
info(sres$reaction$ispecies)[, 1:6]
thermo.refs(sres)
reset()

## ----DEW_Ctot, eval = FALSE---------------------------------------------------
#  ### Based on demo/DEW.R in CHNOSZ
#  
#  # Set up subplots
#  par(mfrow = c(1, 2), mar = c(3.0, 3.5, 2.5, 1.0), mgp = c(1.7, 0.3, 0), las = 1, tcl = 0.3, xaxs = "i", yaxs = "i")
#  
#  # Activate DEW model
#  water("DEW")
#  
#  ###########
#  ## logfO2-pH diagram for aqueous inorganic and organic carbon species at high pressure
#  ## After Figure 1b of Sverjensky et al., 2014b [SSH14]
#  ## (Nature Geoscience, https://doi.org/10.1038/NGEO2291)
#  ###########
#  
#  # Define system
#  basis("CHNOS+")
#  inorganic <- c("CO2", "HCO3-", "CO3-2", "CH4")
#  organic <- c("formic acid", "formate", "acetic acid", "acetate", "propanoic acid", "propanoate")
#  species(c(inorganic, organic), 0)
#  fill <- c(rep("aliceblue", length(inorganic)), rep("aquamarine", length(organic)))
#  
#  # A function to make the diagrams
#  dfun <- function(T = 600, P = 50000, res = 200, Ctot = 0.03) {
#    a <- affinity(pH = c(0, 10, res), O2 = c(-24, -12, res), T = T, P = P)
#    # Set total C concentration to 0.03 molal
#    # (from EQ3NR model for eclogite [Supporting Information of SSH14])
#    e <- equilibrate(a, loga.balance = log10(Ctot))
#    diagram(e, fill = fill)
#    dp <- describe.property(c("     T", "     P"), c(T, P), digits = 0)
#    legend("bottomleft", legend = dp, bty = "n")
#  }
#  
#  water("DEW")
#  add.OBIGT("DEW")
#  dfun(Ctot = 0.03)
#  mtitle(c("Inorganic and organic species", "C[total] = 0.03 molal"))
#  dfun(Ctot = 3)
#  mtitle(c("Inorganic and organic species", "C[total] = 3 molal"))
#  
#  # Restore default settings for the questions below
#  reset()

## ----DEW_Ctot, echo = FALSE, message = FALSE, results = "hide", fig.width = 8, fig.height = 4, out.width = "100%", fig.align = "center", pngquant = pngquant, cache = TRUE----
### Based on demo/DEW.R in CHNOSZ

# Set up subplots
par(mfrow = c(1, 2), mar = c(3.0, 3.5, 2.5, 1.0), mgp = c(1.7, 0.3, 0), las = 1, tcl = 0.3, xaxs = "i", yaxs = "i")

# Activate DEW model
water("DEW")

###########
## logfO2-pH diagram for aqueous inorganic and organic carbon species at high pressure
## After Figure 1b of Sverjensky et al., 2014b [SSH14]
## (Nature Geoscience, https://doi.org/10.1038/NGEO2291)
###########

# Define system
basis("CHNOS+")
inorganic <- c("CO2", "HCO3-", "CO3-2", "CH4")
organic <- c("formic acid", "formate", "acetic acid", "acetate", "propanoic acid", "propanoate")
species(c(inorganic, organic), 0)
fill <- c(rep("aliceblue", length(inorganic)), rep("aquamarine", length(organic)))

# A function to make the diagrams
dfun <- function(T = 600, P = 50000, res = 200, Ctot = 0.03) {
  a <- affinity(pH = c(0, 10, res), O2 = c(-24, -12, res), T = T, P = P)
  # Set total C concentration to 0.03 molal
  # (from EQ3NR model for eclogite [Supporting Information of SSH14])
  e <- equilibrate(a, loga.balance = log10(Ctot))
  diagram(e, fill = fill)
  dp <- describe.property(c("     T", "     P"), c(T, P), digits = 0)
  legend("bottomleft", legend = dp, bty = "n")
}

water("DEW")
add.OBIGT("DEW")
dfun(Ctot = 0.03)
mtitle(c("Inorganic and organic species", "C[total] = 0.03 molal"))
dfun(Ctot = 3)
mtitle(c("Inorganic and organic species", "C[total] = 3 molal"))

# Restore default settings for the questions below
reset()

## ----pyrrhotite_values, message = FALSE---------------------------------------
# The formula of the new mineral and literature reference
formula <- "Fe0.877S"
ref1 <- "PMW87"
# Because the properties from Pankratz et al. (1987) are listed in calories,
# we need to change the output of subcrt() to calories (the default is Joules)
E.units("cal")
# Use temperature in Kelvin for the calculations below
T.units("K")
# Thermodynamic properties of polymorph 1 at 25 °C (298.15 K)
G1 <- -25543
H1 <- -25200
S1 <- 14.531
Cp1 <- 11.922
# Heat capacity coefficients for polymorph 1
a1 <- 7.510
b1 <- 0.014798
# For volume, use the value from Helgeson et al. (1978)
V1 <- V2 <- 18.2
# Transition temperature
Ttr <- 598
# Transition enthalpy (cal/mol)
DHtr <- 95
# Heat capacity coefficients for polymorph 2
a2 <- -1.709
b2 <- 0.011746
c2 <- 3073400
# Maximum temperature of polymorph 2
T2 <- 1800

## ----pyrrhotite_Ttr, message = FALSE------------------------------------------
DGtr <- 0  # DON'T CHANGE THIS
TDStr <- DHtr - DGtr  # TΔS° = ΔH° - ΔG°
DStr <- TDStr / Ttr

## ----pyrrhotite_Cp, results = "hide", message = FALSE-------------------------
mod.OBIGT("pyrrhotite_new", formula = formula, state = "cr", ref1 = ref1,
  E_units = "cal", G = 0, H = 0, S = 0, V = V1, Cp = Cp1,
  a = a1, b = b1, c = 0, d = 0, e = 0, f = 0, lambda = 0, T = Ttr)
mod.OBIGT("pyrrhotite_new", formula = formula, state = "cr2", ref1 = ref1,
  E_units = "cal", G = 0, H = 0, S = 0, V = V2,
  a = a2, b = b2, c = c2, d = 0, e = 0, f = 0, lambda = 0, T = T2)

## ----pyrrhotite_info, results = "hide", collapse = TRUE-----------------------
info(info("pyrrhotite_new", "cr"))
info(info("pyrrhotite_new", "cr2"))

## ----pyrrhotite_S, message = FALSE--------------------------------------------
DS1 <- subcrt("pyrrhotite_new", "cr", P = 1, T = Ttr, use.polymorphs = FALSE)$out[[1]]$S
DS2 <- subcrt("pyrrhotite_new", "cr2", P = 1, T = Ttr)$out[[1]]$S
DS298 <- DS1 + DStr - DS2

## ----pyrrhotite_D1_D2, message = FALSE, results = "hide"----------------------
mod.OBIGT("pyrrhotite_new", state = "cr", S = S1)
mod.OBIGT("pyrrhotite_new", state = "cr2", S = S1 + DS298)
D1 <- subcrt("pyrrhotite_new", "cr", P = 1, T = Ttr, use.polymorphs = FALSE)$out[[1]]
D2 <- subcrt("pyrrhotite_new", "cr2", P = 1, T = Ttr)$out[[1]]

## ----pyrrhotite_check_S, results = "hide"-------------------------------------
stopifnot(all.equal(D2$S - D1$S, DStr))

## ----pyrrhotite_DG298_DH298, results = "hide", message = FALSE----------------
DG298 <- D1$G + DGtr - D2$G
DH298 <- D1$H + DHtr - D2$H
mod.OBIGT("pyrrhotite_new", state = "cr", G = G1, H = H1)
mod.OBIGT("pyrrhotite_new", state = "cr2", G = G1 + DG298, H = H1 + DH298)

## ----pyrrhotite_info2, collapse = TRUE----------------------------------------
cr_parameters <- info(info("pyrrhotite_new", "cr"))
stopifnot(abs(check.GHS(cr_parameters)) < 1)
cr2_parameters <- info(info("pyrrhotite_new", "cr2"))
stopifnot(abs(check.GHS(cr2_parameters)) < 1)

## ----pyrrhotite_parameters----------------------------------------------------
cr_parameters
cr2_parameters

## ----pyrrhotite_T, eval = FALSE-----------------------------------------------
#  T <- seq(550, 650, 1)
#  sout <- subcrt("pyrrhotite_new", T = T, P = 1)$out[[1]]
#  par(mfrow = c(1, 4), mar = c(4, 4.2, 1, 1))
#  labels <- c(G = "DG0f", H = "DH0f", S = "S0", Cp = "Cp0")
#  for(property in c("G", "H", "S", "Cp")) {
#    plot(sout$T, sout[, property], col = sout$polymorph,
#      xlab = axis.label("T"), ylab = axis.label(labels[property])
#    )
#    abline(v = Ttr, lty = 3, col = 8)
#    if(property == "G")
#      legend("bottomleft", c("1", "2"), pch = 1, col = c(1, 2), title = "Polymorph")
#  }
#  
#  # Restore default settings for the questions below
#  reset()

## ----pyrrhotite_T, echo = FALSE, message = FALSE, results = "hide", fig.width = 8, fig.height = 2.5, out.width = "100%", fig.align = "center", pngquant = pngquant----
T <- seq(550, 650, 1)
sout <- subcrt("pyrrhotite_new", T = T, P = 1)$out[[1]]
par(mfrow = c(1, 4), mar = c(4, 4.2, 1, 1))
labels <- c(G = "DG0f", H = "DH0f", S = "S0", Cp = "Cp0")
for(property in c("G", "H", "S", "Cp")) {
  plot(sout$T, sout[, property], col = sout$polymorph,
    xlab = axis.label("T"), ylab = axis.label(labels[property])
  )
  abline(v = Ttr, lty = 3, col = 8)
  if(property == "G")
    legend("bottomleft", c("1", "2"), pch = 1, col = c(1, 2), title = "Polymorph")
}

# Restore default settings for the questions below
reset()

## ----trisulfur, eval = FALSE, echo = FALSE------------------------------------
#  par(mfrow = c(1, 3))
#  
#  ## BLOCK 1
#  T <- 350
#  P <- 5000
#  res <- 200
#  
#  ## BLOCK 2
#  wt_percent_S <- 10
#  wt_permil_S <- wt_percent_S * 10
#  molar_mass_S <- mass("S") # 32.06
#  moles_S <- wt_permil_S / molar_mass_S
#  grams_H2O <- 1000 - wt_permil_S
#  molality_S <- 1000 * moles_S / grams_H2O
#  logm_S <- log10(molality_S)
#  
#  ## BLOCK 3
#  basis(c("Ni", "SiO2", "Fe2O3", "H2S", "H2O", "oxygen", "H+"))
#  species(c("H2S", "HS-", "SO2", "HSO4-", "SO4-2", "S3-"))
#  a <- affinity(pH = c(2, 10, res), O2 = c(-34, -22, res), T = T, P = P)
#  e <- equilibrate(a, loga.balance = logm_S)
#  diagram(e)
#  
#  ## BLOCK 4
#  mod.buffer("NNO", c("nickel", "bunsenite"), state = "cr", logact = 0)
#  for(buffer in c("HM", "QFM", "NNO")) {
#    basis("O2", buffer)
#    logfO2_ <- affinity(return.buffer = TRUE, T = T, P = P)$O2
#    abline(h = logfO2_, lty = 3, col = 2)
#    text(8.5, logfO2_, buffer, adj = c(0, 0), col = 2, cex = 0.9)
#  }
#  
#  ## BLOCK 5
#  pH <- subcrt(c("H2O", "H+", "OH-"), c(-1, 1, 1), T = T, P = P)$out$logK / -2
#  abline(v = pH, lty = 2, col = 4)
#  
#  ## BLOCK 6
#  lTP <- describe.property(c("T", "P"), c(T, P))
#  lS <- paste(wt_percent_S, "wt% S(aq)")
#  ltext <- c(lTP, lS)
#  legend("topright", ltext, bty = "n", bg = "transparent")
#  title(quote("Parameters for"~S[3]^"-"~"from Pokrovski and Dubessy (2015)"), xpd = NA)
#  
#  ######## Plot 2: Modify Gibbs energy
#  
#  oldG <- info(info("S3-"))$G
#  newG <- oldG - 12548
#  mod.OBIGT("S3-", G = newG)
#  a <- affinity(pH = c(2, 10, res), O2 = c(-34, -22, res), T = T, P = P)
#  e <- equilibrate(a, loga.balance = logm_S)
#  diagram(e)
#  legend("topright", ltext, bty = "n", bg = "transparent")
#  title(quote("Modified"~log*italic(K)~"after Pokrovski and Dubrovinsky (2011)"), xpd = NA)
#  OBIGT()
#  
#  ######## Plot 3: Do it with DEW
#  
#  T <- 700
#  P <- 10000 # 10000 bar = 1 GPa
#  oldwat <- water("DEW")
#  add.OBIGT("DEW")
#  info(species()$ispecies)
#  a <- affinity(pH = c(2, 10, res), O2 = c(-18, -10, res), T = T, P = P)
#  e <- equilibrate(a, loga.balance = logm_S)
#  diagram(e)
#  lTP <- describe.property(c("T", "P"), c(T, P))
#  ltext <- c(lTP, lS)
#  legend("topright", ltext, bty = "n", bg = "transparent")
#  title(quote("Deep Earth Water (DEW)"~"model"))
#  water(oldwat)
#  OBIGT()

## ----trisulfur, eval = FALSE, echo = 3:43-------------------------------------
#  par(mfrow = c(1, 3))
#  
#  ## BLOCK 1
#  T <- 350
#  P <- 5000
#  res <- 200
#  
#  ## BLOCK 2
#  wt_percent_S <- 10
#  wt_permil_S <- wt_percent_S * 10
#  molar_mass_S <- mass("S") # 32.06
#  moles_S <- wt_permil_S / molar_mass_S
#  grams_H2O <- 1000 - wt_permil_S
#  molality_S <- 1000 * moles_S / grams_H2O
#  logm_S <- log10(molality_S)
#  
#  ## BLOCK 3
#  basis(c("Ni", "SiO2", "Fe2O3", "H2S", "H2O", "oxygen", "H+"))
#  species(c("H2S", "HS-", "SO2", "HSO4-", "SO4-2", "S3-"))
#  a <- affinity(pH = c(2, 10, res), O2 = c(-34, -22, res), T = T, P = P)
#  e <- equilibrate(a, loga.balance = logm_S)
#  diagram(e)
#  
#  ## BLOCK 4
#  mod.buffer("NNO", c("nickel", "bunsenite"), state = "cr", logact = 0)
#  for(buffer in c("HM", "QFM", "NNO")) {
#    basis("O2", buffer)
#    logfO2_ <- affinity(return.buffer = TRUE, T = T, P = P)$O2
#    abline(h = logfO2_, lty = 3, col = 2)
#    text(8.5, logfO2_, buffer, adj = c(0, 0), col = 2, cex = 0.9)
#  }
#  
#  ## BLOCK 5
#  pH <- subcrt(c("H2O", "H+", "OH-"), c(-1, 1, 1), T = T, P = P)$out$logK / -2
#  abline(v = pH, lty = 2, col = 4)
#  
#  ## BLOCK 6
#  lTP <- describe.property(c("T", "P"), c(T, P))
#  lS <- paste(wt_percent_S, "wt% S(aq)")
#  ltext <- c(lTP, lS)
#  legend("topright", ltext, bty = "n", bg = "transparent")
#  title(quote("Parameters for"~S[3]^"-"~"from Pokrovski and Dubessy (2015)"), xpd = NA)
#  
#  ######## Plot 2: Modify Gibbs energy
#  
#  oldG <- info(info("S3-"))$G
#  newG <- oldG - 12548
#  mod.OBIGT("S3-", G = newG)
#  a <- affinity(pH = c(2, 10, res), O2 = c(-34, -22, res), T = T, P = P)
#  e <- equilibrate(a, loga.balance = logm_S)
#  diagram(e)
#  legend("topright", ltext, bty = "n", bg = "transparent")
#  title(quote("Modified"~log*italic(K)~"after Pokrovski and Dubrovinsky (2011)"), xpd = NA)
#  OBIGT()
#  
#  ######## Plot 3: Do it with DEW
#  
#  T <- 700
#  P <- 10000 # 10000 bar = 1 GPa
#  oldwat <- water("DEW")
#  add.OBIGT("DEW")
#  info(species()$ispecies)
#  a <- affinity(pH = c(2, 10, res), O2 = c(-18, -10, res), T = T, P = P)
#  e <- equilibrate(a, loga.balance = logm_S)
#  diagram(e)
#  lTP <- describe.property(c("T", "P"), c(T, P))
#  ltext <- c(lTP, lS)
#  legend("topright", ltext, bty = "n", bg = "transparent")
#  title(quote("Deep Earth Water (DEW)"~"model"))
#  water(oldwat)
#  OBIGT()

## ----trisulfur_logK, message = FALSE, echo = 1:3------------------------------
species <- c("H2S", "SO4-2", "H+", "S3-", "oxygen", "H2O")
coeffs <- c(-2, -1, -1, 1, 0.75, 2.5)
(calclogK <- subcrt(species, coeffs, T = seq(300, 450, 50), P = 5000)$out$logK)
fcalclogK <- formatC(calclogK, format = "f", digits = 1)
reflogK <- -9.6
dlogK <- calclogK[2] - reflogK
# Put in a test so that we don't get surprised by
# future database updates or changes to this vignette
stopifnot(round(dlogK, 4) == -4.4132)

## ----trisulfur, eval = FALSE, echo = 46:55------------------------------------
#  par(mfrow = c(1, 3))
#  
#  ## BLOCK 1
#  T <- 350
#  P <- 5000
#  res <- 200
#  
#  ## BLOCK 2
#  wt_percent_S <- 10
#  wt_permil_S <- wt_percent_S * 10
#  molar_mass_S <- mass("S") # 32.06
#  moles_S <- wt_permil_S / molar_mass_S
#  grams_H2O <- 1000 - wt_permil_S
#  molality_S <- 1000 * moles_S / grams_H2O
#  logm_S <- log10(molality_S)
#  
#  ## BLOCK 3
#  basis(c("Ni", "SiO2", "Fe2O3", "H2S", "H2O", "oxygen", "H+"))
#  species(c("H2S", "HS-", "SO2", "HSO4-", "SO4-2", "S3-"))
#  a <- affinity(pH = c(2, 10, res), O2 = c(-34, -22, res), T = T, P = P)
#  e <- equilibrate(a, loga.balance = logm_S)
#  diagram(e)
#  
#  ## BLOCK 4
#  mod.buffer("NNO", c("nickel", "bunsenite"), state = "cr", logact = 0)
#  for(buffer in c("HM", "QFM", "NNO")) {
#    basis("O2", buffer)
#    logfO2_ <- affinity(return.buffer = TRUE, T = T, P = P)$O2
#    abline(h = logfO2_, lty = 3, col = 2)
#    text(8.5, logfO2_, buffer, adj = c(0, 0), col = 2, cex = 0.9)
#  }
#  
#  ## BLOCK 5
#  pH <- subcrt(c("H2O", "H+", "OH-"), c(-1, 1, 1), T = T, P = P)$out$logK / -2
#  abline(v = pH, lty = 2, col = 4)
#  
#  ## BLOCK 6
#  lTP <- describe.property(c("T", "P"), c(T, P))
#  lS <- paste(wt_percent_S, "wt% S(aq)")
#  ltext <- c(lTP, lS)
#  legend("topright", ltext, bty = "n", bg = "transparent")
#  title(quote("Parameters for"~S[3]^"-"~"from Pokrovski and Dubessy (2015)"), xpd = NA)
#  
#  ######## Plot 2: Modify Gibbs energy
#  
#  oldG <- info(info("S3-"))$G
#  newG <- oldG - 12548
#  mod.OBIGT("S3-", G = newG)
#  a <- affinity(pH = c(2, 10, res), O2 = c(-34, -22, res), T = T, P = P)
#  e <- equilibrate(a, loga.balance = logm_S)
#  diagram(e)
#  legend("topright", ltext, bty = "n", bg = "transparent")
#  title(quote("Modified"~log*italic(K)~"after Pokrovski and Dubrovinsky (2011)"), xpd = NA)
#  OBIGT()
#  
#  ######## Plot 3: Do it with DEW
#  
#  T <- 700
#  P <- 10000 # 10000 bar = 1 GPa
#  oldwat <- water("DEW")
#  add.OBIGT("DEW")
#  info(species()$ispecies)
#  a <- affinity(pH = c(2, 10, res), O2 = c(-18, -10, res), T = T, P = P)
#  e <- equilibrate(a, loga.balance = logm_S)
#  diagram(e)
#  lTP <- describe.property(c("T", "P"), c(T, P))
#  ltext <- c(lTP, lS)
#  legend("topright", ltext, bty = "n", bg = "transparent")
#  title(quote("Deep Earth Water (DEW)"~"model"))
#  water(oldwat)
#  OBIGT()

## ----trisulfur_DEW, message = FALSE, echo = 8:22, collapse = TRUE, fig.keep = "none"----
# The first four lines are stand-ins to make this block runnable and get the right output from info();
# the diagram actually shown in the vignette is made using the 'trisulfur' block above
b <- basis("CHNOS+")
s <- species(c("H2S", "HS-", "SO2", "HSO4-", "SO4-2", "S3-"))
res <- 10
wt_percent_S <- logm_S <- 0
######## End of stand-in code ########
T <- 700
P <- 10000 # 10000 bar = 1 GPa
oldwat <- water("DEW")
add.OBIGT("DEW")
info(species()$ispecies)[, 1:8]
a <- affinity(pH = c(2, 10, res), O2 = c(-18, -10, res), T = T, P = P)
e <- equilibrate(a, loga.balance = logm_S)
diagram(e)
lTP <- describe.property(c("T", "P"), c(T, P))
lS <- paste(wt_percent_S, "wt% S(aq)")
ltext <- c(lTP, lS)
legend("topright", ltext, bty = "n", bg = "transparent")
title(quote("Deep Earth Water (DEW)"~"model"))
water(oldwat)
OBIGT()

## ----trisulfur, echo = FALSE, message = FALSE, results = "hide", fig.width = 10, fig.height = 3.33, out.width = "100%", out.extra='class="full-width"', pngquant = pngquant, cache = TRUE----
par(mfrow = c(1, 3))

## BLOCK 1
T <- 350
P <- 5000
res <- 200

## BLOCK 2
wt_percent_S <- 10
wt_permil_S <- wt_percent_S * 10
molar_mass_S <- mass("S") # 32.06
moles_S <- wt_permil_S / molar_mass_S
grams_H2O <- 1000 - wt_permil_S
molality_S <- 1000 * moles_S / grams_H2O
logm_S <- log10(molality_S)

## BLOCK 3
basis(c("Ni", "SiO2", "Fe2O3", "H2S", "H2O", "oxygen", "H+"))
species(c("H2S", "HS-", "SO2", "HSO4-", "SO4-2", "S3-"))
a <- affinity(pH = c(2, 10, res), O2 = c(-34, -22, res), T = T, P = P)
e <- equilibrate(a, loga.balance = logm_S)
diagram(e)

## BLOCK 4
mod.buffer("NNO", c("nickel", "bunsenite"), state = "cr", logact = 0)
for(buffer in c("HM", "QFM", "NNO")) {
  basis("O2", buffer)
  logfO2_ <- affinity(return.buffer = TRUE, T = T, P = P)$O2
  abline(h = logfO2_, lty = 3, col = 2)
  text(8.5, logfO2_, buffer, adj = c(0, 0), col = 2, cex = 0.9)
}

## BLOCK 5
pH <- subcrt(c("H2O", "H+", "OH-"), c(-1, 1, 1), T = T, P = P)$out$logK / -2
abline(v = pH, lty = 2, col = 4)

## BLOCK 6
lTP <- describe.property(c("T", "P"), c(T, P))
lS <- paste(wt_percent_S, "wt% S(aq)")
ltext <- c(lTP, lS)
legend("topright", ltext, bty = "n", bg = "transparent")
title(quote("Parameters for"~S[3]^"-"~"from Pokrovski and Dubessy (2015)"), xpd = NA)

######## Plot 2: Modify Gibbs energy

oldG <- info(info("S3-"))$G
newG <- oldG - 12548
mod.OBIGT("S3-", G = newG)
a <- affinity(pH = c(2, 10, res), O2 = c(-34, -22, res), T = T, P = P)
e <- equilibrate(a, loga.balance = logm_S)
diagram(e)
legend("topright", ltext, bty = "n", bg = "transparent")
title(quote("Modified"~log*italic(K)~"after Pokrovski and Dubrovinsky (2011)"), xpd = NA)
OBIGT()

######## Plot 3: Do it with DEW

T <- 700
P <- 10000 # 10000 bar = 1 GPa
oldwat <- water("DEW")
add.OBIGT("DEW")
info(species()$ispecies)
a <- affinity(pH = c(2, 10, res), O2 = c(-18, -10, res), T = T, P = P)
e <- equilibrate(a, loga.balance = logm_S)
diagram(e)
lTP <- describe.property(c("T", "P"), c(T, P))
ltext <- c(lTP, lS)
legend("topright", ltext, bty = "n", bg = "transparent")
title(quote("Deep Earth Water (DEW)"~"model"))
water(oldwat)
OBIGT()

## ----pyrrhotite_polymorphs, collapse = TRUE-----------------------------------
subcrt("pyrrhotite", T = c(25, 150, 350), property = "G")$out

## ----pyrite_limit-------------------------------------------------------------
subcrt("pyrite", T = seq(200, 1000, 200), P = 1)

## ----mineral_Ttr, collapse = TRUE---------------------------------------------
file <- system.file("extdata/OBIGT/inorganic_cr.csv", package = "CHNOSZ")
dat <- read.csv(file)
# Reverse rows so highest-T polymorph for each mineral is listed first
dat <- dat[nrow(dat):1, ]
# Remove low-T polymorphs
dat <- dat[!duplicated(dat$name), ]
# Remove minerals with no T limit for phase stability (+ve) or Cp equation (-ve)
dat <- dat[!is.na(dat$z.T), ]
# Keep minerals with phase stability limit
dat <- dat[dat$z.T > 0, ]
# Get names of minerals and put into original order
rev(dat$name)

## ----muscovite_limit, message = FALSE-----------------------------------------
add.OBIGT("SUPCRT92")
subcrt("muscovite", T = 850, P = 4500)
reset()

## ----KMQ_basis_species, message = FALSE---------------------------------------
basis(c("Al+3", "quartz", "K+", "H+", "H2O", "oxygen"))
species(c("kaolinite", "muscovite", "K-feldspar"))

## ----KMQ_m_K, message = FALSE-------------------------------------------------
# Define temperature, pressure, and molality of Cl- (==IS)
T <- 150
P <- 500
IS <- m_Cl <- 1
# Calculate equilibrium constant for Ab-Kfs reaction, corrected for ionic strength
logK_AK <- subcrt(c("albite", "K+", "K-feldspar", "Na+"), c(-1, -1, 1, 1),
  T = T, P = P, IS = IS)$out$logK
K_AK <- 10 ^ logK_AK
# Calculate molality of K+
(m_K <- m_Cl / (K_AK + 1))

## ----KMQ_diagram, eval = FALSE, echo = 2:10-----------------------------------
#  par(mfrow = c(1, 2))
#  basis("K+", log10(m_K))
#  a <- affinity(pH = c(2, 10), O2 = c(-55, -38), T = T, P = P, IS = IS)
#  diagram(a, srt = 90)
#  lTP <- as.expression(c(lT(T), lP(P)))
#  legend("topleft", legend = lTP, bty = "n", inset = c(-0.05, 0), cex = 0.9)
#  ltxt <- c(quote("Unit molality of Cl"^"-"), "Quartz saturation")
#  legend("topright", legend = ltxt, bty = "n", cex = 0.9)
#  title("Mineral data from Berman (1988)\nand Sverjensky et al. (1991) (OBIGT default)",
#    font.main = 1, cex.main = 0.9)
#  add.OBIGT("SUPCRT92")
#  a <- affinity(pH = c(2, 10), O2 = c(-55, -38), T = T, P = P, IS = IS)
#  diagram(a, srt = 90)
#  title("Mineral data from Helgeson et al. (1978)\n(as used in SUPCRT92)",
#    font.main = 1, cex.main = 0.9)
#  OBIGT()

## ----KMQ_diagram, message = FALSE, fig.width = 8, fig.height = 4, out.width = "100%", results = "hide", echo = FALSE----
par(mfrow = c(1, 2))
basis("K+", log10(m_K))
a <- affinity(pH = c(2, 10), O2 = c(-55, -38), T = T, P = P, IS = IS)
diagram(a, srt = 90)
lTP <- as.expression(c(lT(T), lP(P)))
legend("topleft", legend = lTP, bty = "n", inset = c(-0.05, 0), cex = 0.9)
ltxt <- c(quote("Unit molality of Cl"^"-"), "Quartz saturation")
legend("topright", legend = ltxt, bty = "n", cex = 0.9)
title("Mineral data from Berman (1988)\nand Sverjensky et al. (1991) (OBIGT default)",
  font.main = 1, cex.main = 0.9)
add.OBIGT("SUPCRT92")
a <- affinity(pH = c(2, 10), O2 = c(-55, -38), T = T, P = P, IS = IS)
diagram(a, srt = 90)
title("Mineral data from Helgeson et al. (1978)\n(as used in SUPCRT92)",
  font.main = 1, cex.main = 0.9)
OBIGT()

## ----KMQ_refs, message = FALSE------------------------------------------------
thermo.refs(species()$ispecies)

## ----KMQ_diagram, eval = FALSE, echo = 11:15----------------------------------
#  par(mfrow = c(1, 2))
#  basis("K+", log10(m_K))
#  a <- affinity(pH = c(2, 10), O2 = c(-55, -38), T = T, P = P, IS = IS)
#  diagram(a, srt = 90)
#  lTP <- as.expression(c(lT(T), lP(P)))
#  legend("topleft", legend = lTP, bty = "n", inset = c(-0.05, 0), cex = 0.9)
#  ltxt <- c(quote("Unit molality of Cl"^"-"), "Quartz saturation")
#  legend("topright", legend = ltxt, bty = "n", cex = 0.9)
#  title("Mineral data from Berman (1988)\nand Sverjensky et al. (1991) (OBIGT default)",
#    font.main = 1, cex.main = 0.9)
#  add.OBIGT("SUPCRT92")
#  a <- affinity(pH = c(2, 10), O2 = c(-55, -38), T = T, P = P, IS = IS)
#  diagram(a, srt = 90)
#  title("Mineral data from Helgeson et al. (1978)\n(as used in SUPCRT92)",
#    font.main = 1, cex.main = 0.9)
#  OBIGT()

