## ----setup, include = FALSE---------------------------------------------------
library(CHNOSZ)
options(width = 80)
## Use pngquant to optimize PNG images
library(knitr)
knit_hooks$set(pngquant = hook_pngquant)
pngquant <- "--speed=1 --quality=0-25"
if (!nzchar(Sys.which("pngquant"))) pngquant <- NULL

# Set dpi 20231129
knitr::opts_chunk$set(
  dpi = if(nzchar(Sys.getenv("CHNOSZ_BUILD_LARGE_VIGNETTES"))) 100 else 72
)

## ----HTML, include = FALSE----------------------------------------------------
NOTE <- '<span style="background-color: yellow;">NOTE</span>'
# OBIGT columns
model <- '<span style="color: blue;">model</span>'
name <- '<tt style="color: blue;">name</tt>'
abbrv <- '<tt style="color: blue;">abbrv</tt>'
formula <- '<tt style="color: blue;">formula</tt>'
state <- '<tt style="color: blue;">state</tt>'
ref1 <- '<tt style="color: blue;">ref1</tt>'
ref2 <- '<tt style="color: blue;">ref2</tt>'
date <- '<tt style="color: blue;">date</tt>'
E_units <- '<tt style="color: blue;">E_units</tt>'
b <- '<tt style="color: blue;">b</tt>'
c <- '<tt style="color: blue;">c</tt>'
e <- '<tt style="color: blue;">e</tt>'
f <- '<tt style="color: blue;">f</tt>'
lambda <- '<tt style="color: blue;">lambda</tt>'
c1 <- '<tt style="color: blue;">c1</tt>'
c2 <- '<tt style="color: blue;">c2</tt>'
omega <- '<tt style="color: blue;">omega</tt>'
G_ <- '<tt style="color: blue;">G</tt>'
H_ <- '<tt style="color: blue;">H</tt>'
S_ <- '<tt style="color: blue;">S</tt>'
Cp_ <- '<tt style="color: blue;">Cp</tt>'
V_ <- '<tt style="color: blue;">V</tt>'
T_ <- '<tt style="color: blue;">T</tt>'
# CHNOSZ functions
reset_ <- '<code style="color: red;">reset()</code>'
OBIGT_ <- '<code style="color: red;">OBIGT()</code>'
add.OBIGT_ <- '<code style="color: red;">add.OBIGT()</code>'
mod.OBIGT_ <- '<code style="color: red;">mod.OBIGT()</code>'
logB.to.OBIGT_ <- '<code style="color: red;">logB.to.OBIGT()</code>'
basis_ <- '<code style="color: red;">basis()</code>'
species_ <- '<code style="color: red;">species()</code>'
E.units_ <- '<code style="color: red;">E.units()</code>'
info_ <- '<code style="color: green;">info()</code>'
subcrt_ <- '<code style="color: green;">subcrt()</code>'
affinity_ <- '<code style="color: green;">affinity()</code>'
thermo.refs_ <- '<code style="color: green;">thermo.refs()</code>'
thermo_ <- '<code style="color: green;">thermo()</code>'
check.GHS_ <- '<code style="color: green;">check.GHS()</code>'
check.EOS_ <- '<code style="color: green;">check.EOS()</code>'
# Math stuff
logK <- "log&thinsp;<i>K</i>"
logB <- "log&thinsp;&beta;"
F_ <- "F<sup>-</sup>"
Hplus <- "H<sup>+</sup>"
HWO4_ <- "HWO<sub>4</sub><sup>-</sup>"
H2WO4 <- "H<sub>2</sub>WO<sub>4</sub>"
H3WO4F2_ <- "H<sub>3</sub>WO<sub>4</sub>F<sub>2</sub><sup>-</sup>"
# Thermodynamic properties
Cp_0 <- "<i>C<sub>p</sub></i>&deg;"
DG_0 <- "&Delta;<i>G</i>&deg;"

## ----system.file--------------------------------------------------------------
system.file("extdata/OBIGT", package = "CHNOSZ")

## ----dir.system.file----------------------------------------------------------
dir(system.file("extdata/OBIGT", package = "CHNOSZ"))

## ----reset_-------------------------------------------------------------------
reset()

## ----thermo.OBIGT-------------------------------------------------------------
head(thermo()$OBIGT)

## ----colnames.OBIGT, echo = FALSE---------------------------------------------
paste(1:22, colnames(thermo()$OBIGT))

## ----icr, message = FALSE-----------------------------------------------------
icr <- info(c("orpiment,amorphous", "arsenic,alpha", "tin"))
thermo()$OBIGT[icr, ]

## ----orpiment-----------------------------------------------------------------
subcrt("orpiment,amorphous", T = c(25, 50, 75))$out[[1]]

## ----arsenic------------------------------------------------------------------
subcrt("arsenic,alpha", T = c(25, 50, 75))$out[[1]]

## ----tin----------------------------------------------------------------------
subcrt("tin", T = c(25, 50, 75))$out[[1]]

## ----info_.tin----------------------------------------------------------------
info(info("tin"))

## ----Berman-------------------------------------------------------------------
info(info(c("quartz", "pyrite")))

## ----add.OBIGT_quartz---------------------------------------------------------
add.OBIGT("SUPCRT92", "quartz")
info(info("quartz"))

## ----add.OBIGT_SUPCRT92-------------------------------------------------------
iSUPCRT92 <- add.OBIGT("SUPCRT92")
unique(suppressMessages(info(iSUPCRT92))$name)

## ----BZA10--------------------------------------------------------------------
file <- system.file("extdata/adds/BZA10.csv", package = "CHNOSZ")
read.csv(file, as.is = TRUE)

## ----BZA10_Cd-----------------------------------------------------------------
iCd <- add.OBIGT(file)
subcrt(c("CdCl+", "Cl-", "CdCl2"), c(-1, -1, 1), T = 25, P = c(1, 2000))

## ----SSH97_subcrt-------------------------------------------------------------
reset()
thermo.refs(iCd)[, 1:3]
subcrt(c("CdCl+", "Cl-", "CdCl2"), c(-1, -1, 1), T = 25, P = c(1, 2000))

## ----mod.OBIGT__CoCl4_ghs-----------------------------------------------------
mod.OBIGT("CoCl4-2", formula = "CoCl4-2", state = "aq", ref1 = "LBT+11", E_units = "cal",
  date = as.character(Sys.Date()), G = -134150, H = -171558, S = 19.55, Cp = 72.09, V = 27.74)

## ----mod.OBIGT__CoCl4_eos-----------------------------------------------------
mod.OBIGT("CoCl4-2", a1 = 6.5467, a2 = 8.2069, a3 = 2.0130, a4 = -3.1183,
  c1 = 76.3357, c2 = 11.6389, omega = 2.9159, z = -2)

## ----CoCl4_reaction, message = FALSE, echo = 1:3------------------------------
T <- c(25, seq(50, 350, 50))
sres <- subcrt(c("Co+2", "Cl-", "CoCl4-2"), c(-1, -4, 1), T = T)
round(sres$out$logK, 2)
stopifnot(identical(round(sres$out$logK, 2), c(-3.2, -2.96, -2.02, -0.74, 0.77, 2.5, 4.57, 7.29)))

## ----info__CoCl4, results = "hide"--------------------------------------------
inew <- info("CoCl4-2")
info(inew)

## ----mod.OBIGT__magnesiochromite_ghs------------------------------------------
H <- -1762000
S <- 119.6
V <- 43.56
mod.OBIGT("magnesiochromite", formula = "MgCr2O4", state = "cr", ref1 = "KOSG00",
          date = as.character(Sys.Date()), E_units = "J", H = H, S = S, V = V)

## ----mod.OBIGT__magnesiochromite_eos------------------------------------------
a <- 221.4
b <- -0.00102030
c <- -1757210
d <- -1247.9
mod.OBIGT("magnesiochromite", E_units = "J", a = a, b = b, c = c, d = d,
          e = 0, f = 0, lambda = 0, T = 1500)

## ----subcrt__magnesiochromite-------------------------------------------------
T.units("K")
Tref <- c(250, 300, 340)
(sres <- subcrt("magnesiochromite", property = "Cp", T = Tref, P = 1))

## ----magnesiochromite_check_Cp------------------------------------------------
Cpref <- c(114.3, 129.8, 138.4)
stopifnot(max(abs(sres$out[[1]]$Cp - Cpref)) < 0.3)

## ----restore_units_magnesiochromite-------------------------------------------
T.units("C")

## ----Psat---------------------------------------------------------------------
P <- "Psat"

## ----HWO4_--------------------------------------------------------------------
T <- c(250, 300, 350)
logB <- c(5.58, 6.51, 7.99)
species <- c("WO4-2", "H+", "HWO4-")
coeff <- c(-1, -1, 1)
logB.to.OBIGT(logB, species, coeff, T, P)

## ----H3WO4F2------------------------------------------------------------------
T <- seq(100, 250, 25)
logB <- c(17.00, 17.11, 17.46, 17.75, 18.17, 18.71, 19.23)
# Species and coefficients in the formation reaction
species <- c("H+", "WO4-2", "F-", "H3WO4F2-")
coeff <- c(-3, -1, -2, 1)
logB.to.OBIGT(logB, species, coeff, T, P)

## ----H2WO4--------------------------------------------------------------------
logB <- c(7.12, 7.82, 7.07, 7.76, 7.59, 7.98, 8.28)
species <- c("H+", "WO4-2", "H2WO4")
coeff <- c(-2, -1, 1)
logB.to.OBIGT(logB, species, coeff, T, P, tolerance = 0.3)

## ----diagram1, message = FALSE, results = "hide", fig.width = 6, fig.height = 5, out.width = "75%", fig.align = "center", pngquant = pngquant----
basis(c("H+", "WO4-2", "F-", "H2O", "O2"))
basis("F-", log10(0.1))
iaq <- retrieve("W", c("O", "H", "F"), "aq")
species(iaq)
a <- affinity(pH = c(2, 7), T = 300, IS = 0.9)
e <- equilibrate(a)
col <- c(1, 4, 5, 2)
diagram(e, alpha = TRUE, col = col, lty = 1, lwd = 2, ylab = "Fraction total W")

## ----a_F----------------------------------------------------------------------
T <- 300
pH <- seq(2, 7, 0.1)
logK_HF <- subcrt(c("H+", "F-", "HF"), c(-1, -1, 1), T = T, IS = 0.9)$out$logK
F_tot <- 0.1
a_F <- F_tot / (1 + 10^(logK_HF - pH))

## ----diagram2, message = FALSE, results = "hide", results = "hide", fig.width = 6, fig.height = 5, out.width = "75%", fig.align = "center", pngquant = pngquant----
basis(c("H+", "WO4-2", "F-", "H2O", "O2"))
iaq <- retrieve("W", c("O", "H", "F"), "aq")
species(iaq)
a <- affinity(pH = pH, "F-" = log10(a_F), T = T, IS = 0.9)
e <- equilibrate(a)
diagram(e, alpha = TRUE, col = col, lty = 1, lwd = 2, ylab = "Fraction total W")

