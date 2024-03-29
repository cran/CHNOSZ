## ----width80, include=FALSE---------------------------------------------------
options(width = 80)

## ----digits6, include=FALSE---------------------------------------------------
options(digits = 6)

## ----HTML, include=FALSE------------------------------------------------------
## Some frequently used HTML expressions
V0 <- "<i>V</i>&deg;"
Cp0 <- "<i>C<sub>P</sub></i>&deg;"
c1 <- "<i>c</i><sub>1</sub>"
c2 <- "<i>c</i><sub>2</sub>"
a1 <- "<i>a</i><sub>1</sub>"
a2 <- "<i>a</i><sub>2</sub>"
a3 <- "<i>a</i><sub>3</sub>"
a4 <- "<i>a</i><sub>4</sub>"
h4sio4 <- "H<sub>4</sub>SiO<sub>4</sub>"
sio2 <- "SiO<sub>2</sub>"
h2o <- "H<sub>2</sub>O"
ch4 <- "CH<sub>4</sub>"
wPrTr <- "&omega;<sub><i>P<sub>r</sub></i>,<i>T<sub>r</sub></i></sub>"

## ----setup, include=FALSE-----------------------------------------------------
library(knitr)
# Invalidate cache when the tufte version changes
opts_chunk$set(tidy = FALSE, cache.extra = packageVersion('tufte'))
options(htmltools.dir.version = FALSE)
# Adjust plot margins
knit_hooks$set(small.mar = function(before, options, envir) {
    if (before) par(mar = c(4.2, 4.2, .3, .3))  # smaller margin on top and right
})
# Use pngquant to optimize PNG images
knit_hooks$set(pngquant = hook_pngquant)
# pngquant isn't available on R-Forge ...
if (!nzchar(Sys.which("pngquant"))) {
  pngquant <- NULL 
} else {
  pngquant <- "--speed=1 --quality=0-50"
}

## Colorize messages 20171031
## Adapted from https://gist.github.com/yihui/2629886#file-knitr-color-msg-rnw
color_block = function(color) {
  function(x, options) sprintf('<pre style="color:%s">%s</pre>', color, x)
}
knit_hooks$set(warning = color_block('magenta'), error = color_block('red'), message = color_block('blue'))

# Set dpi 20231129
knitr::opts_chunk$set(
  dpi = if(nzchar(Sys.getenv("CHNOSZ_BUILD_LARGE_VIGNETTES"))) 72 else 50
)

## ----library_CHNOSZ-----------------------------------------------------------
library(CHNOSZ)
reset()

## ----Cpdat--------------------------------------------------------------------
file <- system.file("extdata/cpetc/HW97_Cp.csv", package = "CHNOSZ")
Cpdat <- read.csv(file)
# Use data for CH4
Cpdat <- Cpdat[Cpdat$species == "CH4", ]
Cpdat <- Cpdat[, -1]
Cpdat$P <- convert(Cpdat$P, "bar")

## ----EOSregress---------------------------------------------------------------
var <- c("invTTheta2", "TXBorn")
Cplm_high <- EOSregress(Cpdat, var)
Cplm_low <- EOSregress(Cpdat, var, T.max = 600)
Cplm_low$coefficients

## ----EOSplot, fig.margin = TRUE, fig.cap = "Heat capacity of aqueous methane.", fig.width=3.5, fig.height=3.5, cache=TRUE, results="hide", message=FALSE, echo=FALSE, out.width=672, out.height=336, pngquant=pngquant----
EOSplot(Cpdat, coefficients = round(Cplm_low$coefficients, 1))
EOSplot(Cpdat, coeficients = Cplm_high, add = TRUE, lty = 3)
PS01_data <- convert(EOScoeffs("CH4", "Cp"), "J")
EOSplot(Cpdat, coefficients = PS01_data, add = TRUE, lty = 2, col = "blue")

## ----EOSplot, eval=FALSE------------------------------------------------------
#  EOSplot(Cpdat, coefficients = round(Cplm_low$coefficients, 1))
#  EOSplot(Cpdat, coeficients = Cplm_high, add = TRUE, lty = 3)
#  PS01_data <- convert(EOScoeffs("CH4", "Cp"), "J")
#  EOSplot(Cpdat, coefficients = PS01_data, add = TRUE, lty = 2, col = "blue")

## ----Cpcoeffs, message=FALSE--------------------------------------------------
convert(EOScoeffs("CH4", "Cp"), "J")

## ----Vdat---------------------------------------------------------------------
file <- system.file("extdata/cpetc/HWM96_V.csv", package = "CHNOSZ")
Vdat <- read.csv(file)
# Use data for CH4 near 280 bar
Vdat <- Vdat[Vdat$species == "CH4", ]
Vdat <- Vdat[abs(Vdat$P - 28) < 0.1, ]
Vdat <- Vdat[, -1]
Vdat$P <- convert(Vdat$P, "bar")

## ----Vdat_non-----------------------------------------------------------------
QBorn <- EOSvar("QBorn", T = Vdat$T, P = Vdat$P)
Vomega <- convert(Cplm_low$coefficients[["TXBorn"]], "cm3bar")
V_sol <- Vomega * QBorn
V_non <- Vdat$V - V_sol
Vdat_non <- data.frame(T = Vdat$T, P = Vdat$P, V = V_non)

## ----Vdat_non_regress, message=FALSE------------------------------------------
var <- "invTTheta"
Vnonlm <- EOSregress(Vdat_non, var, T.max = 450)
Vcoeffs <- round(c(Vnonlm$coefficients, QBorn = Vomega), 1)
Vcoeffs_database <- convert(EOScoeffs("CH4", "V"), "J")

## ----Vplot, fig.margin=TRUE, results="hide", message=FALSE, echo=FALSE, fig.width=3.5, fig.height=7, fig.cap="Volume of aqueous methane.", out.width=672, out.height=672, pngquant=pngquant----
par(mfrow = c(2, 1))
# plot 1
EOSplot(Vdat, coefficients = Vcoeffs)
EOSplot(Vdat, coefficients = Vcoeffs_database, add = TRUE, lty = 2)
# plot 2
EOSplot(Vdat, coefficients = Vcoeffs_database, T.plot = 600, lty = 2)
EOSplot(Vdat, coefficients = Vcoeffs, add = TRUE)

## ----Vplot, eval=FALSE--------------------------------------------------------
#  par(mfrow = c(2, 1))
#  # plot 1
#  EOSplot(Vdat, coefficients = Vcoeffs)
#  EOSplot(Vdat, coefficients = Vcoeffs_database, add = TRUE, lty = 2)
#  # plot 2
#  EOSplot(Vdat, coefficients = Vcoeffs_database, T.plot = 600, lty = 2)
#  EOSplot(Vdat, coefficients = Vcoeffs, add = TRUE)

## ----Nadat--------------------------------------------------------------------
T <- convert(seq(0, 600, 50), "K")
P <- 1000
prop.PT <- subcrt("Na+", T = T, P = P, grid = "T", convert = FALSE)$out[[1]]
Nadat <- prop.PT[, c("T", "P", "Cp")]

## ----Nalm, fig.margin=TRUE, fig.width=3.5, fig.height=3.5, fig.cap="Heat capacity of Na<sup>+</sup> (inapplicable: constant &omega;).", out.width=672, out.height=336, pngquant=pngquant----
var <- c("invTTheta2", "TXBorn")
Nalm <- EOSregress(Nadat, var, T.max = 600)
EOSplot(Nadat, coefficients = Nalm$coefficients, fun.legend = NULL)
EOSplot(Nadat, add = TRUE, lty = 3)

## ----Navars1------------------------------------------------------------------
var1 <- c("invTTheta2", "Cp_s_var")
omega.guess <- coef(Nalm)[3]

## ----Nawhile, fig.margin=TRUE, fig.width=3.5, fig.height=3.5, fig.cap="Heat capacity of Na<sup>+</sup> (variable &omega;).", out.width=672, out.height=336, pngquant=pngquant----
diff.omega <- 999
while(abs(diff.omega) > 1) {
  Nalm1 <- EOSregress(Nadat, var1, omega.PrTr = tail(omega.guess, 1), Z = 1)
  omega.guess <- c(omega.guess, coef(Nalm1)[3])
  diff.omega <- tail(diff(omega.guess), 1)
}
EOSplot(Nadat, coefficients = signif(coef(Nalm1), 6),
  omega.PrTr = tail(omega.guess, 1), Z = 1)
convert(EOScoeffs("Na+", "Cp"), "J")

## ----NaVolume, fig.margin=TRUE, fig.width=3.5, fig.height=3.5, fig.cap="Volume of Na<sup>+</sup> (variable &omega;).", results="hide", message=FALSE, echo=FALSE, out.width=672, out.height=336, pngquant=pngquant----
T <- convert(seq(0, 600, 25), "K")
P <- 1000
prop.PT <- subcrt("Na+", T = T, P = P, grid = "T", convert = FALSE)$out[[1]]
NaVdat <- prop.PT[, c("T", "P", "V")]
var1 <- c("invTTheta", "V_s_var")
omega.guess <- 1400000
diff.omega <- 999
while(abs(diff.omega) > 1) {
  NaVlm1 <- EOSregress(NaVdat, var1,
    omega.PrTr = tail(convert(omega.guess, "joules"), 1), Z = 1)
  omega.guess <- c(omega.guess, coef(NaVlm1)[3])
  diff.omega <- tail(diff(omega.guess), 1)
}
EOSplot(NaVdat, coefficients = signif(coef(NaVlm1), 6),
  omega.PrTr = tail(convert(omega.guess, "joules"), 1), Z = 1,
  fun.legend = "bottomleft")
coefficients <- convert(EOScoeffs("Na+", "V", P = 1000), "J")
names(coefficients)[3] <- "V_s_var"
EOSplot(NaVdat, coefficients = coefficients, Z = 1, add = TRUE, lty = 2,
  omega.PrTr = convert(coefficients["V_s_var"], "joules"))

## ----NaVolume, eval=FALSE-----------------------------------------------------
#  T <- convert(seq(0, 600, 25), "K")
#  P <- 1000
#  prop.PT <- subcrt("Na+", T = T, P = P, grid = "T", convert = FALSE)$out[[1]]
#  NaVdat <- prop.PT[, c("T", "P", "V")]
#  var1 <- c("invTTheta", "V_s_var")
#  omega.guess <- 1400000
#  diff.omega <- 999
#  while(abs(diff.omega) > 1) {
#    NaVlm1 <- EOSregress(NaVdat, var1,
#      omega.PrTr = tail(convert(omega.guess, "joules"), 1), Z = 1)
#    omega.guess <- c(omega.guess, coef(NaVlm1)[3])
#    diff.omega <- tail(diff(omega.guess), 1)
#  }
#  EOSplot(NaVdat, coefficients = signif(coef(NaVlm1), 6),
#    omega.PrTr = tail(convert(omega.guess, "joules"), 1), Z = 1,
#    fun.legend = "bottomleft")
#  coefficients <- convert(EOScoeffs("Na+", "V", P = 1000), "J")
#  names(coefficients)[3] <- "V_s_var"
#  EOSplot(NaVdat, coefficients = coefficients, Z = 1, add = TRUE, lty = 2,
#    omega.PrTr = convert(coefficients["V_s_var"], "joules"))

## ----AS04, message=FALSE------------------------------------------------------
add.OBIGT("AS04")

## ----SiO2_2H2O, message=FALSE-------------------------------------------------
s_25C <- subcrt(c("SiO2", "H2O"), c(1, 2), T = 25)$out
s_near25 <- subcrt(c("SiO2", "H2O"), c(1, 2), T = seq(20, 30, length.out=50))$out
s_lowT <- subcrt(c("SiO2", "H2O"), c(1, 2), T = seq(0, 100, 10))$out
s_Psat <- subcrt(c("SiO2", "H2O"), c(1, 2))$out
s_P500 <- subcrt(c("SiO2", "H2O"), c(1, 2), T = seq(0, 1000, 100), P = 500)$out
s_P1000 <- subcrt(c("SiO2", "H2O"), c(1, 2), T = seq(0, 1000, 100), P = 1000)$out

## ----new_H4SiO4---------------------------------------------------------------
mod.OBIGT("calc-H4SiO4", formula = "H4SiO4", ref1 = "this_vignette",
          date = as.character(Sys.Date()), G = s_25C$G, H = s_25C$H, S = s_25C$S,
          Cp = s_25C$Cp, V = s_25C$V, z = 0)

## ----substuff-----------------------------------------------------------------
substuff <- rbind(s_near25, s_lowT, s_Psat, s_P500, s_P1000)
substuff$T <- convert(substuff$T, "K")

## ----Cp_H4SiO4, results="hide"------------------------------------------------
Cpdat <- substuff[, c("T", "P", "Cp")]
var <- c("invTTheta2", "TXBorn")
Cplm <- EOSregress(Cpdat, var) 
Cpcoeffs <- Cplm$coefficients
mod.OBIGT("calc-H4SiO4", c1 = Cpcoeffs[1],
  c2 = Cpcoeffs[2]/10000, omega = Cpcoeffs[3]/100000)

## ----V_H4SiO4_nonsolvation----------------------------------------------------
Vdat <- substuff[, c("T", "P", "V")]
QBorn <- EOSvar("QBorn", T = Vdat$T, P = Vdat$P)
Vomega <- convert(Cplm$coefficients[["TXBorn"]], "cm3bar")
V_sol <- Vomega * QBorn
V_non <- Vdat$V - V_sol
Vdat$V <- V_non

## ----V_H4SiO4, results="hide"-------------------------------------------------
var <- c("invPPsi", "invTTheta", "invPPsiTTheta")
Vlm <- EOSregress(Vdat, var)
Vcoeffs <- convert(Vlm$coefficients, "joules")
mod.OBIGT("calc-H4SiO4", a1 = Vcoeffs[1]*10, a2 = Vcoeffs[2]/100,
  a3 = Vcoeffs[3], a4 = Vcoeffs[4]/10000)

## ----width180, include=FALSE------------------------------------------------------------------------------------------------------------------------------------------------------
options(width = 180)

## ----info_H4SiO4, message=FALSE---------------------------------------------------------------------------------------------------------------------------------------------------
info(info(c("calc-H4SiO4", "H4SiO4")))

## ----width80, include=FALSE---------------------------------------------------
options(width = 80)

## ----subcrt_H4SiO4, fig.margin=TRUE, fig.width=4, fig.height=4, small.mar=TRUE, echo=FALSE, results="hide", message=FALSE, out.width="100%", cache=TRUE, fig.cap="Comparison of H<sub>4</sub>SiO<sub>4</sub> pseudospecies.", pngquant=pngquant----
s1 <- subcrt(c("calc-H4SiO4", "SiO2", "H2O"), c(-1, 1, 2))
plot(s1$out$T, s1$out$G, type = "l", ylim = c(-500, 2000),
  xlab = axis.label("T"), ylab = axis.label("DG0"))
s2 <- subcrt(c("H4SiO4", "SiO2", "H2O"), c(-1, 1, 2))
lines(s2$out$T, s2$out$G, lty = 2)
abline(h = 0, lty = 3)
legend("topright", legend = c("calc-H4SiO4 (this vignette)",
  "H4SiO4 (Stef\u00e1nsson, 2001)"), lty = c(1, 2), bty = "n")
text(225, 1000, describe.reaction(s1$reaction))

## ----subcrt_H4SiO4, eval=FALSE------------------------------------------------
#  s1 <- subcrt(c("calc-H4SiO4", "SiO2", "H2O"), c(-1, 1, 2))
#  plot(s1$out$T, s1$out$G, type = "l", ylim = c(-500, 2000),
#    xlab = axis.label("T"), ylab = axis.label("DG0"))
#  s2 <- subcrt(c("H4SiO4", "SiO2", "H2O"), c(-1, 1, 2))
#  lines(s2$out$T, s2$out$G, lty = 2)
#  abline(h = 0, lty = 3)
#  legend("topright", legend = c("calc-H4SiO4 (this vignette)",
#    "H4SiO4 (Stef\u00e1nsson, 2001)"), lty = c(1, 2), bty = "n")
#  text(225, 1000, describe.reaction(s1$reaction))

## ----activity_diagram, fig.margin=TRUE, fig.width=4, fig.height=4, small.mar=TRUE, echo=TRUE, results="hide", message=FALSE, out.width="100%", cache=TRUE, fig.cap="Activity diagram for K<sub>2</sub>O-Al<sub>2</sub>O<sub>3</sub>-SiO<sub>2</sub>-H<sub>2</sub>O.", pngquant=pngquant----
basis(c("Al+3", "H4SiO4", "K+", "H2O", "H+", "O2"))
species(c("gibbsite", "muscovite", "kaolinite", "pyrophyllite", "K-feldspar"))
a <- affinity(H4SiO4 = c(-8, 0, 300), `K+` = c(-1, 8, 300))
diagram(a, ylab = ratlab("K+"), fill = "terrain", yline = 1.7)
legend("bottomleft", describe.property(c("T", "P"), c(25, 1)), bty = "n")

