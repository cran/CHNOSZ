## ----options, include=FALSE---------------------------------------------------
options(width = 80)
options(digits = 6)

## ----HTML, include=FALSE------------------------------------------------------
## some frequently used HTML expressions
logfO2 <- "log<i>f</i><sub>O<sub>2</sub></sub>"
# use lowercase here because these tend to be variable names in the examples
zc <- "<i>Z</i><sub>C</sub>"
o2 <- "O<sub>2</sub>"
h2o <- "H<sub>2</sub>O"
sio2 <- "SiO<sub>2</sub>"
ch4 <- "CH<sub>4</sub>"

## ----setup, include=FALSE-----------------------------------------------------
library(knitr)

## from "Tufte Handout" example dated 2016-12-27
# invalidate cache when the tufte version changes
opts_chunk$set(tidy = FALSE, cache.extra = packageVersion('tufte'))
options(htmltools.dir.version = FALSE)

## adjust plot margins
## first one from https://yihui.name/knitr/hooks/
knit_hooks$set(small.mar = function(before, options, envir) {
    if (before) par(mar = c(4.2, 4.2, .1, .1))  # smaller margin on top and right
})
knit_hooks$set(tiny.mar = function(before, options, envir) {
    if (before) par(mar = c(.1, .1, .1, .1))  # tiny margin all around
})
knit_hooks$set(smallish.mar = function(before, options, envir) {
    if (before) par(mar = c(4.2, 4.2, 0.9, 0.9))  # smallish margins on top and right
})

## use pngquant to optimize PNG images
knit_hooks$set(pngquant = hook_pngquant)
pngquant <- "--speed=1 --quality=0-25"
# pngquant isn't available on R-Forge ...
if (!nzchar(Sys.which("pngquant"))) pngquant <- NULL 

## use a low resolution to save space in the package
# change this to 72 to make higher-resolution images for the CHNOSZ web page
dpi <- 50

## http://stackoverflow.com/questions/23852753/knitr-with-gridsvg
## Set up a chunk hook for manually saved plots.
knit_hooks$set(custom.plot = hook_plot_custom)

## hook to change <img /> to <embed /> -- required for interactive SVG
hook_plot <- knit_hooks$get("plot")
knit_hooks$set(plot = function(x, options) {
  x <- hook_plot(x, options)
  if (!is.null(options$embed.tag) && options$embed.tag) x <- gsub("<img ", "<embed ", x)
  x
})

## http://stackoverflow.com/questions/30530008/hook-to-time-knitr-chunks
now = Sys.time()
knit_hooks$set(timeit = function(before) {
    if (before) { now <<- Sys.time() }
    else {
        paste("%", sprintf("Chunk rendering time: %s seconds.\n", round(Sys.time() - now, digits = 3))) 
    }
})
timeit <- NULL

## colorize messages 20171031
## adapted from https://gist.github.com/yihui/2629886#file-knitr-color-msg-rnw
color_block = function(color) {
  function(x, options) sprintf('<pre style="color:%s">%s</pre>', color, x)
}
knit_hooks$set(warning = color_block('magenta'), error = color_block('red'), message = color_block('blue'))

## ----install_CHNOSZ, eval=FALSE-----------------------------------------------
#  install.packages("CHNOSZ")

## ----library_CHNOSZ-----------------------------------------------------------
library(CHNOSZ)

## ----reset--------------------------------------------------------------------
reset()

## ----pseudocode, eval=FALSE---------------------------------------------------
#  basis(...)
#  species(...)
#  a <- affinity(...)
#  e <- equilibrate(a)  ## optional
#  diagram(e)           ## or diagram(a)
#  reset()         ## clear settings for next calculation

## ----info_CH4-----------------------------------------------------------------
info("CH4")

## ----info_CH4_gas-------------------------------------------------------------
info("CH4", "gas")

## ----info_names_gas-----------------------------------------------------------
info("methane")
info("oxygen")
info("carbon dioxide")

## ----info_S_S2----------------------------------------------------------------
info("S")
info("S2")

## ----iCH4, message=FALSE------------------------------------------------------
iCH4 <- info("CH4")
info(iCH4)

## ----info_info_water----------------------------------------------------------
info(info("water"))

## ----width180, include=FALSE------------------------------------------------------------------------------------------------------------------------------------------------------
options(width = 180)

## ----info_acid--------------------------------------------------------------------------------------------------------------------------------------------------------------------
info("acid")

## ----width80, include=FALSE---------------------------------------------------
options(width = 80)

## ----info_ribose--------------------------------------------------------------
info(" ribose")

## ----info_CH4_formula, message=FALSE------------------------------------------
info(iCH4)$formula

## ----makeup_iCH4--------------------------------------------------------------
makeup(iCH4)
as.chemical.formula(makeup(iCH4))

## ----ZC_iCH4, message=FALSE---------------------------------------------------
ZC(iCH4)
ZC(info(iCH4)$formula)
ZC(makeup(iCH4))

## ----subcrt_water-------------------------------------------------------------
subcrt("water")

## ----subcrt_water_grid--------------------------------------------------------
subcrt("water", T = c(400, 500, 600), P = c(200, 400, 600), grid = "P")$out$water

## ----subcrt_water_plot, fig.margin=TRUE, fig.width=4, fig.height=4, small.mar=TRUE, dpi=dpi, out.width="100%", echo=FALSE, message=FALSE, fig.cap="Isothermal contours of density (g cm<sup>-3</sup>) and pressure (bar) of water.", cache=TRUE, pngquant=pngquant, timeit=timeit----
substuff <- subcrt("water", T=seq(0,1000,100), P=c(NA, seq(1,500,1)), grid="T")
water <- substuff$out$water
plot(water$P, water$rho, type = "l")

## ----subcrt_water_plot, eval=FALSE--------------------------------------------
#  substuff <- subcrt("water", T=seq(0,1000,100), P=c(NA, seq(1,500,1)), grid="T")
#  water <- substuff$out$water
#  plot(water$P, water$rho, type = "l")

## ----units_CH4, message=FALSE-------------------------------------------------
T.units("K")
P.units("MPa")
E.units("J")
subcrt("CH4", T = 298.15, P = 0.1)$out$CH4$G

## ----convert_G, message=FALSE-------------------------------------------------
convert(info(info("CH4"))$G, "J")

## ----reset--------------------------------------------------------------------
reset()

## ----subcrt_CO2---------------------------------------------------------------
subcrt(c("CO2", "CO2"), c("gas", "aq"), c(-1, 1), T = seq(0, 250, 50))

## ----CO2_logK, echo=FALSE, message=FALSE--------------------------------------
T <- seq(0, 350, 10)
CO2 <- subcrt(c("CO2", "CO2"), c("gas", "aq"), c(-1, 1), T = T)$out$logK
CO <- subcrt(c("CO", "CO"), c("gas", "aq"), c(-1, 1), T = T)$out$logK
CH4 <- subcrt(c("CH4", "CH4"), c("gas", "aq"), c(-1, 1), T = T)$out$logK
logK <- data.frame(T, CO2, CO, CH4)

## ----CO2_plot, fig.margin=TRUE, fig.width=4, fig.height=4, small.mar=TRUE, dpi=dpi, out.width="100%", echo=FALSE, message=FALSE, fig.cap="Calculated equilibrium constants for dissolution of CO<sub>2</sub>, CO, and CH<sub>4</sub>.", cache=TRUE, pngquant=pngquant, timeit=timeit----
matplot(logK[, 1], logK[, -1], type = "l", col = 1, lty = 1,
        xlab = axis.label("T"), ylab = axis.label("logK"))
text(80, -1.7, expr.species("CO2"))
text(240, -2.37, expr.species("CO"))
text(300, -2.57, expr.species("CH4"))

## ----CO2_logK, eval=FALSE-----------------------------------------------------
#  T <- seq(0, 350, 10)
#  CO2 <- subcrt(c("CO2", "CO2"), c("gas", "aq"), c(-1, 1), T = T)$out$logK
#  CO <- subcrt(c("CO", "CO"), c("gas", "aq"), c(-1, 1), T = T)$out$logK
#  CH4 <- subcrt(c("CH4", "CH4"), c("gas", "aq"), c(-1, 1), T = T)$out$logK
#  logK <- data.frame(T, CO2, CO, CH4)

## ----CO2_plot, eval=FALSE-----------------------------------------------------
#  matplot(logK[, 1], logK[, -1], type = "l", col = 1, lty = 1,
#          xlab = axis.label("T"), ylab = axis.label("logK"))
#  text(80, -1.7, expr.species("CO2"))
#  text(240, -2.37, expr.species("CO"))
#  text(300, -2.57, expr.species("CH4"))

## ----subcrt_unbalanced, results="hide"----------------------------------------
subcrt(c("CO2", "CH4"), c(-1, 1))

## ----basis_singular, error=TRUE-----------------------------------------------
basis(c("CO2", "H2", "H2CO2"))

## ----basis_CHO----------------------------------------------------------------
basis(c("CO2", "H2", "H2O"))

## ----basis_CHOZ---------------------------------------------------------------
basis(c("CO2", "H2", "H2O", "H+"))

## ----subcrt_acetoclastic, message=FALSE---------------------------------------
subcrt(c("acetate", "CH4"), c(-1, 1))$reaction

## ----subcrt_methanogenesis, message=FALSE-------------------------------------
acetate_oxidation <- subcrt("acetate", -1)
hydrogenotrophic <- subcrt("CH4", 1)
acetoclastic <- subcrt(c("acetate", "CH4"), c(-1, 1))

## ----describe_reaction_plot, fig.margin=TRUE, fig.width=3.5, fig.height=1.8, tiny.mar=TRUE, dpi=dpi, out.width="100%", pngquant=pngquant, timeit=timeit----
plot(0, 0, type = "n", axes = FALSE, ann=FALSE, xlim=c(0, 5), ylim=c(5.2, -0.2))
text(0, 0, "acetoclastic methanogenesis", adj = 0)
text(5, 1, describe.reaction(acetoclastic$reaction), adj = 1)
text(0, 2, "acetate oxidation", adj = 0)
text(5, 3, describe.reaction(acetate_oxidation$reaction), adj = 1)
text(0, 4, "hydrogenotrophic methanogenesis", adj = 0)
text(5, 5, describe.reaction(hydrogenotrophic$reaction), adj = 1)

## ----basis_mayumi, message=FALSE, results="hide"------------------------------
E.units("J")
basis(c("CO2", "H2", "H2O", "H+"))
basis(c("CO2", "H2"), "gas")
basis(c("H2", "pH"), c(-3.92, 7.3))

## ----affinity_acetoclastic, message=FALSE-------------------------------------
subcrt(c("acetate", "CH4"), c(-1, 1),
       c("aq", "gas"), logact = c(-3.4, -0.18), T = 55, P = 50)$out

## ----affinity_hydrogenotrophic, message=FALSE---------------------------------
subcrt("CH4", 1, "gas", logact = -0.18, T = 55, P = 50)$out

## ----rxnfun, message=FALSE----------------------------------------------------
rxnfun <- function(coeffs) {
  subcrt(c("acetate", "CH4"), coeffs,
         c("aq", "gas"), logact = c(-3.4, -0.18), T = 55, P = 50)$out
}

## ----methanogenesis_plot, fig.margin=TRUE, fig.width=4.1, fig.height=4.1, small.mar=TRUE, dpi=dpi, out.width="100%", echo=FALSE, message=FALSE, fig.cap="Gibbs energies of acetate oxidation and methanogenesis (after Mayumi et al., 2013).", cache=TRUE, pngquant=pngquant, timeit=timeit----
Adat <- lapply(c(-3, 3), function(logfCO2) {
  basis("CO2", logfCO2)
  data.frame(logfCO2,
    rxnfun(c(0, 0))$A,
    rxnfun(c(-1, 0))$A,
    rxnfun(c(-1, 1))$A,
    rxnfun(c(0, 1))$A
  )
})
Adat <- do.call(rbind, Adat)
matplot(Adat[, 1], -Adat[, -1]/1000, type = "l", lty = 1, lwd = 2,
  xlab = axis.label("CO2"), ylab = axis.label("DG", prefix = "k"))
legend("topleft", c("acetate oxidation", "acetoclastic methanogenesis",
  "hydrogenotrophic methanogenesis"), lty = 1, col = 2:4)

## ----methanogenesis_plot, eval=FALSE------------------------------------------
#  Adat <- lapply(c(-3, 3), function(logfCO2) {
#    basis("CO2", logfCO2)
#    data.frame(logfCO2,
#      rxnfun(c(0, 0))$A,
#      rxnfun(c(-1, 0))$A,
#      rxnfun(c(-1, 1))$A,
#      rxnfun(c(0, 1))$A
#    )
#  })
#  Adat <- do.call(rbind, Adat)
#  matplot(Adat[, 1], -Adat[, -1]/1000, type = "l", lty = 1, lwd = 2,
#    xlab = axis.label("CO2"), ylab = axis.label("DG", prefix = "k"))
#  legend("topleft", c("acetate oxidation", "acetoclastic methanogenesis",
#    "hydrogenotrophic methanogenesis"), lty = 1, col = 2:4)

## ----reset, message=FALSE-----------------------------------------------------
reset()

## ----basis_CHNOSZ, results="hide"---------------------------------------------
basis("CHNOSe")

## ----species_sulfur-----------------------------------------------------------
species(c("H2S", "HS-", "HSO4-", "SO4-2"))

## ----affinity-----------------------------------------------------------------
unlist(affinity()$values)

## ----EhpH_plot, fig.margin=TRUE, fig.width=4, fig.height=4, dpi=dpi, out.width="100%", echo = FALSE, message=FALSE, cache=TRUE, fig.cap="Aqueous sulfur species at 25 °C.", pngquant=pngquant, timeit=timeit----
a <- affinity(pH = c(0, 12), Eh = c(-0.5, 1))
diagram(a, limit.water = TRUE)

## ----EhpH_plot, echo=TRUE, eval=FALSE-----------------------------------------
#  a <- affinity(pH = c(0, 12), Eh = c(-0.5, 1))
#  diagram(a, limit.water = TRUE)

## ----EhpH_plot_color, fig.margin=TRUE, fig.width=4, fig.height=4, smallish.mar=TRUE, dpi=dpi, out.width="100%", echo=FALSE, message=FALSE, cache=TRUE, fig.cap="The same plot, with different colors and labels.", pngquant=pngquant, timeit=timeit----
diagram(a, fill = "terrain", lwd = 2, lty = 3,
        names = c("hydrogen sulfide", "bisulfide", "bisulfate", "sulfate"),
        las = 0)
water.lines(a, col = 6, lwd = 2)

## ----EhpH_plot_color, echo=TRUE, eval=FALSE-----------------------------------
#  diagram(a, fill = "terrain", lwd = 2, lty = 3,
#          names = c("hydrogen sulfide", "bisulfide", "bisulfate", "sulfate"),
#          las = 0)
#  water.lines(a, col = 6, lwd = 2)

## ----retrieve-----------------------------------------------------------------
retrieve("Mn", c("O", "H"), "aq")
retrieve("Mn", c("O", "H"), "cr")

## ----retrieve_diagram, fig.margin=TRUE, fig.width=5, fig.height=5, dpi=dpi, out.width="100%", message=FALSE, results = "hide", cache=TRUE, fig.cap="Eh-pH diagram for the Mn-O-H system.", pngquant=pngquant, timeit=timeit----
# Set decimal logarithm of activity of aqueous species,
# temperature and plot resolution
logact <- -4
T <- 100
res <- 400
# Start with the aqueous species
basis(c("Mn+2", "H2O", "H+", "e-"))
iaq <- retrieve("Mn", c("O", "H"), "aq")
species(iaq, logact)
aaq <- affinity(pH = c(4, 16, res), Eh = c(-1.5, 1.5, res), T = T)
# Show names for only the metastable species here
names <- names(iaq)
names[!names(iaq) %in% c("MnOH+", "MnO", "HMnO2-")] <- ""
diagram(aaq, lty = 2, col = "#4169E188", names = names, col.names = 4)
# Overlay mineral stability fields
icr <- retrieve("Mn", c("O", "H"), "cr")
species(icr, add = TRUE)
# Supply the previous result from affinity() to use
# argument recall (for plotted variables and T)
acr <- affinity(aaq)
diagram(acr, add = TRUE, bold = acr$species$state=="cr", limit.water = FALSE)
# Add legend
legend <- c(
  bquote(log * italic(a)["Mn(aq)"] == .(logact)),
  bquote(italic(T) == .(T) ~ degree*C)
)
legend("topright", legend = as.expression(legend), bty = "n")

## ----info_CuCl, results="hide"------------------------------------------------
info(" CuCl")

## ----copper_setup, echo=TRUE, results="hide"----------------------------------
basis(c("Cu", "H2S", "Cl-", "H2O", "H+", "e-"))
basis("H2S", -6)
basis("Cl-", -0.7)
species(c("CuCl", "CuCl2-", "CuCl3-2", "CuCl+", "CuCl2", "CuCl3-", "CuCl4-2"))
species(c("chalcocite", "tenorite", "cuprite", "copper"), add = TRUE)

## ----info_chalcocite, message=FALSE-------------------------------------------
info(info("chalcocite", c("cr", "cr2", "cr3")))$T

## ----copper_mosaic, fig.margin=TRUE, fig.width=4, fig.height=4, dpi=dpi, out.width="100%", message=FALSE, cache=TRUE, fig.cap="Copper minerals and aqueous complexes with chloride, 200 °C.", pngquant=pngquant, timeit=timeit----
T <- 200
res <- 200
bases <- c("H2S", "HS-", "HSO4-", "SO4-2")
m1 <- mosaic(bases, pH = c(0, 12, res), Eh=c(-1.2, 0.75, res), T=T)
diagram(m1$A.species, lwd = 2)
diagram(m1$A.bases, add = TRUE, col = 4, col.names = 4, lty = 3,
        italic = TRUE)
water.lines(m1$A.species, col = "blue1")

## ----rainbow_data-------------------------------------------------------------
file <- system.file("extdata/cpetc/SC10_Rainbow.csv", package = "CHNOSZ")
rb <- read.csv(file, check.names = FALSE)

## ----rainbow_species, results="hide"------------------------------------------
basis(c("CO2", "H2", "NH4+", "H2O", "H2S", "H+"))
species("CH4", -3)
species(c("adenine", "cytosine", "aspartic acid", "deoxyribose",
          "CH4", "leucine", "tryptophan", "n-nonanoic acid"), -6)

## ----rainbow_affinity, message=FALSE------------------------------------------
a <- affinity(T = rb$T, CO2 = rb$CO2, H2 = rb$H2,
              `NH4+` = rb$`NH4+`, H2S = rb$H2S, pH = rb$pH)
T <- convert(a$vals[[1]], "K")
a$values <- lapply(a$values, convert, "G", T)
a$values <- lapply(a$values, `*`, -0.001)

## ----rainbow_diagram, fig.margin=TRUE, fig.width=4, fig.height=4, dpi=dpi, out.width="100%", echo=FALSE, message=FALSE, cache=TRUE, fig.cap="Affinities of organic synthesis in a hydrothermal system, after Shock and Canovas (2010).", pngquant=pngquant, timeit=timeit----
diagram(a, balance = 1, ylim = c(-100, 100), ylab = axis.label("A", prefix="k"),
        col = rainbow(8), lwd = 2, bg = "slategray3")
abline(h = 0, lty = 2, lwd = 2)

## ----rainbow_diagram, eval=FALSE----------------------------------------------
#  diagram(a, balance = 1, ylim = c(-100, 100), ylab = axis.label("A", prefix="k"),
#          col = rainbow(8), lwd = 2, bg = "slategray3")
#  abline(h = 0, lty = 2, lwd = 2)

## ----PPM_basis, results="hide", message=FALSE---------------------------------
basis(c("FeS2", "H2S", "O2", "H2O"))
species(c("pyrite", "magnetite"))
species("pyrrhotite", "cr2")

## ----PPM_affinity, message=FALSE, echo=1--------------------------------------
unlist(affinity(T = 300, P = 100)$values)
## 2031 1999 2036 
##    0    0    0

## ----PPM_setup, results="hide"------------------------------------------------
mod.buffer("PPM", "pyrrhotite", "cr2")
basis(c("H2S", "O2"), c("PPM", "PPM"))

## ----PPM_activities, message=FALSE--------------------------------------------
unlist(affinity(T = 300, P = 100, return.buffer = TRUE)[1:3])

## ----demo_buffer_noecho, fig.margin=TRUE, fig.width=4, fig.height=4, dpi=dpi, out.width="100%", message=FALSE, echo=FALSE, cache=TRUE, fig.cap="Values of log<i>f</i><sub>H<sub>2</sub></sub> corresponding to mineral buffers or to given activities of aqueous species.", pngquant=pngquant, timeit=timeit----
demo(buffer, echo = FALSE)

## ----PPM_affinity, eval=FALSE-------------------------------------------------
#  unlist(affinity(T = 300, P = 100)$values)
#  ## 2031 1999 2036
#  ##    0    0    0

## ----demo_buffer, eval=FALSE--------------------------------------------------
#  demo(buffer)

## ----bjerrum_diagram, fig.margin=TRUE, fig.width=3, fig.height=6, dpi=dpi, out.width="100%", echo=FALSE, results="hide", message=FALSE, cache=TRUE, fig.cap="Three views of carbonate speciation: affinity, activity, degree of formation.", pngquant=pngquant, timeit=timeit----
par(mfrow = c(3, 1))
basis("CHNOS+")
species(c("CO2", "HCO3-", "CO3-2"))
a25 <- affinity(pH = c(4, 13))
a150 <- affinity(pH = c(4, 13), T = 150)
diagram(a25, dy = 0.4)
diagram(a150, add = TRUE, names = FALSE, col = "red")
e25 <- equilibrate(a25, loga.balance = -3)
e150 <- equilibrate(a150, loga.balance = -3)
diagram(e25, ylim = c(-6, 0), dy = 0.15)
diagram(e150, add = TRUE, names = FALSE, col = "red")
diagram(e25, alpha = TRUE, dy = -0.25)
diagram(e150, alpha = TRUE, add = TRUE, names = FALSE, col = "red")

## ----bjerrum_diagram, echo=1:7, eval=FALSE------------------------------------
#  par(mfrow = c(3, 1))
#  basis("CHNOS+")
#  species(c("CO2", "HCO3-", "CO3-2"))
#  a25 <- affinity(pH = c(4, 13))
#  a150 <- affinity(pH = c(4, 13), T = 150)
#  diagram(a25, dy = 0.4)
#  diagram(a150, add = TRUE, names = FALSE, col = "red")
#  e25 <- equilibrate(a25, loga.balance = -3)
#  e150 <- equilibrate(a150, loga.balance = -3)
#  diagram(e25, ylim = c(-6, 0), dy = 0.15)
#  diagram(e150, add = TRUE, names = FALSE, col = "red")
#  diagram(e25, alpha = TRUE, dy = -0.25)
#  diagram(e150, alpha = TRUE, add = TRUE, names = FALSE, col = "red")

## ----bjerrum_diagram, echo=8:11, eval=FALSE-----------------------------------
#  par(mfrow = c(3, 1))
#  basis("CHNOS+")
#  species(c("CO2", "HCO3-", "CO3-2"))
#  a25 <- affinity(pH = c(4, 13))
#  a150 <- affinity(pH = c(4, 13), T = 150)
#  diagram(a25, dy = 0.4)
#  diagram(a150, add = TRUE, names = FALSE, col = "red")
#  e25 <- equilibrate(a25, loga.balance = -3)
#  e150 <- equilibrate(a150, loga.balance = -3)
#  diagram(e25, ylim = c(-6, 0), dy = 0.15)
#  diagram(e150, add = TRUE, names = FALSE, col = "red")
#  diagram(e25, alpha = TRUE, dy = -0.25)
#  diagram(e150, alpha = TRUE, add = TRUE, names = FALSE, col = "red")

## ----bjerrum_diagram, echo=12:13, eval=FALSE----------------------------------
#  par(mfrow = c(3, 1))
#  basis("CHNOS+")
#  species(c("CO2", "HCO3-", "CO3-2"))
#  a25 <- affinity(pH = c(4, 13))
#  a150 <- affinity(pH = c(4, 13), T = 150)
#  diagram(a25, dy = 0.4)
#  diagram(a150, add = TRUE, names = FALSE, col = "red")
#  e25 <- equilibrate(a25, loga.balance = -3)
#  e150 <- equilibrate(a150, loga.balance = -3)
#  diagram(e25, ylim = c(-6, 0), dy = 0.15)
#  diagram(e150, add = TRUE, names = FALSE, col = "red")
#  diagram(e25, alpha = TRUE, dy = -0.25)
#  diagram(e150, alpha = TRUE, add = TRUE, names = FALSE, col = "red")

## ----corundum, fig.margin=TRUE, fig.width=4, fig.height=4, dpi=dpi, out.width="100%", results="hide", message=FALSE, cache=TRUE, fig.cap="Solubility of corundum (green line) and equilibrium concentrations of aqueous species (black lines).", pngquant=pngquant, timeit=timeit----
add.OBIGT("SLOP98")
basis(c("Al+3", "H2O", "H+", "O2"))
species("corundum")
iaq <- c("Al+3", "AlO2-", "AlOH+2", "AlO+", "HAlO2")
s <- solubility(iaq, pH = c(0, 10), IS = 0, in.terms.of = "Al+3")
diagram(s, type = "loga.balance", ylim = c(-10, 0), lwd = 4, col = "green3")
diagram(s, add = TRUE, adj = c(0, 1, 2.1, -0.2, -1.5), dy = c(0, 0, 4, -0.3, 0.1))
legend("topright", c("25 °C", "1 bar"), text.font = 2, bty = "n")
reset()

## ----Alberty------------------------------------------------------------------
oldnon <- nonideal("Alberty")

## ----subcrt_IS----------------------------------------------------------------
subcrt(c("MgATP-2", "MgHATP-", "MgH2ATP"),
       T = c(25, 100), IS = c(0, 0.25), property = "G")$out

## ----info_ATP, results="hide"-------------------------------------------------
info(" ATP")

## ----T_100--------------------------------------------------------------------
T <- 100

## ----ATP, eval=FALSE, echo=2:6------------------------------------------------
#  par(mfrow = c(1, 4), mar = c(3.1, 3.6, 2.1, 1.6), mgp = c(1.8, 0.5, 0))
#  basis("MgCHNOPS+")
#  species(c("ATP-4", "HATP-3", "H2ATP-2", "H3ATP-", "H4ATP"))
#  a <- affinity(pH = c(3, 9), T = T)
#  e <- equilibrate(a)
#  d <- diagram(e, alpha = TRUE, tplot = FALSE)
#  title(main = describe.property("T", T))
#  alphas <- do.call(rbind, d$plotvals)
#  nH <- alphas * 0:4
#  Hlab <- substitute(italic(N)[H^`+`])
#  plot(a$vals[[1]], colSums(nH), type = "l", xlab = "pH", ylab=Hlab, lty=2, col=2)
#  a <- affinity(pH = c(3, 9), IS = 0.25, T = T)
#  e <- equilibrate(a)
#  d <- diagram(e, alpha = TRUE, plot.it = FALSE)
#  alphas <- do.call(rbind, d$plotvals)
#  nH <- alphas * 0:4
#  lines(a$vals[[1]], colSums(nH))
#  legend("topright", legend = c("I = 0 M", "I = 0.25 M"), lty = 2:1, col = 2:1, cex = 0.8)
#  ATP.H <- substitute("ATP and H"^`+`)
#  title(main = ATP.H)
#  species(c("MgATP-2", "MgHATP-", "MgH2ATP", "Mg2ATP"), add = TRUE)
#  Hplot <- function(pMg, IS = 0.25) {
#    basis("Mg+2", -pMg)
#    a <- affinity(pH = c(3, 9), IS = IS, T = T)
#    e <- equilibrate(a)
#    d <- diagram(e, alpha = TRUE, plot.it = FALSE)
#    alphas <- do.call(rbind, d$plotvals)
#    NH <- alphas * c(0:4, 0, 1, 2, 0)
#    lines(a$vals[[1]], colSums(NH), lty = 7 - pMg, col = 7 - pMg)
#  }
#  plot(c(3, 9), c(0, 2), type = "n", xlab = "pH", ylab = Hlab)
#  lapply(2:6, Hplot)
#  legend("topright", legend = paste("pMg = ", 2:6), lty = 5:1, col = 5:1, cex = 0.8)
#  ATP.H.Mg <- substitute("ATP and H"^`+`~"and Mg"^`+2`)
#  title(main = ATP.H.Mg)
#  Mgplot <- function(pH, IS = 0.25) {
#    basis("pH", pH)
#    a <- affinity(`Mg+2` = c(-2, -7), IS = IS, T = T)
#    e <- equilibrate(a)
#    d <- diagram(e, alpha = TRUE, plot.it = FALSE)
#    alphas <- do.call(rbind, d$plotvals)
#    NMg <- alphas * species()$`Mg+`
#    lines(-a$vals[[1]], colSums(NMg), lty = 10 - pH, col = 10 - pH)
#  }
#  Mglab <- substitute(italic(N)[Mg^`+2`])
#  plot(c(2, 7), c(0, 1.2), type = "n", xlab = "pMg", ylab = Mglab)
#  lapply(3:9, Mgplot)
#  legend("topright", legend = paste("pH = ", 3:9), lty = 7:1, col = 7:1, cex = 0.8)
#  title(main = ATP.H.Mg)

## ----ATP, eval=FALSE, echo=8:11-----------------------------------------------
#  par(mfrow = c(1, 4), mar = c(3.1, 3.6, 2.1, 1.6), mgp = c(1.8, 0.5, 0))
#  basis("MgCHNOPS+")
#  species(c("ATP-4", "HATP-3", "H2ATP-2", "H3ATP-", "H4ATP"))
#  a <- affinity(pH = c(3, 9), T = T)
#  e <- equilibrate(a)
#  d <- diagram(e, alpha = TRUE, tplot = FALSE)
#  title(main = describe.property("T", T))
#  alphas <- do.call(rbind, d$plotvals)
#  nH <- alphas * 0:4
#  Hlab <- substitute(italic(N)[H^`+`])
#  plot(a$vals[[1]], colSums(nH), type = "l", xlab = "pH", ylab=Hlab, lty=2, col=2)
#  a <- affinity(pH = c(3, 9), IS = 0.25, T = T)
#  e <- equilibrate(a)
#  d <- diagram(e, alpha = TRUE, plot.it = FALSE)
#  alphas <- do.call(rbind, d$plotvals)
#  nH <- alphas * 0:4
#  lines(a$vals[[1]], colSums(nH))
#  legend("topright", legend = c("I = 0 M", "I = 0.25 M"), lty = 2:1, col = 2:1, cex = 0.8)
#  ATP.H <- substitute("ATP and H"^`+`)
#  title(main = ATP.H)
#  species(c("MgATP-2", "MgHATP-", "MgH2ATP", "Mg2ATP"), add = TRUE)
#  Hplot <- function(pMg, IS = 0.25) {
#    basis("Mg+2", -pMg)
#    a <- affinity(pH = c(3, 9), IS = IS, T = T)
#    e <- equilibrate(a)
#    d <- diagram(e, alpha = TRUE, plot.it = FALSE)
#    alphas <- do.call(rbind, d$plotvals)
#    NH <- alphas * c(0:4, 0, 1, 2, 0)
#    lines(a$vals[[1]], colSums(NH), lty = 7 - pMg, col = 7 - pMg)
#  }
#  plot(c(3, 9), c(0, 2), type = "n", xlab = "pH", ylab = Hlab)
#  lapply(2:6, Hplot)
#  legend("topright", legend = paste("pMg = ", 2:6), lty = 5:1, col = 5:1, cex = 0.8)
#  ATP.H.Mg <- substitute("ATP and H"^`+`~"and Mg"^`+2`)
#  title(main = ATP.H.Mg)
#  Mgplot <- function(pH, IS = 0.25) {
#    basis("pH", pH)
#    a <- affinity(`Mg+2` = c(-2, -7), IS = IS, T = T)
#    e <- equilibrate(a)
#    d <- diagram(e, alpha = TRUE, plot.it = FALSE)
#    alphas <- do.call(rbind, d$plotvals)
#    NMg <- alphas * species()$`Mg+`
#    lines(-a$vals[[1]], colSums(NMg), lty = 10 - pH, col = 10 - pH)
#  }
#  Mglab <- substitute(italic(N)[Mg^`+2`])
#  plot(c(2, 7), c(0, 1.2), type = "n", xlab = "pMg", ylab = Mglab)
#  lapply(3:9, Mgplot)
#  legend("topright", legend = paste("pH = ", 3:9), lty = 7:1, col = 7:1, cex = 0.8)
#  title(main = ATP.H.Mg)

## ----ATP, eval=FALSE, echo=12:17----------------------------------------------
#  par(mfrow = c(1, 4), mar = c(3.1, 3.6, 2.1, 1.6), mgp = c(1.8, 0.5, 0))
#  basis("MgCHNOPS+")
#  species(c("ATP-4", "HATP-3", "H2ATP-2", "H3ATP-", "H4ATP"))
#  a <- affinity(pH = c(3, 9), T = T)
#  e <- equilibrate(a)
#  d <- diagram(e, alpha = TRUE, tplot = FALSE)
#  title(main = describe.property("T", T))
#  alphas <- do.call(rbind, d$plotvals)
#  nH <- alphas * 0:4
#  Hlab <- substitute(italic(N)[H^`+`])
#  plot(a$vals[[1]], colSums(nH), type = "l", xlab = "pH", ylab=Hlab, lty=2, col=2)
#  a <- affinity(pH = c(3, 9), IS = 0.25, T = T)
#  e <- equilibrate(a)
#  d <- diagram(e, alpha = TRUE, plot.it = FALSE)
#  alphas <- do.call(rbind, d$plotvals)
#  nH <- alphas * 0:4
#  lines(a$vals[[1]], colSums(nH))
#  legend("topright", legend = c("I = 0 M", "I = 0.25 M"), lty = 2:1, col = 2:1, cex = 0.8)
#  ATP.H <- substitute("ATP and H"^`+`)
#  title(main = ATP.H)
#  species(c("MgATP-2", "MgHATP-", "MgH2ATP", "Mg2ATP"), add = TRUE)
#  Hplot <- function(pMg, IS = 0.25) {
#    basis("Mg+2", -pMg)
#    a <- affinity(pH = c(3, 9), IS = IS, T = T)
#    e <- equilibrate(a)
#    d <- diagram(e, alpha = TRUE, plot.it = FALSE)
#    alphas <- do.call(rbind, d$plotvals)
#    NH <- alphas * c(0:4, 0, 1, 2, 0)
#    lines(a$vals[[1]], colSums(NH), lty = 7 - pMg, col = 7 - pMg)
#  }
#  plot(c(3, 9), c(0, 2), type = "n", xlab = "pH", ylab = Hlab)
#  lapply(2:6, Hplot)
#  legend("topright", legend = paste("pMg = ", 2:6), lty = 5:1, col = 5:1, cex = 0.8)
#  ATP.H.Mg <- substitute("ATP and H"^`+`~"and Mg"^`+2`)
#  title(main = ATP.H.Mg)
#  Mgplot <- function(pH, IS = 0.25) {
#    basis("pH", pH)
#    a <- affinity(`Mg+2` = c(-2, -7), IS = IS, T = T)
#    e <- equilibrate(a)
#    d <- diagram(e, alpha = TRUE, plot.it = FALSE)
#    alphas <- do.call(rbind, d$plotvals)
#    NMg <- alphas * species()$`Mg+`
#    lines(-a$vals[[1]], colSums(NMg), lty = 10 - pH, col = 10 - pH)
#  }
#  Mglab <- substitute(italic(N)[Mg^`+2`])
#  plot(c(2, 7), c(0, 1.2), type = "n", xlab = "pMg", ylab = Mglab)
#  lapply(3:9, Mgplot)
#  legend("topright", legend = paste("pH = ", 3:9), lty = 7:1, col = 7:1, cex = 0.8)
#  title(main = ATP.H.Mg)

## ----ATP, eval=FALSE, echo=21-------------------------------------------------
#  par(mfrow = c(1, 4), mar = c(3.1, 3.6, 2.1, 1.6), mgp = c(1.8, 0.5, 0))
#  basis("MgCHNOPS+")
#  species(c("ATP-4", "HATP-3", "H2ATP-2", "H3ATP-", "H4ATP"))
#  a <- affinity(pH = c(3, 9), T = T)
#  e <- equilibrate(a)
#  d <- diagram(e, alpha = TRUE, tplot = FALSE)
#  title(main = describe.property("T", T))
#  alphas <- do.call(rbind, d$plotvals)
#  nH <- alphas * 0:4
#  Hlab <- substitute(italic(N)[H^`+`])
#  plot(a$vals[[1]], colSums(nH), type = "l", xlab = "pH", ylab=Hlab, lty=2, col=2)
#  a <- affinity(pH = c(3, 9), IS = 0.25, T = T)
#  e <- equilibrate(a)
#  d <- diagram(e, alpha = TRUE, plot.it = FALSE)
#  alphas <- do.call(rbind, d$plotvals)
#  nH <- alphas * 0:4
#  lines(a$vals[[1]], colSums(nH))
#  legend("topright", legend = c("I = 0 M", "I = 0.25 M"), lty = 2:1, col = 2:1, cex = 0.8)
#  ATP.H <- substitute("ATP and H"^`+`)
#  title(main = ATP.H)
#  species(c("MgATP-2", "MgHATP-", "MgH2ATP", "Mg2ATP"), add = TRUE)
#  Hplot <- function(pMg, IS = 0.25) {
#    basis("Mg+2", -pMg)
#    a <- affinity(pH = c(3, 9), IS = IS, T = T)
#    e <- equilibrate(a)
#    d <- diagram(e, alpha = TRUE, plot.it = FALSE)
#    alphas <- do.call(rbind, d$plotvals)
#    NH <- alphas * c(0:4, 0, 1, 2, 0)
#    lines(a$vals[[1]], colSums(NH), lty = 7 - pMg, col = 7 - pMg)
#  }
#  plot(c(3, 9), c(0, 2), type = "n", xlab = "pH", ylab = Hlab)
#  lapply(2:6, Hplot)
#  legend("topright", legend = paste("pMg = ", 2:6), lty = 5:1, col = 5:1, cex = 0.8)
#  ATP.H.Mg <- substitute("ATP and H"^`+`~"and Mg"^`+2`)
#  title(main = ATP.H.Mg)
#  Mgplot <- function(pH, IS = 0.25) {
#    basis("pH", pH)
#    a <- affinity(`Mg+2` = c(-2, -7), IS = IS, T = T)
#    e <- equilibrate(a)
#    d <- diagram(e, alpha = TRUE, plot.it = FALSE)
#    alphas <- do.call(rbind, d$plotvals)
#    NMg <- alphas * species()$`Mg+`
#    lines(-a$vals[[1]], colSums(NMg), lty = 10 - pH, col = 10 - pH)
#  }
#  Mglab <- substitute(italic(N)[Mg^`+2`])
#  plot(c(2, 7), c(0, 1.2), type = "n", xlab = "pMg", ylab = Mglab)
#  lapply(3:9, Mgplot)
#  legend("topright", legend = paste("pH = ", 3:9), lty = 7:1, col = 7:1, cex = 0.8)
#  title(main = ATP.H.Mg)

## ----ATP, eval=FALSE, echo=22:30----------------------------------------------
#  par(mfrow = c(1, 4), mar = c(3.1, 3.6, 2.1, 1.6), mgp = c(1.8, 0.5, 0))
#  basis("MgCHNOPS+")
#  species(c("ATP-4", "HATP-3", "H2ATP-2", "H3ATP-", "H4ATP"))
#  a <- affinity(pH = c(3, 9), T = T)
#  e <- equilibrate(a)
#  d <- diagram(e, alpha = TRUE, tplot = FALSE)
#  title(main = describe.property("T", T))
#  alphas <- do.call(rbind, d$plotvals)
#  nH <- alphas * 0:4
#  Hlab <- substitute(italic(N)[H^`+`])
#  plot(a$vals[[1]], colSums(nH), type = "l", xlab = "pH", ylab=Hlab, lty=2, col=2)
#  a <- affinity(pH = c(3, 9), IS = 0.25, T = T)
#  e <- equilibrate(a)
#  d <- diagram(e, alpha = TRUE, plot.it = FALSE)
#  alphas <- do.call(rbind, d$plotvals)
#  nH <- alphas * 0:4
#  lines(a$vals[[1]], colSums(nH))
#  legend("topright", legend = c("I = 0 M", "I = 0.25 M"), lty = 2:1, col = 2:1, cex = 0.8)
#  ATP.H <- substitute("ATP and H"^`+`)
#  title(main = ATP.H)
#  species(c("MgATP-2", "MgHATP-", "MgH2ATP", "Mg2ATP"), add = TRUE)
#  Hplot <- function(pMg, IS = 0.25) {
#    basis("Mg+2", -pMg)
#    a <- affinity(pH = c(3, 9), IS = IS, T = T)
#    e <- equilibrate(a)
#    d <- diagram(e, alpha = TRUE, plot.it = FALSE)
#    alphas <- do.call(rbind, d$plotvals)
#    NH <- alphas * c(0:4, 0, 1, 2, 0)
#    lines(a$vals[[1]], colSums(NH), lty = 7 - pMg, col = 7 - pMg)
#  }
#  plot(c(3, 9), c(0, 2), type = "n", xlab = "pH", ylab = Hlab)
#  lapply(2:6, Hplot)
#  legend("topright", legend = paste("pMg = ", 2:6), lty = 5:1, col = 5:1, cex = 0.8)
#  ATP.H.Mg <- substitute("ATP and H"^`+`~"and Mg"^`+2`)
#  title(main = ATP.H.Mg)
#  Mgplot <- function(pH, IS = 0.25) {
#    basis("pH", pH)
#    a <- affinity(`Mg+2` = c(-2, -7), IS = IS, T = T)
#    e <- equilibrate(a)
#    d <- diagram(e, alpha = TRUE, plot.it = FALSE)
#    alphas <- do.call(rbind, d$plotvals)
#    NMg <- alphas * species()$`Mg+`
#    lines(-a$vals[[1]], colSums(NMg), lty = 10 - pH, col = 10 - pH)
#  }
#  Mglab <- substitute(italic(N)[Mg^`+2`])
#  plot(c(2, 7), c(0, 1.2), type = "n", xlab = "pMg", ylab = Mglab)
#  lapply(3:9, Mgplot)
#  legend("topright", legend = paste("pH = ", 3:9), lty = 7:1, col = 7:1, cex = 0.8)
#  title(main = ATP.H.Mg)

## ----ATP, eval=FALSE, echo=31:32----------------------------------------------
#  par(mfrow = c(1, 4), mar = c(3.1, 3.6, 2.1, 1.6), mgp = c(1.8, 0.5, 0))
#  basis("MgCHNOPS+")
#  species(c("ATP-4", "HATP-3", "H2ATP-2", "H3ATP-", "H4ATP"))
#  a <- affinity(pH = c(3, 9), T = T)
#  e <- equilibrate(a)
#  d <- diagram(e, alpha = TRUE, tplot = FALSE)
#  title(main = describe.property("T", T))
#  alphas <- do.call(rbind, d$plotvals)
#  nH <- alphas * 0:4
#  Hlab <- substitute(italic(N)[H^`+`])
#  plot(a$vals[[1]], colSums(nH), type = "l", xlab = "pH", ylab=Hlab, lty=2, col=2)
#  a <- affinity(pH = c(3, 9), IS = 0.25, T = T)
#  e <- equilibrate(a)
#  d <- diagram(e, alpha = TRUE, plot.it = FALSE)
#  alphas <- do.call(rbind, d$plotvals)
#  nH <- alphas * 0:4
#  lines(a$vals[[1]], colSums(nH))
#  legend("topright", legend = c("I = 0 M", "I = 0.25 M"), lty = 2:1, col = 2:1, cex = 0.8)
#  ATP.H <- substitute("ATP and H"^`+`)
#  title(main = ATP.H)
#  species(c("MgATP-2", "MgHATP-", "MgH2ATP", "Mg2ATP"), add = TRUE)
#  Hplot <- function(pMg, IS = 0.25) {
#    basis("Mg+2", -pMg)
#    a <- affinity(pH = c(3, 9), IS = IS, T = T)
#    e <- equilibrate(a)
#    d <- diagram(e, alpha = TRUE, plot.it = FALSE)
#    alphas <- do.call(rbind, d$plotvals)
#    NH <- alphas * c(0:4, 0, 1, 2, 0)
#    lines(a$vals[[1]], colSums(NH), lty = 7 - pMg, col = 7 - pMg)
#  }
#  plot(c(3, 9), c(0, 2), type = "n", xlab = "pH", ylab = Hlab)
#  lapply(2:6, Hplot)
#  legend("topright", legend = paste("pMg = ", 2:6), lty = 5:1, col = 5:1, cex = 0.8)
#  ATP.H.Mg <- substitute("ATP and H"^`+`~"and Mg"^`+2`)
#  title(main = ATP.H.Mg)
#  Mgplot <- function(pH, IS = 0.25) {
#    basis("pH", pH)
#    a <- affinity(`Mg+2` = c(-2, -7), IS = IS, T = T)
#    e <- equilibrate(a)
#    d <- diagram(e, alpha = TRUE, plot.it = FALSE)
#    alphas <- do.call(rbind, d$plotvals)
#    NMg <- alphas * species()$`Mg+`
#    lines(-a$vals[[1]], colSums(NMg), lty = 10 - pH, col = 10 - pH)
#  }
#  Mglab <- substitute(italic(N)[Mg^`+2`])
#  plot(c(2, 7), c(0, 1.2), type = "n", xlab = "pMg", ylab = Mglab)
#  lapply(3:9, Mgplot)
#  legend("topright", legend = paste("pH = ", 3:9), lty = 7:1, col = 7:1, cex = 0.8)
#  title(main = ATP.H.Mg)

## ----ATP, eval=FALSE, echo=36:44----------------------------------------------
#  par(mfrow = c(1, 4), mar = c(3.1, 3.6, 2.1, 1.6), mgp = c(1.8, 0.5, 0))
#  basis("MgCHNOPS+")
#  species(c("ATP-4", "HATP-3", "H2ATP-2", "H3ATP-", "H4ATP"))
#  a <- affinity(pH = c(3, 9), T = T)
#  e <- equilibrate(a)
#  d <- diagram(e, alpha = TRUE, tplot = FALSE)
#  title(main = describe.property("T", T))
#  alphas <- do.call(rbind, d$plotvals)
#  nH <- alphas * 0:4
#  Hlab <- substitute(italic(N)[H^`+`])
#  plot(a$vals[[1]], colSums(nH), type = "l", xlab = "pH", ylab=Hlab, lty=2, col=2)
#  a <- affinity(pH = c(3, 9), IS = 0.25, T = T)
#  e <- equilibrate(a)
#  d <- diagram(e, alpha = TRUE, plot.it = FALSE)
#  alphas <- do.call(rbind, d$plotvals)
#  nH <- alphas * 0:4
#  lines(a$vals[[1]], colSums(nH))
#  legend("topright", legend = c("I = 0 M", "I = 0.25 M"), lty = 2:1, col = 2:1, cex = 0.8)
#  ATP.H <- substitute("ATP and H"^`+`)
#  title(main = ATP.H)
#  species(c("MgATP-2", "MgHATP-", "MgH2ATP", "Mg2ATP"), add = TRUE)
#  Hplot <- function(pMg, IS = 0.25) {
#    basis("Mg+2", -pMg)
#    a <- affinity(pH = c(3, 9), IS = IS, T = T)
#    e <- equilibrate(a)
#    d <- diagram(e, alpha = TRUE, plot.it = FALSE)
#    alphas <- do.call(rbind, d$plotvals)
#    NH <- alphas * c(0:4, 0, 1, 2, 0)
#    lines(a$vals[[1]], colSums(NH), lty = 7 - pMg, col = 7 - pMg)
#  }
#  plot(c(3, 9), c(0, 2), type = "n", xlab = "pH", ylab = Hlab)
#  lapply(2:6, Hplot)
#  legend("topright", legend = paste("pMg = ", 2:6), lty = 5:1, col = 5:1, cex = 0.8)
#  ATP.H.Mg <- substitute("ATP and H"^`+`~"and Mg"^`+2`)
#  title(main = ATP.H.Mg)
#  Mgplot <- function(pH, IS = 0.25) {
#    basis("pH", pH)
#    a <- affinity(`Mg+2` = c(-2, -7), IS = IS, T = T)
#    e <- equilibrate(a)
#    d <- diagram(e, alpha = TRUE, plot.it = FALSE)
#    alphas <- do.call(rbind, d$plotvals)
#    NMg <- alphas * species()$`Mg+`
#    lines(-a$vals[[1]], colSums(NMg), lty = 10 - pH, col = 10 - pH)
#  }
#  Mglab <- substitute(italic(N)[Mg^`+2`])
#  plot(c(2, 7), c(0, 1.2), type = "n", xlab = "pMg", ylab = Mglab)
#  lapply(3:9, Mgplot)
#  legend("topright", legend = paste("pH = ", 3:9), lty = 7:1, col = 7:1, cex = 0.8)
#  title(main = ATP.H.Mg)

## ----ATP, eval=FALSE, echo=45:47----------------------------------------------
#  par(mfrow = c(1, 4), mar = c(3.1, 3.6, 2.1, 1.6), mgp = c(1.8, 0.5, 0))
#  basis("MgCHNOPS+")
#  species(c("ATP-4", "HATP-3", "H2ATP-2", "H3ATP-", "H4ATP"))
#  a <- affinity(pH = c(3, 9), T = T)
#  e <- equilibrate(a)
#  d <- diagram(e, alpha = TRUE, tplot = FALSE)
#  title(main = describe.property("T", T))
#  alphas <- do.call(rbind, d$plotvals)
#  nH <- alphas * 0:4
#  Hlab <- substitute(italic(N)[H^`+`])
#  plot(a$vals[[1]], colSums(nH), type = "l", xlab = "pH", ylab=Hlab, lty=2, col=2)
#  a <- affinity(pH = c(3, 9), IS = 0.25, T = T)
#  e <- equilibrate(a)
#  d <- diagram(e, alpha = TRUE, plot.it = FALSE)
#  alphas <- do.call(rbind, d$plotvals)
#  nH <- alphas * 0:4
#  lines(a$vals[[1]], colSums(nH))
#  legend("topright", legend = c("I = 0 M", "I = 0.25 M"), lty = 2:1, col = 2:1, cex = 0.8)
#  ATP.H <- substitute("ATP and H"^`+`)
#  title(main = ATP.H)
#  species(c("MgATP-2", "MgHATP-", "MgH2ATP", "Mg2ATP"), add = TRUE)
#  Hplot <- function(pMg, IS = 0.25) {
#    basis("Mg+2", -pMg)
#    a <- affinity(pH = c(3, 9), IS = IS, T = T)
#    e <- equilibrate(a)
#    d <- diagram(e, alpha = TRUE, plot.it = FALSE)
#    alphas <- do.call(rbind, d$plotvals)
#    NH <- alphas * c(0:4, 0, 1, 2, 0)
#    lines(a$vals[[1]], colSums(NH), lty = 7 - pMg, col = 7 - pMg)
#  }
#  plot(c(3, 9), c(0, 2), type = "n", xlab = "pH", ylab = Hlab)
#  lapply(2:6, Hplot)
#  legend("topright", legend = paste("pMg = ", 2:6), lty = 5:1, col = 5:1, cex = 0.8)
#  ATP.H.Mg <- substitute("ATP and H"^`+`~"and Mg"^`+2`)
#  title(main = ATP.H.Mg)
#  Mgplot <- function(pH, IS = 0.25) {
#    basis("pH", pH)
#    a <- affinity(`Mg+2` = c(-2, -7), IS = IS, T = T)
#    e <- equilibrate(a)
#    d <- diagram(e, alpha = TRUE, plot.it = FALSE)
#    alphas <- do.call(rbind, d$plotvals)
#    NMg <- alphas * species()$`Mg+`
#    lines(-a$vals[[1]], colSums(NMg), lty = 10 - pH, col = 10 - pH)
#  }
#  Mglab <- substitute(italic(N)[Mg^`+2`])
#  plot(c(2, 7), c(0, 1.2), type = "n", xlab = "pMg", ylab = Mglab)
#  lapply(3:9, Mgplot)
#  legend("topright", legend = paste("pH = ", 3:9), lty = 7:1, col = 7:1, cex = 0.8)
#  title(main = ATP.H.Mg)

## ----ATP, fig.fullwidth=TRUE, fig.width=10, fig.height=2.5, dpi=ifelse(dpi==50, 50, 100), out.width="100%", echo=FALSE, message=FALSE, results="hide", fig.cap="Binding of H<sup>+</sup> and Mg<sup>+2</sup> to ATP at 100 °C and *I* = 0 M (first plot) or *I* = 0.25 M (third and fourth plots).", cache=TRUE, pngquant=pngquant, timeit=timeit----
par(mfrow = c(1, 4), mar = c(3.1, 3.6, 2.1, 1.6), mgp = c(1.8, 0.5, 0))
basis("MgCHNOPS+")
species(c("ATP-4", "HATP-3", "H2ATP-2", "H3ATP-", "H4ATP"))
a <- affinity(pH = c(3, 9), T = T)
e <- equilibrate(a)
d <- diagram(e, alpha = TRUE, tplot = FALSE)
title(main = describe.property("T", T))
alphas <- do.call(rbind, d$plotvals)
nH <- alphas * 0:4
Hlab <- substitute(italic(N)[H^`+`])
plot(a$vals[[1]], colSums(nH), type = "l", xlab = "pH", ylab=Hlab, lty=2, col=2)
a <- affinity(pH = c(3, 9), IS = 0.25, T = T)
e <- equilibrate(a)
d <- diagram(e, alpha = TRUE, plot.it = FALSE)
alphas <- do.call(rbind, d$plotvals)
nH <- alphas * 0:4
lines(a$vals[[1]], colSums(nH))
legend("topright", legend = c("I = 0 M", "I = 0.25 M"), lty = 2:1, col = 2:1, cex = 0.8)
ATP.H <- substitute("ATP and H"^`+`)
title(main = ATP.H)
species(c("MgATP-2", "MgHATP-", "MgH2ATP", "Mg2ATP"), add = TRUE)
Hplot <- function(pMg, IS = 0.25) {
  basis("Mg+2", -pMg)
  a <- affinity(pH = c(3, 9), IS = IS, T = T)
  e <- equilibrate(a)
  d <- diagram(e, alpha = TRUE, plot.it = FALSE)
  alphas <- do.call(rbind, d$plotvals)
  NH <- alphas * c(0:4, 0, 1, 2, 0)
  lines(a$vals[[1]], colSums(NH), lty = 7 - pMg, col = 7 - pMg)
}
plot(c(3, 9), c(0, 2), type = "n", xlab = "pH", ylab = Hlab)
lapply(2:6, Hplot)
legend("topright", legend = paste("pMg = ", 2:6), lty = 5:1, col = 5:1, cex = 0.8)
ATP.H.Mg <- substitute("ATP and H"^`+`~"and Mg"^`+2`)
title(main = ATP.H.Mg)
Mgplot <- function(pH, IS = 0.25) {
  basis("pH", pH)
  a <- affinity(`Mg+2` = c(-2, -7), IS = IS, T = T)
  e <- equilibrate(a)
  d <- diagram(e, alpha = TRUE, plot.it = FALSE)
  alphas <- do.call(rbind, d$plotvals)
  NMg <- alphas * species()$`Mg+`
  lines(-a$vals[[1]], colSums(NMg), lty = 10 - pH, col = 10 - pH)
}
Mglab <- substitute(italic(N)[Mg^`+2`])
plot(c(2, 7), c(0, 1.2), type = "n", xlab = "pMg", ylab = Mglab)
lapply(3:9, Mgplot)
legend("topright", legend = paste("pH = ", 3:9), lty = 7:1, col = 7:1, cex = 0.8)
title(main = ATP.H.Mg)

## ----oldnon-------------------------------------------------------------------
nonideal(oldnon)

## ----pinfo_LYSC_CHICK---------------------------------------------------------
p1 <- pinfo("LYSC_CHICK")
p2 <- pinfo(c("SHH", "OLIG2"), "HUMAN")
pinfo(c(p1, p2))

## ----formula_LYSC_CHICK-------------------------------------------------------
pl <- protein.length("LYSC_CHICK")
pf <- protein.formula("LYSC_CHICK")
list(length = pl, protein = pf, residue = pf / pl,
     ZC_protein = ZC(pf), ZC_residue = ZC(pf / pl))

## ----subcrt_LYSC_CHICK, message=FALSE-----------------------------------------
subcrt("LYSC_CHICK")$out[[1]][1:6, ]

## ----protein_Cp, fig.margin=TRUE, fig.width=4, fig.height=4, small.mar=TRUE, dpi=dpi, out.width="100%", echo=FALSE, message=FALSE, fig.cap='The heat capacity calculated by group additivity closely approximates experimental values for aqueous proteins. For a related figure showing the effects of ionization in the calculations, see <span style="color:blue">`?ionize.aa`</span>.', cache=TRUE, pngquant=pngquant, timeit=timeit----
PM90 <- read.csv(system.file("extdata/cpetc/PM90.csv", package = "CHNOSZ"))
plength <- protein.length(colnames(PM90)[2:5])
Cp_expt <- t(t(PM90[, 2:5]) / plength)
matplot(PM90[, 1], convert(Cp_expt, "cal"), type = "p", pch = 19,
        xlab = axis.label("T"), ylab = axis.label("Cp0"), ylim = c(28, 65))
for(i in 1:4) {
  pname <- colnames(Cp_expt)[i]
  aq <- subcrt(pname, "aq", T = seq(0, 150))$out[[1]]
  cr <- subcrt(pname, "cr", T = seq(0, 150))$out[[1]]
  lines(aq$T, aq$Cp / plength[i], col = i)
  lines(cr$T, cr$Cp / plength[i], col = i, lty = 2)
}
legend("right", legend = colnames(Cp_expt),
       col = 1:4, pch = 19, lty = 1, bty = "n", cex = 0.9)
legend("bottomright", legend = c("experimental", "calculated (aq)",
       "calculated (cr)"), lty = c(NA, 1, 2), pch = c(19, NA, NA), bty = "n")

## ----protein_Cp, eval=FALSE, echo=1:5-----------------------------------------
#  PM90 <- read.csv(system.file("extdata/cpetc/PM90.csv", package = "CHNOSZ"))
#  plength <- protein.length(colnames(PM90)[2:5])
#  Cp_expt <- t(t(PM90[, 2:5]) / plength)
#  matplot(PM90[, 1], convert(Cp_expt, "cal"), type = "p", pch = 19,
#          xlab = axis.label("T"), ylab = axis.label("Cp0"), ylim = c(28, 65))
#  for(i in 1:4) {
#    pname <- colnames(Cp_expt)[i]
#    aq <- subcrt(pname, "aq", T = seq(0, 150))$out[[1]]
#    cr <- subcrt(pname, "cr", T = seq(0, 150))$out[[1]]
#    lines(aq$T, aq$Cp / plength[i], col = i)
#    lines(cr$T, cr$Cp / plength[i], col = i, lty = 2)
#  }
#  legend("right", legend = colnames(Cp_expt),
#         col = 1:4, pch = 19, lty = 1, bty = "n", cex = 0.9)
#  legend("bottomright", legend = c("experimental", "calculated (aq)",
#         "calculated (cr)"), lty = c(NA, 1, 2), pch = c(19, NA, NA), bty = "n")

## ----protein_Cp, eval=FALSE, echo=-(1:5)--------------------------------------
#  PM90 <- read.csv(system.file("extdata/cpetc/PM90.csv", package = "CHNOSZ"))
#  plength <- protein.length(colnames(PM90)[2:5])
#  Cp_expt <- t(t(PM90[, 2:5]) / plength)
#  matplot(PM90[, 1], convert(Cp_expt, "cal"), type = "p", pch = 19,
#          xlab = axis.label("T"), ylab = axis.label("Cp0"), ylim = c(28, 65))
#  for(i in 1:4) {
#    pname <- colnames(Cp_expt)[i]
#    aq <- subcrt(pname, "aq", T = seq(0, 150))$out[[1]]
#    cr <- subcrt(pname, "cr", T = seq(0, 150))$out[[1]]
#    lines(aq$T, aq$Cp / plength[i], col = i)
#    lines(cr$T, cr$Cp / plength[i], col = i, lty = 2)
#  }
#  legend("right", legend = colnames(Cp_expt),
#         col = 1:4, pch = 19, lty = 1, bty = "n", cex = 0.9)
#  legend("bottomright", legend = c("experimental", "calculated (aq)",
#         "calculated (cr)"), lty = c(NA, 1, 2), pch = c(19, NA, NA), bty = "n")

## ----protein_ionization, fig.margin=TRUE, fig.width=4, fig.height=4, small.mar=TRUE, dpi=dpi, out.width="100%", echo=FALSE, results="hide", message=FALSE, fig.cap='Affinity of ionization of proteins. See [<span style="color:blue">demo(ionize)</span>](../demo) for ionization properties calculated as a function of temperature and pH.', cache=TRUE, pngquant=pngquant, timeit=timeit----
ip <- pinfo(c("CYC_BOVIN", "LYSC_CHICK", "MYG_PHYCA", "RNAS1_BOVIN"))
basis("CHNOS+")
a_ion <- affinity(pH = c(0, 14), iprotein = ip)
basis("CHNOS")
a_nonion <- affinity(iprotein = ip)
plot(c(0, 14), c(50, 300), xlab = "pH", ylab = axis.label("A"), type = "n")
for(i in 1:4) {
  A_ion <- as.numeric(a_ion$values[[i]])
  A_nonion <- as.numeric(a_nonion$values[[i]])
  lines(a_ion$vals[[1]], A_ion - A_nonion, col=i)
}
legend("topright", legend = a_ion$species$name,
       col = 1:4, lty = 1, bty = "n", cex = 0.9)

## ----protein_ionization, eval=FALSE-------------------------------------------
#  ip <- pinfo(c("CYC_BOVIN", "LYSC_CHICK", "MYG_PHYCA", "RNAS1_BOVIN"))
#  basis("CHNOS+")
#  a_ion <- affinity(pH = c(0, 14), iprotein = ip)
#  basis("CHNOS")
#  a_nonion <- affinity(iprotein = ip)
#  plot(c(0, 14), c(50, 300), xlab = "pH", ylab = axis.label("A"), type = "n")
#  for(i in 1:4) {
#    A_ion <- as.numeric(a_ion$values[[i]])
#    A_nonion <- as.numeric(a_nonion$values[[i]])
#    lines(a_ion$vals[[1]], A_ion - A_nonion, col=i)
#  }
#  legend("topright", legend = a_ion$species$name,
#         col = 1:4, lty = 1, bty = "n", cex = 0.9)

## ----basis_CHNOS, echo=FALSE, results="hide"----------------------------------
basis("CHNOS")

## ----species_protein, results="hide", message=FALSE, echo=1:2-----------------
species(c("CYC_BOVIN", "LYSC_CHICK", "MYG_PHYCA", "RNAS1_BOVIN"))
a_nonion_species <- affinity()

## ----nonion_species_values----------------------------------------------------
unlist(a_nonion_species$values)

## ----nonion_values------------------------------------------------------------
unlist(a_nonion$values)

## ----rubisco_svg, echo=FALSE, results="hide", fig.margin=TRUE, fig.width=4, fig.height=4, small.mar=TRUE, out.width="100%", fig.ext='svg', custom.plot=TRUE, embed.tag=TRUE, fig.cap='Average oxidation state of carbon in Rubisco compared with optimal growth temperature of organisms. **This is an interactive image.** Move the mouse over the points to show the names of the organisms, and click to open a reference in a new window. (Made with [**RSVGTipsDevice**](https://cran.r-project.org/package=RSVGTipsDevice) using code that can be found in the source of this document.)'----
## copies the premade SVG image to the knitr figure path
file.copy("rubisco.svg", fig_path(".svg"))
## the code for making the SVG image -- not "live" in the vignette because RSVGTipsDevice isn't available on Windows
#if(require(RSVGTipsDevice)) {
#  datfile <- system.file("extdata/cpetc/rubisco.csv", package = "CHNOSZ")
#  fastafile <- system.file("extdata/protein/rubisco.fasta", package = "CHNOSZ")
#  dat <- read.csv(datfile)
#  aa <- read.fasta(fastafile)
#  Topt <- (dat$T1 + dat$T2) / 2
#  idat <- match(dat$ID, substr(aa$protein, 4, 9))
#  aa <- aa[idat, ]
#  ZC <- ZC(protein.formula(aa))
#  pch <- match(dat$domain, c("E", "B", "A")) - 1
#  col <- match(dat$domain, c("A", "B", "E")) + 1
#  # because the tooltip titles in the SVG file are shown by recent browsers,
#  # we do not need to draw the tooltips explicitly, so set toolTipMode=0
#  devSVGTips("rubisco.svg", toolTipMode=0, title="Rubisco")
#  par(cex=1.4)
#  # unfortunately, plotmath can't be used with devSVGTips,
#  # so axis labels here don't contain italics.
#  plot(Topt, ZC, type="n", xlab="T, &#176;C", ylab="ZC")
#  n <- rep(1:9, 3)
#  for(i in seq_along(Topt)) {
#    # adjust cex to make the symbols look the same size
#    cex <- ifelse(pch[i]==1, 2.5, 3.5)
#    points(Topt[i], ZC[i], pch=pch[i], cex=cex, col=col[i])
#    URL <- dat$URL[i]
#    setSVGShapeURL(URL, target="_blank")
#    setSVGShapeContents(paste0("<title>", dat$species[i], "</title>"))
#    text(Topt[i], ZC[i], n[i], cex = 1.2)
#  }
#  abline(v = c(36, 63), lty = 2, col = "grey")
#  legend("topright", legend = c("Archaea", "Bacteria", "Eukaryota"),
#         pch = c(2, 1, 0), col = 2:4, cex=1.5, pt.cex = c(3, 2.3, 3), bty="n")
#  dev.off()
#}

## ----rubisco_ZC, fig.keep="none", message=FALSE-------------------------------
datfile <- system.file("extdata/cpetc/rubisco.csv", package = "CHNOSZ")
fastafile <- system.file("extdata/protein/rubisco.fasta", package = "CHNOSZ")
dat <- read.csv(datfile)
aa <- read.fasta(fastafile)
Topt <- (dat$T1 + dat$T2) / 2
idat <- match(dat$ID, substr(aa$protein, 4, 9))
aa <- aa[idat, ]
ZC <- ZC(protein.formula(aa))
pch <- match(dat$domain, c("E", "B", "A")) - 1
col <- match(dat$domain, c("A", "B", "E")) + 1
plot(Topt, ZC, pch = pch, cex = 2, col = col,
     xlab = expression(list(italic(T)[opt], degree*C)),
     ylab = expression(italic(Z)[C]))
text(Topt, ZC, rep(1:9, 3), cex = 0.8)
abline(v = c(36, 63), lty = 2, col = "grey")
legend("topright", legend = c("Archaea", "Bacteria", "Eukaryota"),
       pch = c(2, 1, 0), col = 2:4, pt.cex = 2)

## ----rubisco_O2, fig.margin=TRUE, fig.width=4, fig.height=4, small.mar=TRUE, dpi=dpi, out.width="100%", echo=FALSE, results="hide", message=FALSE, fig.cap="Compositions of proteins projected into different sets of basis species.", cache=TRUE, pngquant=pngquant, timeit=timeit----
layout(matrix(1:4, nrow = 2))
par(mgp = c(1.8, 0.5, 0))
pl <- protein.length(aa)
ZClab <- axis.label("ZC")
nO2lab <- expression(bar(italic(n))[O[2]])
nH2Olab <- expression(bar(italic(n))[H[2]*O])
lapply(c("CHNOS", "QEC"), function(thisbasis) {
  basis(thisbasis)
  pb <- protein.basis(aa)
  nO2 <- pb[, "O2"] / pl
  plot(ZC, nO2, pch = pch, col = col, xlab = ZClab, ylab = nO2lab)
  nH2O <- pb[, "H2O"] / pl
  plot(ZC, nH2O, pch = pch, col = col, xlab = ZClab, ylab = nH2Olab)
  mtext(thisbasis, font = 2)
})

## ----rubisco_O2, eval=FALSE---------------------------------------------------
#  layout(matrix(1:4, nrow = 2))
#  par(mgp = c(1.8, 0.5, 0))
#  pl <- protein.length(aa)
#  ZClab <- axis.label("ZC")
#  nO2lab <- expression(bar(italic(n))[O[2]])
#  nH2Olab <- expression(bar(italic(n))[H[2]*O])
#  lapply(c("CHNOS", "QEC"), function(thisbasis) {
#    basis(thisbasis)
#    pb <- protein.basis(aa)
#    nO2 <- pb[, "O2"] / pl
#    plot(ZC, nO2, pch = pch, col = col, xlab = ZClab, ylab = nO2lab)
#    nH2O <- pb[, "H2O"] / pl
#    plot(ZC, nH2O, pch = pch, col = col, xlab = ZClab, ylab = nH2Olab)
#    mtext(thisbasis, font = 2)
#  })

## ----yeastgfp-----------------------------------------------------------------
protein <- c("YDL195W", "YHR098C", "YIL109C", "YLR208W", "YNL049C", "YPL085W")
abundance <- c(1840, 12200, NA, 21400, 1720, 358)
ina <- is.na(abundance)

## ----add_protein_yeast, message=FALSE-----------------------------------------
ip <- match(protein[!ina], thermo()$protein$protein)

## ----unitize------------------------------------------------------------------
pl <- protein.length(ip)
logact <- unitize(numeric(5), pl)
logabundance <- unitize(log10(abundance[!ina]), pl)

## ----yeastplot, eval=FALSE, echo=1:6------------------------------------------
#  par(mfrow = c(1, 3))
#  basis("CHNOS+")
#  a <- affinity(O2 = c(-80, -73), iprotein = ip, loga.protein = logact)
#  e <- equilibrate(a)
#  diagram(e, ylim = c(-5, -2), col = 1:5, lwd = 2)
#  e <- equilibrate(a, normalize = TRUE)
#  diagram(e, ylim = c(-5, -2.5), col = 1:5, lwd = 2)
#  abline(h = logabundance, lty = 1:5, col = 1:5)
#  revisit(e, "DGinf", logabundance)

## ----yeastplot, eval=FALSE, echo=7:9------------------------------------------
#  par(mfrow = c(1, 3))
#  basis("CHNOS+")
#  a <- affinity(O2 = c(-80, -73), iprotein = ip, loga.protein = logact)
#  e <- equilibrate(a)
#  diagram(e, ylim = c(-5, -2), col = 1:5, lwd = 2)
#  e <- equilibrate(a, normalize = TRUE)
#  diagram(e, ylim = c(-5, -2.5), col = 1:5, lwd = 2)
#  abline(h = logabundance, lty = 1:5, col = 1:5)
#  revisit(e, "DGinf", logabundance)

## ----yeastplot, eval=FALSE, echo=10-------------------------------------------
#  par(mfrow = c(1, 3))
#  basis("CHNOS+")
#  a <- affinity(O2 = c(-80, -73), iprotein = ip, loga.protein = logact)
#  e <- equilibrate(a)
#  diagram(e, ylim = c(-5, -2), col = 1:5, lwd = 2)
#  e <- equilibrate(a, normalize = TRUE)
#  diagram(e, ylim = c(-5, -2.5), col = 1:5, lwd = 2)
#  abline(h = logabundance, lty = 1:5, col = 1:5)
#  revisit(e, "DGinf", logabundance)

## ----yeastplot, fig.fullwidth=TRUE, fig.width=7.5, fig.height=2.5, dpi=ifelse(dpi==50, 50, 100), out.width="85%", echo=FALSE, message=FALSE, results="hide", cache=TRUE, fig.cap="ER-to-Golgi proteins: calculations without and with length normalization, and free energy difference between experimental and calculated abundances in metastable equilibrium with normalization.", pngquant=pngquant, timeit=timeit----
par(mfrow = c(1, 3))
basis("CHNOS+")
a <- affinity(O2 = c(-80, -73), iprotein = ip, loga.protein = logact)
e <- equilibrate(a)
diagram(e, ylim = c(-5, -2), col = 1:5, lwd = 2)
e <- equilibrate(a, normalize = TRUE)
diagram(e, ylim = c(-5, -2.5), col = 1:5, lwd = 2)
abline(h = logabundance, lty = 1:5, col = 1:5)
revisit(e, "DGinf", logabundance)

## ----read_csv-----------------------------------------------------------------
file <- system.file("extdata/protein/POLG.csv", package = "CHNOSZ")
aa_POLG <- read.csv(file, as.is = TRUE, nrows = 5)

## ----read_fasta, message=FALSE------------------------------------------------
file <- system.file("extdata/protein/EF-Tu.aln", package = "CHNOSZ")
aa_Ef <- read.fasta(file, iseq = 1:2)

## ----seq2aa-------------------------------------------------------------------
aa_PRIO <- seq2aa("PRIO_HUMAN", "
MANLGCWMLVLFVATWSDLGLCKKRPKPGGWNTGGSRYPGQGSPGGNRYPPQGGGGWGQP
HGGGWGQPHGGGWGQPHGGGWGQPHGGGWGQGGGTHSQWNKPSKPKTNMKHMAGAAAAGA
VVGGLGGYMLGSAMSRPIIHFGSDYEDRYYRENMHRYPNQVYYRPMDEYSNQNNFVHDCV
NITIKQHTVTTTTKGENFTETDVKMMERVVEQMCITQYERESQAYYQRGSSMVLFSSPPV
ILLISFLIFLIVG
")

## ----uniprot_aa, eval=FALSE---------------------------------------------------
#  IDs <- c("ALAT1_HUMAN", "P02452")
#  aa <- lapply(IDs, uniprot.aa)
#  ## uniprot.aa: trying http://www.uniprot.org/uniprot/ALAT1_HUMAN ... accession P24298 ...
#  ## >sp|P24298|ALAT1_HUMAN Alanine aminotransferase 1 OS=Homo sapiens GN=GPT PE=1 SV=3 (length 496)
#  ## uniprot.aa: trying http://www.uniprot.org/uniprot/P02452 ... accession P02452 ...
#  ## >sp|P02452|CO1A1_HUMAN Collagen alpha-1(I) chain OS=Homo sapiens GN=COL1A1 PE=1 SV=5 (length 1464)
#  aa_UniProt <- do.call(rbind, aa)

## ----uniprot_aa_offline, echo=FALSE-------------------------------------------
aa_ALAT1 <- seq2aa("ALAT1_HUMAN", "
MASSTGDRSQAVRHGLRAKVLTLDGMNPRVRRVEYAVRGPIVQRALELEQELRQGVKKPF
TEVIRANIGDAQAMGQRPITFLRQVLALCVNPDLLSSPNFPDDAKKRAERILQACGGHSL
GAYSVSSGIQLIREDVARYIERRDGGIPADPNNVFLSTGASDAIVTVLKLLVAGEGHTRT
GVLIPIPQYPLYSATLAELGAVQVDYYLDEERAWALDVAELHRALGQARDHCRPRALCVI
NPGNPTGQVQTRECIEAVIRFAFEERLFLLADEVYQDNVYAAGSQFHSFKKVLMEMGPPY
AGQQELASFHSTSKGYMGECGFRGGYVEVVNMDAAVQQQMLKLMSVRLCPPVPGQALLDL
VVSPPAPTDPSFAQFQAEKQAVLAELAAKAKLTEQVFNEAPGISCNPVQGAMYSFPRVQL
PPRAVERAQELGLAPDMFFCLRLLEETGICVVPGSGFGQREGTYHFRMTILPPLEKLRLL
LEKLSRFHAKFTLEYS
")
aa_CO1A1 <- seq2aa("CO1A1_HUMAN", "
MFSFVDLRLLLLLAATALLTHGQEEGQVEGQDEDIPPITCVQNGLRYHDRDVWKPEPCRI
CVCDNGKVLCDDVICDETKNCPGAEVPEGECCPVCPDGSESPTDQETTGVEGPKGDTGPR
GPRGPAGPPGRDGIPGQPGLPGPPGPPGPPGPPGLGGNFAPQLSYGYDEKSTGGISVPGP
MGPSGPRGLPGPPGAPGPQGFQGPPGEPGEPGASGPMGPRGPPGPPGKNGDDGEAGKPGR
PGERGPPGPQGARGLPGTAGLPGMKGHRGFSGLDGAKGDAGPAGPKGEPGSPGENGAPGQ
MGPRGLPGERGRPGAPGPAGARGNDGATGAAGPPGPTGPAGPPGFPGAVGAKGEAGPQGP
RGSEGPQGVRGEPGPPGPAGAAGPAGNPGADGQPGAKGANGAPGIAGAPGFPGARGPSGP
QGPGGPPGPKGNSGEPGAPGSKGDTGAKGEPGPVGVQGPPGPAGEEGKRGARGEPGPTGL
PGPPGERGGPGSRGFPGADGVAGPKGPAGERGSPGPAGPKGSPGEAGRPGEAGLPGAKGL
TGSPGSPGPDGKTGPPGPAGQDGRPGPPGPPGARGQAGVMGFPGPKGAAGEPGKAGERGV
PGPPGAVGPAGKDGEAGAQGPPGPAGPAGERGEQGPAGSPGFQGLPGPAGPPGEAGKPGE
QGVPGDLGAPGPSGARGERGFPGERGVQGPPGPAGPRGANGAPGNDGAKGDAGAPGAPGS
QGAPGLQGMPGERGAAGLPGPKGDRGDAGPKGADGSPGKDGVRGLTGPIGPPGPAGAPGD
KGESGPSGPAGPTGARGAPGDRGEPGPPGPAGFAGPPGADGQPGAKGEPGDAGAKGDAGP
PGPAGPAGPPGPIGNVGAPGAKGARGSAGPPGATGFPGAAGRVGPPGPSGNAGPPGPPGP
AGKEGGKGPRGETGPAGRPGEVGPPGPPGPAGEKGSPGADGPAGAPGTPGPQGIAGQRGV
VGLPGQRGERGFPGLPGPSGEPGKQGPSGASGERGPPGPMGPPGLAGPPGESGREGAPGA
EGSPGRDGSPGAKGDRGETGPAGPPGAPGAPGAPGPVGPAGKSGDRGETGPAGPTGPVGP
VGARGPAGPQGPRGDKGETGEQGDRGIKGHRGFSGLQGPPGPPGSPGEQGPSGASGPAGP
RGPPGSAGAPGKDGLNGLPGPIGPPGPRGRTGDAGPVGPPGPPGPPGPPGPPSAGFDFSF
LPQPPQEKAHDGGRYYRADDANVVRDRDLEVDTTLKSLSQQIENIRSPEGSRKNPARTCR
DLKMCHSDWKSGEYWIDPNQGCNLDAIKVFCNMETGETCVYPTQPSVAQKNWYISKNPKD
KRHVWFGESMTDGFQFEYGGQGSDPADVAIQLTFLRLMSTEASQNITYHCKNSVAYMDQQ
TGNLKKALLLQGSNEIEIRAEGNSRFTYSVTVDGCTSHTGAWGKTVIEYKTTKTSRLPII
DVAPLDVGAPDQEFGFDVGPVCFL
")
aa_UniProt <- rbind(aa_ALAT1, aa_CO1A1)
aa_UniProt$abbrv <- c("ALAT1", "CO1A1")

## ----aa_UniProt, cache=TRUE---------------------------------------------------
aa_UniProt

## ----protein_length-----------------------------------------------------------
myaa <- rbind(aa_Ef, aa_PRIO, aa_ALAT1)
protein.length(myaa)

## ----thermo_refs_table, eval=FALSE--------------------------------------------
#  thermo.refs()  ## shows table in a browser

## ----width180, include=FALSE------------------------------------------------------------------------------------------------------------------------------------------------------
options(width = 180)

## ----thermo_refs_numeric----------------------------------------------------------------------------------------------------------------------------------------------------------
iATP <- info("ATP-4")
iMgATP <- info("MgATP-2")
thermo.refs(c(iATP, iMgATP))

## ----thermo_refs_character--------------------------------------------------------------------------------------------------------------------------------------------------------
thermo.refs(c("HDNB78", "MGN03"))

## ----thermo_refs_subcrt, message=FALSE--------------------------------------------------------------------------------------------------------------------------------------------
substuff <- subcrt(c("C2H5OH", "O2", "CO2", "H2O"), c(-1, -3, 2, 3))
thermo.refs(substuff)

## ----width80, include=FALSE---------------------------------------------------
options(width = 80)

## ----thermo_refs_browse, eval=FALSE-------------------------------------------
#  iFo <- info("forsterite")
#  ref <- thermo.refs(iFo)
#  browseURL(ref$URL)  ## opens a link to worldcat.org

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

## ----mod_OBIGT_CoCl4_ghs------------------------------------------------------
mod.OBIGT("CoCl4-2", formula = "CoCl4-2", state = "aq", ref1 = "LBT+11",
  date = as.character(Sys.Date()), G = -134150, H = -171558, S = 19.55, Cp = 72.09, V = 27.74)

## ----mod_OBIGT_CoCl4_eos------------------------------------------------------
mod.OBIGT("CoCl4-2", a1 = 6.5467, a2 = 8.2069, a3 = 2.0130, a4 = -3.1183,
  c1 = 76.3357, c2 = 11.6389, omega = 2.9159, z = -2)

## ----CoCl4_reaction, message = FALSE, echo = 1:3------------------------------
T <- c(25, seq(50, 350, 50))
sres <- subcrt(c("Co+2", "Cl-", "CoCl4-2"), c(-1, -4, 1), T = T)
round(sres$out$logK, 2)
stopifnot(identical(round(sres$out$logK, 2), c(-3.2, -2.96, -2.02, -0.74, 0.77, 2.5, 4.57, 7.29)))

## ----mod_OBIGT_magnesiochromite_ghs-------------------------------------------
H <- -1762000
S <- 119.6
V <- 43.56
mod.OBIGT("magnesiochromite", formula = "MgCr2O4", state = "cr", ref1 = "KOSG00",
          date = as.character(Sys.Date()), E_units = "J", H = H, S = S, V = V)

## ----mod_OBIGT_magnesiochromite_eos-------------------------------------------
a <- 221.4
b <- -0.00102030 * 1e3
c <- -1757210 * 1e-5
d <- -1247.9
mod.OBIGT("magnesiochromite", E_units = "J", a = a, b = b, c = c, d = d,
          e = 0, f = 0, lambda = 0, T = 1500)

## ----subcrt_magnesiochromite--------------------------------------------------
T.units("K")
E.units("J")
subcrt("magnesiochromite", property = "Cp", T = c(250, 300, 340), P = 1)

## ----restore_units_magnesiochromite-------------------------------------------
T.units("C")
E.units("cal")

## ----info_CoCl4, results = "hide"---------------------------------------------
inew <- info("CoCl4-2")
info(inew)

## ----info_S3, results="hide"--------------------------------------------------
info(info("S3-"))

## ----check_OBIGT--------------------------------------------------------------
file <- system.file("extdata/adds/OBIGT_check.csv", package = "CHNOSZ")
dat <- read.csv(file, as.is = TRUE)
nrow(dat)

## ----citation_CHNOSZ, results="asis", echo = FALSE----------------------------
cref <- citation("CHNOSZ")
print(cref[1], style = "html")

## ----citation_multimetal, results="asis", echo = FALSE------------------------
print(cref[2], style = "html")

## ----maintainer_CHNOSZ--------------------------------------------------------
maintainer("CHNOSZ")

## ----file_edit_anintro, eval=FALSE--------------------------------------------
#  file.edit(system.file("doc/anintro.Rmd", package = "CHNOSZ"))

## ----the_end------------------------------------------------------------------
   ######    ##   ##    ##   ##    ######     #####  #####
 ##         ##===##    ## \\##   ##    ##     \\       //
 ######    ##   ##    ##   ##    ######    #####      #####

## ----sessionInfo--------------------------------------------------------------
sessionInfo()

