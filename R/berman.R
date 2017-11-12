# CHNOSZ/berman.R 20170930
# calculate thermodynamic properties of minerals using equations from:
#   Berman, R. G. (1988) Internally-consistent thermodynamic data for minerals
#      in the system Na2O-K2O-CaO-MgO-FeO-Fe2O3-Al2O3-SiO2-TiO2-H2O-CO2.
#      J. Petrol. 29, 445-522. https://doi.org/10.1093/petrology/29.2.445

berman <- function(name, T = 298.15, P = 1, thisinfo=NULL, check.G=FALSE, calc.transition=TRUE, calc.disorder=TRUE, units="cal") {
  # reference temperature and pressure
  Pr <- 1
  Tr <- 298.15
  # get T and P to be same length
  ncond <- max(length(T), length(P))
  T <- rep(T, length.out=ncond)
  P <- rep(P, length.out=ncond)
  # get thermodynamic parameters
  dir <- system.file("extdata/Berman/", package="CHNOSZ")
  Ber88 <- read.csv(paste0(dir, "/Ber88.csv"), as.is=TRUE)
  Ber90 <- read.csv(paste0(dir, "/Ber90.csv"), as.is=TRUE)
  SHD91 <- read.csv(paste0(dir, "/SHD91.csv"), as.is=TRUE)
  ZS92 <- read.csv(paste0(dir, "/ZS92.csv"), as.is=TRUE)
  JUN92 <- read.csv(paste0(dir, "/JUN92.csv"), as.is=TRUE)
  DS10 <- read.csv(paste0(dir, "/DS10.csv"), as.is=TRUE)
  FDM14 <- read.csv(paste0(dir, "/FDM+14.csv"), as.is=TRUE)
  BDat17 <- read.csv(paste0(dir, "/BDat17.csv"), as.is=TRUE)
  # assemble the files in reverse chronological order
  dat <- rbind(BDat17, FDM14, DS10, JUN92, ZS92, SHD91, Ber90, Ber88)
  # remove duplicates (only the first, i.e. latest entry is kept)
  dat <- dat[!duplicated(dat$name), ]
  # remove the multipliers
  multexp <- c(0, 0, 0, 0,          # Ber88 Table 2
               0, -2, -5, -7,             # Table 3a
               6, 12, 6, 10,              # Table 4
               0, 0, 0, 2, 5, 0,          # Table 3b
               0, 0, 0, -3, -5, 2, 6, -4  # Table 5
               )
  dat[, 2:27] <- t(t(dat[, 2:27]) / 10^multexp)
  # which row has data for this mineral?
  irow <- which(dat$name == name)
  # the function works fine with just the following assign() call,
  # but an explicit dummy assignment here is used to avoid "Undefined global functions or variables" in R CMD check
  GfPrTr <- HfPrTr <- SPrTr <- Tlambda <- Tmax <- Tmin <- Tref <- VPrTr <-
    d0 <- d1 <- d2 <- d3 <- d4 <- d5 <- dTdP <- k0 <- k1 <- k2 <- k3 <- l1 <- l2 <- v1 <- v2 <- v3 <- v4 <- NA
  # assign values to the variables used below
  for(i in 1:ncol(dat)) assign(colnames(dat)[i], dat[irow, i])
  # get the entropy of the elements using the chemical formula in thermo$obigt
  if(is.null(thisinfo)) thisinfo <- info(info(name, "cr_Berman", check.it=FALSE))
  SPrTr_elements <- convert(entropy(thisinfo$formula), "J")
  # check that G in data file is the G of formation from the elements --> Benson-Helgeson convention (DG = DH - T*DS)
  if(check.G) {
    GfPrTr_calc <- HfPrTr - Tr * (SPrTr - SPrTr_elements)
    Gdiff <- GfPrTr_calc - GfPrTr
    if(is.na(GfPrTr)) warning(paste0(name, ": GfPrTr(table) is NA"), call.=FALSE)
    else if(abs(Gdiff) >= 1000) warning(paste0(name, ": GfPrTr(calc) - GfPrTr(table) is too big! == ",
                                          round(GfPrTr_calc - GfPrTr), " J/mol"), call.=FALSE)
    # (the tabulated GfPrTr is unused below)
  }

  ### thermodynamic properties ###
  # calculate Cp and V (Berman, 1988 Eqs. 4 and 5)
  Cp <- k0 + k1 * T^-0.5 + k2 * T^-2 + k3 * T^-3
  P_Pr <- P - Pr
  T_Tr <- T - Tr
  V <- VPrTr * (1 + v1 * P_Pr + v2 * P_Pr^2 + v3 * T_Tr + v4 * T_Tr^2)
  # calculate Ga (Ber88 Eq. 6) --> Berman-Brown convention (DG = DH - T*S)
  Ga <- HfPrTr - T * SPrTr + k0 * ( (T - Tr) - T * (log(T) - log(Tr)) ) +
    2 * k1 * ( (T^0.5 - Tr^0.5) + T*(T^-0.5 - Tr^-0.5) ) -
    k2 * ( (T^-1 - Tr^-1) - T / 2 * (T^-2 - Tr^-2) ) -
    k3 * ( (T^-2 - Tr^-2) / 2 - T / 3 * (T^-3 - Tr^-3) ) +
    VPrTr * ( (v1 / 2 - v2) * (P^2 - Pr^2) + v2 / 3 * (P^3 - Pr^3) +
      (1 - v1 + v2 + v3 * (T - Tr) + v4 * (T - Tr)^2) * (P - Pr) )
  # calculate Ha (symbolically integrated using sympy - expressions not simplified)
  intCp <- T*k0 - Tr*k0 + k2/Tr - k2/T + k3/(2*Tr^2) - k3/(2*T^2) + 2.0*k1*T^0.5 - 2.0*k1*Tr^0.5
  intVminusTdVdT <- -VPrTr + P*(VPrTr + VPrTr*v2 - VPrTr*v1 - Tr*VPrTr*v3 + VPrTr*v4*Tr^2 - VPrTr*v4*T^2) +
    P^2*(VPrTr*v1/2 - VPrTr*v2) + VPrTr*v1/2 - VPrTr*v2/3 + Tr*VPrTr*v3 + VPrTr*v4*T^2 - VPrTr*v4*Tr^2 + VPrTr*v2*P^3/3
  Ha <- HfPrTr + intCp + intVminusTdVdT
  # calculate S (also symbolically integrated)
  intCpoverT <- k0*log(T) - k0*log(Tr) - k3/(3*T^3) + k3/(3*Tr^3) + k2/(2*Tr^2) - k2/(2*T^2) + 2.0*k1*Tr^-0.5 - 2.0*k1*T^-0.5
  intdVdT <- -VPrTr*(v3 + v4*(-2*Tr + 2*T)) + P*VPrTr*(v3 + v4*(-2*Tr + 2*T))
  S <- SPrTr + intCpoverT - intdVdT

  ### polymorphic transition properties ***
  if(!is.na(Tlambda) & !is.na(Tref) & any(T > Tref) & calc.transition) {
    # starting transition contributions are 0
    Cptr <- Htr <- Str <- Gtr <- numeric(ncond)
    ## Ber88 Eq. 8: Cp at 1 bar
    #Cplambda_1bar <- T * (l1 + l2 * T)^2
    # Eq. 9: Tlambda at P
    Tlambda_P <- Tlambda + dTdP * (P - 1)
    # Eq. 8a: Cp at P
    Td <- Tlambda - Tlambda_P
    Tprime <- T + Td
    # with the condition that Tref < Tprime < Tlambda(1bar)
    iTprime <- Tref < Tprime & Tprime < Tlambda
    Tprime <- Tprime[iTprime]
    Cptr[iTprime] <- Tprime * (l1 + l2 * Tprime)^2
    # we got Cp, now calculate the integrations for H and S
    # the lower integration limit is Tref
    iTtr <- T > Tref
    Ttr <- T[iTtr]
    Tlambda_P <- Tlambda_P[iTtr]
    Td <- Td[iTtr]
    # the upper integration limit is Tlambda_P
    Ttr[Ttr >= Tlambda_P] <- Tlambda_P[Ttr >= Tlambda_P]
    # derived variables
    tref <- Tref - Td
    x1 <- l1^2 * Td + 2 * l1 * l2 * Td^2 + l2^2 * Td^3
    x2 <- l1^2 + 4 * l1 * l2 * Td + 3 * l2^2 * Td^2
    x3 <- 2 * l1 * l2 + 3 * l2^2 * Td
    x4 <- l2 ^ 2
    # Eqs. 10, 11, 12
    Htr[iTtr] <- x1 * (Ttr - tref) + x2 / 2 * (Ttr^2 - tref^2) + x3 / 3 * (Ttr^3 - tref^3) + x4 / 4 * (Ttr^4 - tref^4)
    Str[iTtr] <- x1 * (log(Ttr) - log(tref)) + x2 * (Ttr - tref) + x3 / 2 * (Ttr^2 - tref^2) + x4 / 3 * (Ttr^3 - tref^3)
    Gtr <- Htr - T * Str
    # apply the transition contributions
    Ga <- Ga + Gtr
    Ha <- Ha + Htr
    S <- S + Str
    Cp <- Cp + Cptr
  }

  ### disorder thermodynamic properties ###
  if(!is.na(Tmin) & !is.na(Tmax) & any(T > Tmin) & calc.disorder) {
    # starting disorder contributions are 0
    Cpds <- Hds <- Sds <- Vds <- Gds <- 0
    # the lower integration limit is Tmin
    iTds <- T > Tmin
    Tds <- T[iTds]
    # the upper integration limit is Tmax
    Tds[Tds > Tmax] <- Tmax
    # Ber88 Eqs. 15, 16, 17
    Cpds[iTds] <- d0 + d1*Tds^-0.5 + d2*Tds^-2 + d3*Tds + d4*Tds^2
    Hds[iTds] <- d0*(Tds - Tmin) + d1*(Tds^0.5 - Tmin^0.5)/0.5 +
      d2*(Tds^-1 - Tmin^-1)/-1 + d3*(Tds^2 - Tmin^2)/2 + d4*(Tds^3 - Tmin^3)/3
    Sds[iTds] <- d0*(log(Tds) - log(Tmin)) + d1*(Tds^-0.5 - Tmin^-0.5)/-0.5 +
      d2*(Tds^-2 - Tmin^-2)/-2 + d3*(Tds - Tmin) + d4*(Tds^2 - Tmin^2)/2
    # Eq. 18; we can't do this if d5 == 0 (dolomite and gehlenite)
    # "d5 is a constant computed in such as way as to scale the disordring enthalpy to the volume of disorder" (Berman, 1988)
    if(d5 != 0) Vds <- Hds / d5
    # Berman puts the Vds term directly into Eq. 19 (commented below), but that necessarily makes Gds != Hds - T * Sds
    #Gds <- Hds - T * Sds + Vds * (P - Pr)
    # instead, we include the Vds term with Hds
    Hds <- Hds + Vds * (P - Pr)
    # disordering properties above Tmax (Eq. 20)
    ihigh <- T > Tmax
    # again, Berman put the Sds term (for T > Tmax) into Eq. 20 for Gds (commented below), which would also make Gds != Hds - T * Sds
    #Gds[ihigh] <- Gds[ihigh] - (T[ihigh] - Tmax) * Sds[ihigh]
    # instead, we add the Sds[ihigh] term to Hds
    Hds[ihigh] <- Hds[ihigh] - (T[ihigh] - Tmax) * Sds[ihigh]
    # by writing Gds = Hds - T * Sds, the above two changes w.r.t. Berman's
    # equations affect the computed values only for Hds, not Gds
    Gds <- Hds - T * Sds
    # apply the disorder contributions
    Ga <- Ga + Gds
    Ha <- Ha + Hds
    S <- S + Sds
    V <- V + Vds
    Cp <- Cp + Cpds
  }

  ### (for testing) use G = H - TS to check that integrals for H and S are written correctly
  Ga_fromHminusTS <- Ha - T * S
  if(!isTRUE(all.equal(Ga_fromHminusTS, Ga))) stop(paste0(name, ": incorrect integrals detected using DG = DH - T*S"))

  ### thermodynamic and unit conventions used in SUPCRT ###
  # use entropy of the elements in calculation of G --> Benson-Helgeson convention (DG = DH - T*DS)
  Gf <- Ga + Tr * SPrTr_elements
  # the output will just have "G" and "H"
  G <- Gf
  H <- Ha
  # convert J to cal
  if(grepl("cal", units)) {
    G <- convert(Gf, "cal")
    H <- convert(Ha, "cal")
    S <- convert(S, "cal")
    Cp <- convert(Cp, "cal")
  }
  # convert J/bar to cm^3/mol
  V <- V * 10

  data.frame(T=T, P=P, G=G, H=H, S=S, Cp=Cp, V=V)
}
