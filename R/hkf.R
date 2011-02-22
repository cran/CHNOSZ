# CHNOSZ/hkf.R
# calculate thermodynamic properties using equations of state
# 11/17/03 jmd

hkf <- function(property=NULL,T=298.15,P=1,ghs=NULL,eos=NULL,contrib=c('n','s','o'),
  H2O.PT=NULL,H2O.PrTr=NULL,domega=TRUE) {
  # calculate G, H, S, Cp, V, kT, and/or E using
  # the revised HKF equations of state
  # constants
  Tr <- thermo$opt$Tr
  Pr <- thermo$opt$Pr
  Theta <- thermo$opt$Theta
  Psi <- thermo$opt$Psi
  # argument handling
  eargs <- eos.args('hkf',property)
  property <- eargs$prop
  props <- eargs$props
  Prop <- eargs$Prop
  domega <- rep(domega,length.out=nrow(ghs))
  # nonsolvation, solvation, and origination contribution
  contribs <- c('n','s','o')
  notcontrib <- ! contrib %in% contribs
  if(TRUE %in% notcontrib)
    stop(paste('argument',c2s(contrib[notcontrib]),'not in',c2s(contribs),'n'))
  # get water properties, if they weren't supplied in arguments (and we want solvation props)
  if('s' %in% contrib) {
    H2O.props <- c('Q','X','Y','epsilon')
    # only take these ones if we're in SUPCRT92 compatibility mode
    dosupcrt <- length(agrep(tolower(thermo$opt$water),'supcrt9',max.distance=0.3))!=0
    if(dosupcrt) {
      # (E, daldT, V - for partial derivatives of omega (g function))
      H2O.props <- c(H2O.props,'E','daldT','kT','Z')
    } else {
      # (N, UBorn - for compressibility, expansibility)
      H2O.props <- c(H2O.props,'N','UBorn')
    }
    if(is.null(H2O.PT)) H2O.PT <- water(H2O.props,T=T,P=P)
    if(is.null(H2O.PrTr)) H2O.PrTr <- water(H2O.props,T=thermo$opt$Tr,P=thermo$opt$Pr)
    ZBorn <- -1/H2O.PT$epsilon
    ZBorn.PrTr <- -1/H2O.PrTr$epsilon
  }
 # a list to store the result
 x <- list()
 for(k in 1:nrow(ghs)) {
  # loop over each species
  GHS <- ghs[k,]
  EOS <- eos[k,]
  # compute values of omega(P,T) from those of omega(Pr,Tr)
  # using g function etc. (Shock et al., 1992 and others)
  omega <- EOS$omega  # omega.PrTr
  Z <- EOS$Z
  omega.PT <- rep(EOS$omega,length(T))
  d2wdT2 <- dwdT <- dwdP <- rep(0,length(T))
  rhohat <- H2O.PT$rho/1000  # just converting kg/m3 to g/cm3
  do.g <- rhohat < 1
  # temporarily turn off warnings
  warn <- options(warn=-1)
  if(!is.na(Z)) if(Z != 0) if(domega[k]) if(dosupcrt) if(any(do.g)) {
    # g and f function stuff (Shock et al., 1992; Johnson et al., 1992)
    eta <- 1.66027E5
    ag1 <- -2.037662; ag2 <- 5.747000E-3; ag3 <- -6.557892E-6
    bg1 <- 6.107361; bg2 <- -1.074377E-2; bg3 <- 1.268348E-5
    af1 <- 0.3666666E2
    af2 <- -0.1504956E-9
    af3 <- 0.5017997E-13
    ag <- function(T) ag1 + ag2 * T + ag3 * T ^ 2
    bg <- function(T) bg1 + bg2 * T + bg3 * T ^ 2
    g <- function(rho,T,P,ag,bg,do.f) {
      g <- ag*(1-rho)^bg
      f <- ( ((T-155)/300)^4.8 + af1*((T-155)/300)^16 ) *
           ( af2*(1000-P)^3 + af3*(1000-P)^4 )
      g[do.f] <- g[do.f] - f[do.f]
      # why is this necessary - 
      # should we warn if this happens?
      g[is.nan(g)] <- 0
      return(g)
    }
    # on to the calculation
    Tc <- convert(T,'C')
    do.f <- !(Tc < 155 | P > 1000 | Tc > 355)
    # replace the functions with their numerical results
    ag <- ag(Tc)
    bg <- bg(Tc)
    g <- g(rhohat,Tc,P,ag,bg,do.f)
    # partial derivatives of omega - equation numbers from JOH92
    alpha <- H2O.PT$E  # expansivity - dVdT
    daldT <- H2O.PT$daldT  # d2VdT2
    beta <- H2O.PT$kT  # isothermal compressibility
    # Eqn. 76
    d2fdT2 <- (0.0608/300*((Tc-155)/300)^2.8 + af1/375*((Tc-155)/300)^14) * (af2*(1000-P)^3 + af3*(1000-P)^4)
    # Eqn. 75
    dfdT <- (0.016*((Tc-155)/300)^3.8 + 16*af1/300*((Tc-155)/300)^15) * (af2*(1000-P)^3 + af3*(1000-P)^4)
    # Eqn. 74
    dfdP <- -(((Tc-155)/300)^4.8 + af1*((Tc-155)/300)^16) * (3*af2*(1000-P)^2 + 4*af3*(1000-P)^3)
    d2bdT2 <- 2 * bg3  # Eqn. 73
    d2adT2 <- 2 * ag3  # Eqn. 72
    dbdT <- bg2 + 2*bg3*Tc  # Eqn. 71
    dadT <- ag2 + 2*ag3*Tc  # Eqn. 70
    # Eqn. 69
    dgadT <- bg*rhohat*alpha*(1-rhohat)^(bg-1) + log(1-rhohat)*g/ag*dbdT  
    D <- rhohat
    # transcribed from SUPCRT92/reac92.f
    dDdT <- -D * alpha
    dDdP <- D * beta
    dDdTT <- -D * (daldT - alpha^2)
    Db <- (1-D)^bg
    dDbdT <- -bg*(1-D)^(bg-1)*dDdT + log(1-D)*Db*dbdT
    dDbdTT <- -(bg*(1-D)^(bg-1)*dDdTT + (1-D)^(bg-1)*dDdT*dbdT + 
      bg*dDdT*(-(bg-1)*(1-D)^(bg-2)*dDdT + log(1-D)*(1-D)^(bg-1)*dbdT)) +
      log(1-D)*(1-D)^bg*d2bdT2 - (1-D)^bg*dbdT*dDdT/(1-D) + log(1-D)*dbdT*dDbdT
    d2gdT2 <- ag*dDbdTT + 2*dDbdT*dadT + Db*d2adT2
    d2gdT2[do.f] <- d2gdT2[do.f] - d2fdT2[do.f]
    dgdT <- g/ag*dadT + ag*dgadT  # Eqn. 67
    dgdT[do.f] <- dgdT[do.f] - dfdT[do.f]
    dgdP <- -bg*rhohat*beta*g*(1-rhohat)^-1  # Eqn. 66
    dgdP[do.f] <- dgdP[do.f] - dfdP[do.f]
    # after reac92.f
    reref <- Z^2 / (omega/eta + Z/(3.082 + 0))
    re <- reref + abs(Z) * g
    omega.PT <- eta * (Z^2/re - Z/(3.082 + g))
    Z3 <- abs(Z^3)/re^2 - Z/(3.082 + g)^2
    Z4 <- abs(Z^4)/re^3 - Z/(3.082 + g)^3
    dwdP[do.g] <- (-eta * Z3 * dgdP)[do.g]
    dwdT[do.g] <- (-eta * Z3 * dgdT)[do.g]
    d2wdT2[do.g] <- (2 * eta * Z4 * dgdT^2 - eta * Z3 * d2gdT2)[do.g]
  }

  # re-enable warnings
  options(warn)
  # loop over each property
  w <- NULL
  for(i in 1:length(property)) {
    prop <- property[i]
    # over nonsolvation, solvation, or origination contributions
    hkf.p <- numeric(length(T))
    for(icontrib in contrib) {
      # various contributions to the properties
      if( icontrib=="n") {
        # nonsolvation ghs equations
        if(prop=="h") {
          p.c <- EOS$c1*(T-Tr) - EOS$c2*(1/(T-Theta)-1/(Tr-Theta))
          p.a <- EOS$a1*(P-Pr) + EOS$a2*log((Psi+P)/(Psi+Pr)) + ((2*T-Theta)/(T-Theta)^2)*(EOS$a3*(P-Pr)+EOS$a4*log((Psi+P)/(Psi+Pr)))
          p <- p.c + p.a
        }
        if(prop=="s") {
          p.c <- EOS$c1*log(T/Tr) - (EOS$c2/Theta)*( 1/(T-Theta)-1/(Tr-Theta) + log( (Tr*(T-Theta))/(T*(Tr-Theta)) )/Theta )
          p.a <- (T-Theta)^(-2)*(EOS$a3*(P-Pr)+EOS$a4*log((Psi+P)/(Psi+Pr)))
          p <- p.c + p.a
        }
        if(prop=="g") {
          p.c <- -EOS$c1*(T*log(T/Tr)-T+Tr) - EOS$c2*( (1/(T-Theta)-1/(Tr-Theta))*((Theta-T)/Theta) - (T/Theta^2)*log((Tr*(T-Theta))/(T*(Tr-Theta))) )
          p.a <- EOS$a1*(P-Pr) + EOS$a2*log((Psi+P)/(Psi+Pr)) + (EOS$a3*(P-Pr) + EOS$a4*log((Psi+P)/(Psi+Pr)))/(T-Theta)
          p <- p.c + p.a
        }
        # nonsolvation cp v kt e equations
        if(prop=='cp') p <- EOS$c1 + EOS$c2 * ( T - Theta ) ^ (-2)        
        if(prop=='v') p <- convert(EOS$a1,'cm3bar') + convert(EOS$a2,'cm3bar') / ( Psi + P) +
          ( convert(EOS$a3,'cm3bar') + convert(EOS$a4,'cm3bar') / ( Psi + P ) ) / ( T - Theta)
        if(prop=='kt') p <- ( convert(EOS$a2,'cm3bar') + convert(EOS$a4,'cm3bar') / (T - Theta) ) * (Psi + P) ^ (-2)
        if(prop=='e') p <- convert( - ( EOS$a3 + EOS$a4 / convert((Psi + P),'calories') ) * (T - Theta) ^ (-2),'cm3bar')
      }
      if( icontrib=="s") {
        # solvation ghs equations
        if(prop=="g") 
          p <- -omega.PT*(ZBorn+1) + omega*(ZBorn.PrTr+1) + omega*H2O.PrTr$Y*(T-Tr)
        if(prop=="h") 
          p <- -omega.PT*(ZBorn+1) + omega.PT*T*H2O.PT$Y + T*(ZBorn+1)*dwdT +
                 omega*(ZBorn.PrTr+1) - omega*Tr*H2O.PrTr$Y
        if(prop=="s") 
          p <- omega.PT*H2O.PT$Y + (ZBorn+1)*dwdT - omega*H2O.PrTr$Y 
        # solvation cp v kt e equations
        if(prop=='cp') p <- omega.PT*T*H2O.PT$X + 2*T*H2O.PT$Y*dwdT + T*(ZBorn+1)*d2wdT2
        if(prop=='v') p <- -convert(omega.PT,'cm3bar') * H2O.PT$Q + convert(dwdP,'cm3bar') * (-ZBorn - 1)
        # WARNING: the partial derivatives of omega are not implemented for kt and e
        # (to do it, see p. 820 of SOJ+92)
        if(prop=='kt') p <- convert(omega,'cm3bar') * H2O.PT$N
        if(prop=='e') p <- -convert(omega,'cm3bar') * H2O.PT$UBorn
      }
      if( icontrib=='o') {
        # origination ghs equations
        if(prop=='g') p <- GHS$G - GHS$S * (T-Tr)
        else if(prop=='h') p <- GHS$H
        else if(prop=='s') p <- GHS$S
        # origination eos equations: senseless
        else if(prop %in% tolower(props)) p <- 0 * T
      }
      p <- rep(p,length.out=length(hkf.p))
      ip <- 1:length(p)
      if(prop %in% c('g','h','s')) ip <- which(!is.na(p))
      hkf.p[ip] <- hkf.p[ip] + p[ip]
    }
    wnew <- data.frame(hkf.p)
    if(i>1) w <- cbind(w,wnew) else w <- wnew
  }
  colnames(w) <- Prop
  x[[k]] <- w
 }
 return(x)
}

