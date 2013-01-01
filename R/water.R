# CHNOSZ/water.R
# calculate thermodynamic and electrostatic properties of H2O
# 20061016 jmd

water.AW90 <- function(T=298.15,rho=1000,P=0.1) {
  # Equations for the dielectric constant of water
  # from Archer and Wang, 1990
  # T in K
  # rho in kg m-3
  # p in MPa

  # Table 2
  b <- c(-4.044525E-2, 103.6180   , 75.32165   ,
         -23.23778   ,-3.548184   ,-1246.311   ,
         263307.7    ,-6.928953E-1,-204.4473)
  alpha <- 18.1458392E-30 # m^3
  #alpha <- 14.7E-30
  mu <- 6.1375776E-30 # C m
  N.A <- 6.0221367E23 # mol-1
  k <- 1.380658E-23 # Boltzmann constant, J K-1
  M <- 0.0180153 # kg mol-1
  rho.0 <- 1000 # kg m-3
  # Equation 1
  epsilon.0 <- 8.8541878E-12 # permittivity of vacuum, C^2 J-1 m-1
  epsfun.lhs <- function(e) (e-1)*(2*e+1)/(9*e)
  epsfun.rhs <- function(T,V.m) N.A*(alpha+mufun()/(3*epsilon.0*k*T))/(3*V.m)
  epsfun <- function(e,T,V.m) epsfun.lhs(e) - epsfun.rhs(T,V.m)
  mufun <- function() gfun()*mu^2
  gfun <- function() rhofun()*rho/rho.0 + 1
  # Equation 3
  rhofun <- function() b[1]*P*T^-1 + b[2]*T^-0.5 + b[3]*(T-215)^-1 +
    b[4]*(T-215)^-0.5 + b[5]*(T-215)^-0.25 +
    exp(b[6]*T^-1 + b[7]*T^-2 + b[8]*P*T^-1 + b[9]*P*T^-2)
  epsilon <- function(T,rho) {
    tu <- try(uniroot(epsfun,c(1E-1,1E3),T=T,V.m=M/rho)$root,TRUE)
    if(!is.numeric(tu)) {
      warning('water.AW90: no root for density at ',T,' K and ',rho,' kg m-3.',call.=FALSE,immediate.=TRUE)
      tu <- NA
    }
    return(tu)
  }
  # get things the right length
  our.T <- T; our.rho <- rho; our.P <- P
  t <- numeric()
  for(i in 1:length(our.T)) {
    T <- our.T[i]
    rho <- our.rho[i]
    P <- our.P[i]
    t <- c(t,epsilon(T,rho))
  }
  return(t)
}

idealgas.IAPWS95 <- function(p, delta, tau) {
  ## the ideal gas part in the IAPWS-95 formulation
  # from Table 6.1 of Wagner and Pruss, 2002
  n <- c( -8.32044648201, 6.6832105268, 3.00632, 0.012436,
           0.97315, 1.27950, 0.96956, 0.24873 )
  gamma <- c( NA, NA, NA, 1.28728967, 
              3.53734222, 7.74073708, 9.24437796, 27.5075105 )
  # Equation 6.5
  phi <- function() log(delta) + n[1] + n[2]*tau + n[3]*log(tau) +
    sum( n[4:8] * log(1-exp(-gamma[4:8]*tau)) )
  # derivatives from Table 6.4
  phi.delta <- function() 1/delta+0+0+0+0
  phi.delta.delta <- function() -1/delta^2+0+0+0+0
  phi.tau <- function() 0+0+n[2]+n[3]/tau+sum(n[4:8]*gamma[4:8]*((1-exp(-gamma[4:8]*tau))^-1-1))
  phi.tau.tau <- function() 0+0+0-n[3]/tau^2-sum(n[4:8]*gamma[4:8]^2 * 
    exp(-gamma[4:8]*tau)*(1-exp(-gamma[4:8]*tau))^-2)
  phi.delta.tau <- function() 0+0+0+0+0
  return(get(p)())
}

residual.IAPWS95 <- function(p, delta, tau) {
  ## the residual part in the IAPWS-95 formulation
  # from Table 6.2 of Wagner and Pruss, 2002
  c <- c(rep(NA,7),rep(1,15),rep(2,20),rep(3,4),4,rep(6,4),rep(NA,5))
  d <- c(1,1,1,2,2,3,4,1,1,1,2,2,3,4,
         4,5,7,9,10,11,13,15,1,2,2,2,3,4,
         4,4,5,6,6,7,9,9,9,9,9,10,10,12,
         3,4,4,5,14,3,6,6,6,3,3,3,NA,NA)
  t <- c(-0.5,0.875,1,0.5,0.75,0.375,1,4,6,12,1,5,4,2,
         13,9,3,4,11,4,13,1,7,1,9,10,10,3,
         7,10,10,6,10,10,1,2,3,4,8,6,9,8,
         16,22,23,23,10,50,44,46,50,0,1,4,NA,NA)
  n <- c( 0.12533547935523E-1, 0.78957634722828E1 ,-0.87803203303561E1 ,
          0.31802509345418   ,-0.26145533859358   ,-0.78199751687981E-2,
          0.88089493102134E-2,-0.66856572307965   , 0.20433810950965   ,
         -0.66212605039687E-4,-0.19232721156002   ,-0.25709043003438   ,
          0.16074868486251   ,-0.40092828925807E-1, 0.39343422603254E-6,
         -0.75941377088144E-5, 0.56250979351888E-3,-0.15608652257135E-4,
          0.11537996422951E-8, 0.36582165144204E-6,-0.13251180074668E-11,
         -0.62639586912454E-9,-0.10793600908932   , 0.17611491008752E-1,
          0.22132295167546   ,-0.40247669763528   , 0.58083399985759   ,
          0.49969146990806E-2,-0.31358700712549E-1,-0.74315929710341   ,
          0.47807329915480   , 0.20527940895948E-1,-0.13636435110343   ,
          0.14180634400617E-1, 0.83326504880713E-2,-0.29052336009585E-1,
          0.38615085574206E-1,-0.20393486513704E-1,-0.16554050063734E-2,
          0.19955571979541E-2, 0.15870308324157E-3,-0.16388568342530E-4,
          0.43613615723811E-1, 0.34994005463765E-1,-0.76788197844621E-1,
          0.22446277332006E-1,-0.62689710414685E-4,-0.55711118565645E-9,
         -0.19905718354408   , 0.31777497330738   ,-0.11841182425981   ,
         -0.31306260323435E2 , 0.31546140237781E2 ,-0.25213154341695E4 ,
         -0.14874640856724   , 0.31806110878444)
  alpha <- c(rep(NA,51),20,20,20,NA,NA)
  beta <- c(rep(NA,51),150,150,250,0.3,0.3)
  gamma <- c(rep(NA,51),1.21,1.21,1.25,NA,NA)
  epsilon <- c(rep(NA,51),1,1,1,NA,NA)
  a <- c(rep(NA,54),3.5,3.5)
  b <- c(rep(NA,54),0.85,0.95)
  B <- c(rep(NA,54),0.2,0.2)
  C <- c(rep(NA,54),28,32)
  D <- c(rep(NA,54),700,800)
  A <- c(rep(NA,54),0.32,0.32)
  # from Table 6.5
  i1 <- 1:7
  i2 <- 8:51
  i3 <- 52:54
  i4 <- 55:56
  # deriviatives of distance function
  Delta <- function(i) { Theta(i)^2 + B[i] * ((delta-1)^2)^a[i] }
  Theta <- function(i) { (1-tau) + A[i] * ((delta-1)^2)^(1/(2*beta[i])) }
  Psi <- function(i) { exp ( -C[i]*(delta-1)^2 - D[i]*(tau-1)^2 ) }
  dDelta.bi.ddelta <- function(i) { b[i]*Delta(i)^(b[i]-1)*dDelta.ddelta(i) }
  d2Delta.bi.ddelta2 <- function(i) { b[i]*( Delta(i)^(b[i]-1) * d2Delta.ddelta2(i) + 
    (b[i]-1)*Delta(i)^(b[i]-2)*dDelta.ddelta(i)^2 ) }
  dDelta.bi.dtau <- function(i) { -2*Theta(i)*b[i]*Delta(i)^(b[i]-1) }
  d2Delta.bi.dtau2 <- function(i) { 2*b[i]*Delta(i)^(b[i]-1) + 4*Theta(i)^2*b[i]*(b[i]-1)*Delta(i)^(b[i]-2) }
  d2Delta.bi.ddelta.dtau <- function(i) { -A[i]*b[i]*2/beta[i]*Delta(i)^(b[i]-1)*(delta-1) * 
    ((delta-1)^2)^(1/(2*beta[i])-1) - 2*Theta(i)*b[i]*(b[i]-1)*Delta(i)^(b[i]-2)*dDelta.ddelta(i)  }
  dDelta.ddelta <- function(i) { (delta-1) * ( A[i]*Theta(i)*2/beta[i]*((delta-1)^2)^(1/(2*beta[i])-1) +
    2*B[i]*a[i]*((delta-1)^2)^(a[i]-1) ) }
  d2Delta.ddelta2 <- function(i) { 1/(delta-1)*dDelta.ddelta(i) + (delta-1)^2 * (
    4*B[i]*a[i]*(a[i]-1)*((delta-1)^2)^(a[i]-2) + 2*A[i]^2*(1/beta[i])^2 *
      (((delta-1)^2)^(1/(2*B[i])-1))^2 + A[i]*Theta(i)*4/beta[i]*(1/(2*B[i])-1) *
        ((delta-1)^2)^(1/(2*beta[i])-2) ) }
  # derivatives of exponential function
  dPsi.ddelta <- function(i) { -2*C[i]*(delta-1)*Psi(i) }
  d2Psi.ddelta2 <- function(i) { ( 2*C[i]*(delta-1)^2 - 1 ) * 2*C[i]*Psi(i) }
  dPsi.dtau <- function(i) { -2*D[i]*(tau-1)*Psi(i) }
  d2Psi.dtau2 <- function(i) { (2*D[i]*(tau-1)^2 - 1) * 2*D[i]*Psi(i) }
  d2Psi.ddelta.dtau <- function(i) { 4*C[i]*D[i]*(delta-1)*(tau-1)*Psi(i) }
  # dimensionless Helmholtz free energy and derivatives
  phi <- function() {
    sum(n[i1]*delta^d[i1]*tau^t[i1]) +
    sum(n[i2]*delta^d[i2]*tau^t[i2]*exp(-delta^c[i2])) +
    sum(n[i3]*delta^d[i3]*tau^t[i3] *
      exp( -alpha[i3]*(delta-epsilon[i3])^2 - beta[i3]*(tau-gamma[i3])^2 ) ) +
    sum(n[i4]*Delta(i4)^b[i4]*delta*Psi(i4))
  }
  phi.delta <- function() {
    sum(n[i1]*d[i1]*delta^(d[i1]-1)*tau^t[i1]) +
    sum(n[i2]*exp(-delta^c[i2])*(delta^(d[i2]-1)*tau^t[i2]*(d[i2]-c[i2]*delta^c[i2]))) +
    sum(n[i3]*delta^d[i3]*tau^t[i3] *
      exp( -alpha[i3]*(delta-epsilon[i3])^2 - beta[i3]*(tau-gamma[i3])^2 ) * 
        (d[i3]/delta - 2 * alpha[i3]*(delta-epsilon[i3])) ) +
    sum(n[i4] * ( Delta(i4)^b[i4] * (Psi(i4)+delta*dPsi.ddelta(i4)) + dDelta.bi.ddelta(i4)*delta*Psi(i4) ) ) 
  }
  phi.delta.delta <- function() {
    sum(n[i1]*d[i1]*(d[i1]-1)*delta^(d[i1]-2)*tau^t[i1]) +
    sum(n[i2]*exp(-delta^c[i2])*(delta^(d[i2]-2)*tau^t[i2]*((d[i2]-c[i2]*delta^c[i2]) * 
      (d[i2]-1-c[i2]*delta^c[i2])-c[i2]^2*delta^c[i2]))) +
    sum(n[i3]*tau^t[i3]*exp(-alpha[i3]*(delta-epsilon[i3])^2 - beta[i3]*(tau-gamma[i3])^2) * (
      -2*alpha[i3]*delta^d[i3]+4*alpha[i3]^2*delta^d[i3]*(delta-epsilon[i3])^2 -
       4*d[i3]*alpha[i3]*delta^(d[i3]-1)*(delta-epsilon[i3])+d[i3]*(d[i3]-1)*delta^(d[i3]-2) ) ) +
    sum(n[i4]*( Delta(i4)^b[i4]*(2*dPsi.ddelta(i4)+delta*d2Psi.ddelta2(i4)) + 
      2*dDelta.bi.ddelta(i4)*(Psi(i4)+delta*dPsi.ddelta(i4)) + d2Delta.bi.ddelta2(i4)*delta*Psi(i4) ) )
  }
  phi.tau <- function() {
    sum(n[i1]*t[i1]*delta^d[i1]*tau^(t[i1]-1)) +
    sum(n[i2]*t[i2]*delta^d[i2]*tau^(t[i2]-1)*exp(-delta^c[i2])) +
    sum(n[i3]*delta^d[i3]*tau^t[i3]*exp(-alpha[i3]*(delta-epsilon[i3])^2-beta[i3]*(tau-gamma[i3])^2) * 
      (t[i3]/tau-2*beta[i3]*(tau-gamma[i3]))) +
    sum(n[i4]*delta*(dDelta.bi.dtau(i4)*Psi(i4)+Delta(i4)^b[i4]*dPsi.dtau(i4)))
  }
  phi.tau.tau <- function() {
    sum(n[i1]*t[i1]*(t[i1]-1)*delta^d[i1]*tau^(t[i1]-2)) +
    sum(n[i2]*t[i2]*(t[i2]-1)*delta^d[i2]*tau^(t[i2]-2)*exp(-delta^c[i2])) +
    sum(n[i3]*delta^d[i3]*tau^t[i3]*exp(-alpha[i3]*(delta-epsilon[i3])^2-beta[i3]*(tau-gamma[i3])^2) * 
      (((t[i3]/tau)-2*beta[i3]*(tau-gamma[i3]))^2-t[i3]/tau^2-2*beta[i3])) +
    sum(n[i4]*delta*(d2Delta.bi.dtau2(i4)*Psi(i4)+2*dDelta.bi.dtau(i4)*dPsi.dtau(i4) + 
      Delta(i4)^b[i4]*d2Psi.dtau2(i4)))
  }
  phi.delta.tau <- function() {
    sum(n[i1]*d[i1]*t[i1]*delta^(d[i1]-1)*tau^(t[i1]-1)) +
    sum(n[i2]*t[i2]*delta^(d[i2]-1)*tau^(t[i2]-1)*(d[i2]-c[i2]*delta^c[i2])*exp(-delta^c[i2])) +
    sum(n[i3]*delta^d[i3]*tau^t[i3]*exp(-alpha[i3]*(delta-epsilon[i3])^2-beta[i3]*(tau-gamma[i3])^2) * 
      ((d[i3]/delta)-2*alpha[i3]*(delta-epsilon[i3]))*(t[i3]/tau-2*beta[i3]*(tau-gamma[i3])) ) +
    sum(n[i4]*(Delta(i4)^b[i4]*(dPsi.dtau(i4)+delta*d2Psi.ddelta.dtau(i4)) + 
      delta*dDelta.bi.ddelta(i4)*dPsi.dtau(i4)+dDelta.bi.dtau(i4) * (Psi(i4)+delta*dPsi.ddelta(i4)) +
      d2Delta.bi.ddelta.dtau(i4)*delta*Psi(i4) )) 
  }
  return(get(p)())
}

water.IAPWS95 <- function(property,T=298.15,rho=1000) {
  ## the IAPWS-95 formulation for pure H2O from 
  ## Wagner and Pruss, 2002
  property <- tolower(property)
  # triple point
  T.triple <- 273.16 # K
  P.triple <- 611.657 # Pa
  rho.triple.liquid <- 999.793
  rho.triple.vapor <- 0.00485458
  # normal boiling point
  T.boiling <- 373.124
  P.boiling <- 0.101325
  rho.boiling.liquid <- 958.367
  rho.boiling.vapor <- 0.597657
  # critical point constants
  T.critical <- 647.096 # K
  rho.critical <- 322 # kg m-3
  # specific and molar gas constants
  R <- 0.46151805 # kJ kg-1 K-1
  # R.M <- 8.314472 # J mol-1 K-1
  # molar mass
  M <- 18.015268 # g mol-1
  ## define functions idealgas and residual, supplying arguments delta and tau
  idealgas <- function(p) idealgas.IAPWS95(p, delta, tau)
  residual <- function(p) residual.IAPWS95(p, delta, tau)
  ## relation of thermodynamic properties to Helmholtz free energy
  a <- function() {
    x <- idealgas('phi')+residual('phi')
    return(x*R*T)
  }
  # Table 6.3 
  p <- function() {
    x <- 1 + delta*residual('phi.delta')
    return(x*rho*R*T/1000)  # for MPa
  }
  s <- function() {
    x <- tau * (idealgas('phi.tau')+residual('phi.tau'))-idealgas('phi')-residual('phi')
    return(x*R)
  }
  u <- function() {
    x <- tau * (idealgas('phi.tau')+residual('phi.tau'))
    return(x*R*T)
  }
  h <- function() {
    x <- 1 + tau * (idealgas('phi.tau')+residual('phi.tau')) + delta*residual('phi.delta')
    return(x*R*T)
  }
  g <- function() {
    x <- 1 + idealgas('phi') + residual('phi') + delta*residual('phi.delta')
    return(x*R*T)
  }
  cv <- function() {
    x <- -tau^2*(idealgas('phi.tau.tau')+residual('phi.tau.tau'))
    return(x*R)
  }
  cp <- function() {
    x <- -tau^2*(idealgas('phi.tau.tau')+residual('phi.tau.tau')) +
         (1+delta*residual('phi.delta')-delta*tau*residual('phi.delta.tau'))^2 /
         (1+2*delta*residual('phi.delta')+delta^2*residual('phi.delta.delta'))
    return(x*R)
  }
# 20090420 speed of sound calculation is incomplete
# (delta.liquid and drhos.dT not visible)
#  cs <- function() {
#    x <- -tau^2*(idealgas('phi.tau.tau')+residual('phi.tau.tau')) +
#         (1+delta*residual('phi.delta')-delta*tau*residual('phi.delta.tau'))^2 /
#         (1+2*delta*residual('phi.delta')+delta^2*residual('phi.delta.delta')) *
#         ((1+delta.liquid*residual('phi.delta')-delta.liquid*tau*residual('phi.tau.tau'))-rho.critical/(R*delta.liquid)*drhos.dT)
#    return(x*R)
#  }
  w <- function() {
    x <- 1 + 2*delta*residual('phi.delta') + delta^2*residual('phi.delta.delta') - 
         (1+delta*residual('phi.delta')-delta*tau*residual('phi.delta.tau'))^2 /
         tau^2*(idealgas('phi.tau.tau')+residual('phi.tau.tau'))
    return(sqrt(x*R*T))
  }
  mu <- function() {
    x <- -(delta*residual('phi.delta')+delta^2*residual('phi.delta.delta')+delta*tau*residual('phi.delta.tau')) /
          ( ( 1+delta*residual('phi.delta')-delta*tau*residual('phi.delta.tau')^2 ) - tau^2 *
          (idealgas('phi.tau.tau')+residual('phi.tau.tau'))*(1+2*delta*residual('phi.delta')+delta^2*residual('phi.delta.delta')) ) 
    return(x/(R*rho))
  }
  ## run the calculations
  ww <- NULL
  my.T <- T
  my.rho <- rho
  for(j in 1:length(property)) {
    t <- numeric()
    for(i in 1:length(my.T)) {
      T <- my.T[i]
      rho <- my.rho[i]
      # Equation 6.4
      delta <- rho / rho.critical
      tau <- T.critical / T
      t <- c(t,get(property[j])())
    }
    t <- data.frame(t)
    if(j==1) ww <- t else ww <- cbind(ww,t)
  }
  colnames(ww) <- property
  return(ww)
}

water.WP02 <- function(property='rho.liquid',T=298.15) {
  # auxiliary equations for liquid-vapor phase boundary
  # from Wagner and Pruss, 2002
  # critical point
  T.critical <- 647.096 # K
  P.critical <- 22.064 # MPa
  rho.critical <- 322 # kg m-3

  if(property %in% c("P.sigma","dP.sigma.dT")) {
    # vapor pressure
    V <- 1 - T / T.critical # theta (dimensionless)
    a1 <- -7.85951783
    a2 <- 1.84408259
    a3 <- -11.7866497
    a4 <- 22.6807411
    a5 <- -15.9618719
    a6 <- 1.80122502
    ln.P.sigma.P.critical <- T.critical / T *
      ( a1*V + a2*V^1.5 + a3*V^3 + a4*V^3.5 + a5*V^4 + a6*V^7.5 ) 
    P.sigma <- P.critical * exp(ln.P.sigma.P.critical)
    if(property=="dP.sigma.dT") out <- - P.sigma / T * ( ln.P.sigma.P.critical +
      a1 + 1.5*a2*V^0.5 + 3*a3*V^2 + 3.5*a4*V^2.5 + 4*a5*V^3 + 7.5*a6*V^6.5 )
    else out <- P.sigma
  } else if(property=="rho.liquid") {
    # saturated liquid density
    V <- 1 - T / T.critical
    b1 <- 1.99274064
    b2 <- 1.09965342
    b3 <- -0.510839303
    b4 <- -1.75493479
    b5 <- -45.5170352
    b6 <- -6.74694450E5
    rho.liquid <- rho.critical * ( 
      1 + b1*V^(1/3) + b2*V^(2/3) + b3*V^(5/3) + b4*V^(16/3) + b5*V^(43/3) + b6*V^(110/3) )
    out <- rho.liquid
  } else if(property=="rho.vapor") {
  # saturated vapor density
    V <- 1 - T / T.critical
    c1 <- -2.03150240
    c2 <- -2.68302940
    c3 <- -5.38626492
    c4 <- -17.2991605
    c5 <- -44.7586581
    c6 <- -63.9201063
    rho.vapor <- rho.critical * exp (
      c1*V^(2/6) + c2*V^(4/6) + c3*V^(8/6) + c4*V^(18/6) + c5*V^(37/6) + c6*V^(71/6) )
    out <- rho.vapor
  } else stop(paste('i can not calculate',property))
  return(out)
}

water <- function(property = NULL,T = thermo$opt$Tr, P = 'Psat') {
  # calculate the properties of liquid H2O as a function of T and P
  # T in Kelvin, P in bar
  if(is.null(property)) stop('property was NULL')
  # this tells us to do the calculations using code taken from SUPCRT
  do.supcrt <- length(agrep(tolower(thermo$opt$water),'supcrt9',max.distance=0.3)) > 0
  eargs <- eos.args('water',property=property,T=T,P=P)
  property <- eargs$prop; Property <- eargs$Prop
  # working out the arguments
  tpargs <- TP.args(T=T,P=P)
  P <- tpargs$P; T <- tpargs$T
  M <- 18.015268 # g mol-1
  # Psat stuff
  psat <- function() {
    p <- numeric()
    if(do.supcrt) {
      p <- water.SUPCRT92('',T=T,rep(0,length(T)),isat=1)
      p[p==0] <- NaN
      return(p)
    } else {
      for(i in 1:length(T)) {
        if(T[i] < 373.124) p <- c(p,0.1)
        else p <- c(p,water.WP02('P.sigma',T[i]))
      }
      return(convert(p,'bar'))
    }
  }
  # a quick return if property = 'Psat', for use by the TP.args() function
  if(length(property)==1 & property[1]=='psat') return(data.frame(Psat=psat()))
  ### maybe we are using the SUPCRT calculations ###
  if(do.supcrt) {
    names.SUPCRT <- c('Speed','alpha','beta','alpha','beta','diel','ZBorn','YBorn','QBorn','XBorn')
    names.CHNOSZ <- c('w','alpha','beta','E','kT','epsilon','Z','Y','Q','X')
    Property.new <- character()
    # convert names to SUPCRT
    for(i in 1:length(Property)) if(Property[i] %in% names.CHNOSZ) 
      Property.new[i] <- names.SUPCRT[match(Property[i],names.CHNOSZ)]
      else Property.new[i] <- Property[i]
    # deal with compressibility and expansivity 20091203
    iE <- which(Property=="E")
    ikT <- which(Property=="kT")
    iV <- numeric()
    if("kT" %in% Property | "E" %in% Property) iV <- length(Property.new <- c(Property.new,"V"))
    # get the value of the property
    w.out <- water.SUPCRT92(Property.new,T=T,P=P)
    # finish dealing with compressibility and expansivity
    if("E" %in% Property) w.out[,iE] <- w.out$V*w.out$alpha
    if("kT" %in% Property) w.out[,ikT] <- w.out$V*w.out$beta
    if(length(iV) > 0) w.out <- w.out[,-iV,drop=FALSE]
    colnames(w.out) <- Property
    return(w.out)
  } else {
    # check for same T, P as last time and keep values
    keep.prop <- character()
    if(length(T)>25) {
      quiet <- FALSE 
      if(identical(thermo$water$T,T) & identical(thermo$water$P,P)) {
        keep.prop <- which(Property %in% names(thermo$water))
        if(length(keep.prop) > 0)
          cat(paste('water: keeping ',length(T),' values for ',c2s(Property[keep.prop]),'.\n',sep=''))
      } else thermo$water <<- NULL
      if(is.null(thermo$water)) thermo$water <<- list(T=T,P=P)
      if(length(keep.prop)!=length(property)) cat(paste('water: calculating',length(T),'values for'))
    } else quiet <- TRUE
  }
  ### the rest is all about working the properties from IAPWS ###
  rho <- function() {
    # return a density in kg m-3
    # corresponding to the given pressure (MPa) and temperature (K)
    pfun <- function(rho,T,P) {
      P <- convert(P,'MPa')
      t <- water.IAPWS95('p',rho=rho,T=T)[,1] - P
      return(t)
    }
    t <- numeric() 
    for(i in 1:length(T)) {
      if(T[i] < 647.096) {
        rho.lower <- water.WP02('rho.liquid',T=T[i])-2
        rho.upper <- rho.lower + 400
        if(P[i] < 5000) rho.upper <- rho.lower + 300
        if(P[i] < 1000) rho.upper <- rho.lower + 200
        if(P[i] < 300) rho.upper <- rho.lower + 30
      }
      else { rho.lower <- 0.01; rho.upper <- 1200}
      tu <- try(uniroot(pfun,c(rho.lower,rho.upper),T=T[i],P=P[i])$root,TRUE)
      if(!is.numeric(tu)) {
        warning('water: no root for density between ',round(rho.lower,1),
        ' and ',round(rho.upper,1),' kg m-3 at ',T[i],' K and ',P[i],' bar.',call.=FALSE,immediate.=TRUE)
        tu <- NA
      }
      t <- c(t,tu)
    }
    return(t)
  }
  v <- function() return(M*1000/my.rho)
  p <- function() return(P)
  ## thermodynamic properties
  # convert to SUPCRT reference state
  # at the triple point
  # I2S = SUPCRT - IAPWS ( + entropy in G )
  dH <- -68316.76 - 451.75437
  dS <- 16.7123 - 1.581072
  dG <- -56687.71 + 19.64228 - dS * (T - thermo$opt$Tr)
  # does the reference state used for GHS also go here?
  dU <- -67434.5 - 451.3229
  dA <- -55814.06 + 20.07376 - dS * (T - thermo$opt$Tr)
  # convert IAPWS95() (specific, joule) to (molar, cal) 
  s <- function()
    return(convert(water.IAPWS95('s',T=T,rho=my.rho)$s*M,'cal')+dS) 
  # u (internal energy) is not here because the letter
  # is used to denote one of the Born functions
  # scratch that! let's put u here and call the other one uborn
  u <- function()
    return(convert(water.IAPWS95('u',T=T,rho=my.rho)$u*M,'cal')+dU)
  a <- function()
    return(convert(water.IAPWS95('a',T=T,rho=my.rho)$a*M,'cal')+dA)
  h <- function() 
    return(convert(water.IAPWS95('h',T=T,rho=my.rho)$h*M,'cal')+dH) 
  g <- function() 
    return(convert(water.IAPWS95('g',T=T,rho=my.rho)$g*M,'cal')+dG) 
  cv <- function() 
    return(convert(water.IAPWS95('cv',T=T,rho=my.rho)$cv*M,'cal')) 
  cp <- function() 
    return(convert(water.IAPWS95('cp',T=T,rho=my.rho)$cp*M,'cal')) 
  w <- function()
    return(water.IAPWS95('w',T=T,rho=my.rho)$w*100) # to cm/s
  ## electrostatic properties
  epsilon <- function() return(water.AW90(T=T,rho=my.rho,P=convert(P,'MPa')))
  de.dt <- function() {
    p <- numeric()
    for(i in 1:length(T)) {
      this.T <- T[i]
      this.P <- P[i]
      this.rho <- my.rho[i]
      dt <- 0.001; t1 <- this.T-dt; t2 <- this.T+dt
      rho <- water('rho',T=c(t1,t2),P=this.P)[,1]
      e <- water.AW90(T=c(t1,t2),rho=rho,rep(this.P,2))
      p <- c(p,(e[2]-e[1])/(2*dt))
    }
    return(p)
  }
  de.dp <- function() {
    p <- numeric()
    for(i in 1:length(T)) {
      this.T <- T[i]
      this.P <- P[i]
      this.rho <- my.rho[i]
      dp <- 0.001; p1 <- this.P-dp; p2 <- this.P+dp
      rho <- water('rho',T=this.T,P=c(p1,p2))[,1]
      e <- water.AW90(P=c(p1,p2),rho=rho,T=rep(this.T,2))
      p <- c(p,(e[2]-e[1])/(2*dp))
    }
    return(p)
  }
  ## Born functions
  q <- function() {
    p <- numeric()
    for(i in 1:length(T)) {
      this.T <- T[i]; this.P <- P[i]; this.rho <- my.rho[i]
      dp <- 0.01; p1 <- this.P-dp; p2 <- this.P+dp
      rho <- water('rho',T=rep(this.T,2),P=c(p1,p2))[,1]
      e <- water.AW90(T=rep(this.T,2),rho=rho,P=convert(c(p1,p2),'MPa'))
      #p <- c(p,convert(-(1/e[2]-1/e[1])/(2*dp),'cm3bar'))
      p <- c(p,-(1/e[2]-1/e[1])/(2*dp))
    }
    return(p)
  }
  n <- function() {
    p <- numeric()
    for(i in 1:length(T)) {
      this.T <- T[i]; this.P <- P[i]; this.rho <- my.rho[i]
      dp <- 0.01; p1 <- this.P-dp; p2 <- this.P+dp
      rho <- water('rho',T=rep(this.T,3),P=c(p1,this.P,p2))[,1]
      e <- water.AW90(T=rep(this.T,3),rho=rho,P=convert(c(p1,this.P,p2),'MPa'))
      #p <- c(p,convert(convert((-(1/e[3]-1/e[2])/dp+(1/e[2]-1/e[1])/dp)/dp,'cm3bar'),'cm3bar'))
      p <- c(p,(-(1/e[3]-1/e[2])/dp+(1/e[2]-1/e[1])/dp)/dp)
    }
    return(p)
  }
  y <- function() {
    p <- numeric()
    for(i in 1:length(T)) {
      this.T <- T[i]; this.P <- P[i]; this.rho <- my.rho[i]
      dt <- 0.001; t1 <- this.T-dt; t2 <- this.T+dt
      rho <- water('rho',T=c(t1,t2),P=rep(this.P,2))[,1]
      e <- water.AW90(T=c(t1,t2),rho=rho,P=convert(rep(this.P,2),'MPa'))
      p <- c(p,-(1/e[2]-1/e[1])/(2*dt))
    }
    return(p)
  }
  x <- function() {
    p <- numeric()
    for(i in 1:length(T)) {
      this.T <- T[i]; this.P <- P[i]; this.rho <- my.rho[i]
      dt <- 0.001; t1 <- this.T-dt; t2 <- this.T+dt
      rho <- water('rho',T=c(t1,this.T,t2),P=rep(this.P,3))[,1]
      e <- water.AW90(T=c(t1,this.T,t2),rho=rho,P=convert(rep(this.P,3),'MPa'))
      p <- c(p,(-(1/e[3]-1/e[2])/dt+(1/e[2]-1/e[1])/dt)/dt)
    }
    return(p)
  }
  uborn <- function() {
    p <- numeric()
    for(i in 1:length(T)) {
      this.T <- T[i]; this.P <- P[i]; this.rho <- my.rho[i]
      dt <- 0.001; this.T1 <- this.T - dt; this.T2 <- this.T + dt
      dp <- 0.001; p1 <- this.P-dp; p2 <- this.P+dp
      rho1 <- water('rho',T=rep(this.T1,2),P=c(p1,p2))[,1]
      rho2 <- water('rho',T=rep(this.T2,2),P=c(p1,p2))[,1]
      e1 <- water.AW90(T=rep(this.T1,2),rho=rho1,P=convert(c(p1,p2),'MPa'))
      e2 <- water.AW90(T=rep(this.T2,2),rho=rho2,P=convert(c(p1,p2),'MPa'))
      #p1 <- convert(-(1/e1[2]-1/e1[1])/(2*dp),'cm3bar')
      #p2 <- convert(-(1/e2[2]-1/e2[1])/(2*dp),'cm3bar')
      p1 <- -(1/e1[2]-1/e1[1])/(2*dp)
      p2 <- -(1/e2[2]-1/e2[1])/(2*dp)
      p <- c(p,(p2-p1)/(2*dt))
    }
    return(p)
  }
  ### main loop; init dataframe output and density holders
  w.out <- NULL; my.rho <- NULL
  # get densities and tell about it
  if(!quiet) {
    if('rho' %in% names(thermo$water)) my.rho <- thermo$water[(names(thermo$water)=='rho')][[1]]
    else { 
      cat(' rho'); my.rho <- rho() 
      thermo$water <<- c(thermo$water,list(tmp=my.rho))
      names(thermo$water)[length(thermo$water)] <<- 'rho'
    }
  } else my.rho <- rho()
  for(i in 1:length(property)) {
    if(property[i] %in% c('e','kt')) {
      # expansivity isn't in the table yet... set it to zero
      warning('water: values of ',Property[i],' are NA\n',call.=FALSE)
      inew <- rep(NA,length(T))
    } else {
      if(!quiet) {
        if(Property[i] %in% names(thermo$water)) inew <- thermo$water[(names(thermo$water)==Property[i])][[1]]
        else { 
          cat(paste(' ',Property[i],sep='')); inew <- get(property[i])() 
          thermo$water <<- c(thermo$water,list(tmp=inew))
          names(thermo$water)[length(thermo$water)] <<- Property[i]
        }
      } else inew <- get(property[i])()
    }
    #if(NA %in% inew) na.h2o <- TRUE
    wnew <- data.frame(inew)
    if(i > 1) w.out <- cbind(w.out,wnew) else w.out <- wnew
  }  
  colnames(w.out) <- Property
  if(!quiet & length(keep.prop)!=length(property)) cat('.\n')
  return(w.out)
}


water.SUPCRT92 <- function(property,T=298.15,P=1,isat=0) {
  ### interface to H2O92D.f : FORTRAN subroutine taken from 
  ### SUPCRT92 for calculating the thermodynamic and 
  ### electrostatic properties of H2O. 
  ## we restrict the calculations to liquid water
  ## except for getting Psat (vapor-liquid saturation 
  ## pressure as a function of T>100 C). 20071213 jmd
  # H2O92 doesn't output Born functions N or U
  if('n' %in% tolower(property) | 'uborn' %in% tolower(property))
    stop('I can\'t tell you the Born functions N or U (used in calculating compressibilities and expansibilities of aqueous species).')
  # pressure setting
  if(is.null(P)) P <- rep(0,length(T))
  # values to use here gleaned from H2O92D.f and SUP92D.f
  # it, id, ip, ih, itripl, isat, iopt, useLVS, epseqn, icrit
  if(isat) iopt <- 1 else iopt <- 2  # for Psat(T) (1) or T-P (2)
  specs <- c(2,2,2,5,1,isat,iopt,1,4,0)
  states <- rep(0,4)
  # match up properties with the output
  props <- c('a','g','s','u','h','cv','cp','Speed','alpha',
    'beta','diel','visc','tcond','surten','tdiff','Prndtl',
    'visck','albe','ZBorn','YBorn','QBorn','daldT','XBorn')
  iprop <- seq(1,45,length.out=23)
  # check for same T, P as last time
  length.keep <- 25; do.keep <- FALSE
  if(length(T) > length.keep) {
    if(isat) {
      if(identical(thermo$Psat$T,T)) do.keep <- TRUE
      else thermo$Psat <<- list(T=T)
    } else {
      if(identical(thermo$water2$P,P) & identical(thermo$water2$T,T)) do.keep <- TRUE
      else thermo$water2 <<- list(T=T,P=P)
    }
  } 
  # now to the actual calculations
  if(!do.keep) {
    Tc <- convert(T,'C')
    # initialize the output matrix
    w.out <- matrix(NA,nrow=length(T),ncol=23,byrow=TRUE) 
    err.out <- numeric(length(T))
    rho.out <- numeric(length(T))
    p.out <- numeric(length(T))
    # 20091022 TODO: parallelize this
    for(i in 1:length(T)) {
      states[1] <- Tc[i]
      states[2] <- P[i]
      if(any(is.na(c(Tc[i],P[i])))) {
        # if T or P is NA, all properties are NA
        w <- matrix(rep(NA,23),nrow=1)
        w.out[i,] <- w
        p.out[i] <- NA
        rho.out[i] <- NA
      } else {
        inc <- 0
        h2o <- .Fortran('H2O92',as.integer(specs),as.double(states),
          as.double(rep(0,46)),as.integer(0),PACKAGE='CHNOSZ')
        # errors
        err <- h2o[[4]]
        err.out[i] <- err
        # density
        rho <- h2o[[2]][3]
        rho2 <- h2o[[2]][4]
        if(rho2 > rho) {
          # liquid is denser than vapor
          rho <- rho2 
          # for selecting the liquid properties later
          inc <- 1
        }
        rho.out[i] <- rho
        # most of the properties we're interested in
        w <- t(h2o[[3]][iprop+inc])
        if(err==1) w[1,] <- NA
        # update the ith row of the output matrix
        w.out[i,] <- w
        # Psat
        if(isat | 'psat' %in% tolower(property)) {
          p <- h2o[[2]][2]
          p[p==0] <- NA
          # Psat specifies P=1 below 100 degC
          p[p < 1] <- 1
          p.out[i] <- p
        } else {
          p.out[i] <- P[i]
        }
      }
    }
    # convert output to dataframe
    w.out <- as.data.frame(w.out)
    names(w.out) <- props
    # assemble the properties
    mwH2O <- 18.0152 # SUP92.f
    w.out <- cbind(w.out,V=mwH2O/rho.out,rho=rho.out*1000)
    if(isat | 'psat' %in% tolower(property)) w.out <- cbind(w.out,Psat=p.out)
    # keep the calculated properties around for future reference
    if(length(T) > length.keep) {
      if(isat) thermo$Psat <<- c(thermo$Psat,list(out=w.out))
      else thermo$water2 <<- c(thermo$water2,list(out=w.out))
    }
    # tell the user about any problems
    if(any(err.out==1)) {
      if(length(T) > 1) plural <- "s" else plural <- ""
      nerr <- length(which(err.out==1))
      if(nerr > 1) plural2 <- "s" else plural2 <- ""
      if(isat) msgout(paste("water.SUPCRT92: error",plural2," calculating ",
        nerr," of ",length(T)," point",plural,"; for Psat we need T < 647.067 K\n",sep=""))
      else msgout(paste("water.SUPCRT92: error",plural2," calculating ",nerr,
        " of ",length(T)," point",plural,
        "; T < Tfusion@P, T > 2250 degC, or P > 30kb.\n",sep=""))
        # that last bit is taken from SUP92D.f in the SUPCRT92 distribution
    }
  } else {
    if(isat) w.out <- thermo$Psat$out
    else w.out <- thermo$water2$out
  }
  # if isat is 1, just return the calculated pressures
  if(isat) return(w.out$Psat)
  # return only the selected properties
  icol <- match(tolower(property),tolower(colnames(w.out)))
  return(w.out[,icol,drop=FALSE])
}

