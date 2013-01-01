# CHNOSZ/EOSregress.R  
# model volumes and heat capacities of aqueous species
# 20091105 first version
# 20110429 revise and merge with CHNOSZ package

EOSvar <- function(var,T,P) {
  # get the variables of a term in a regression equation
  # T (K), P (bar)
  out <- switch(EXPR = var,
    "(Intercept)" = rep(1,length(T)),
    "T" = T,
    "P" = P,
    "TTheta" = T-thermo$opt$Theta,                 # T-Theta
    "invTTheta" = (T-thermo$opt$Theta)^-1,         # 1/(T-Theta)
    "TTheta2" = (T-thermo$opt$Theta)^2,            # (T-Theta)^2
    "invTTheta2" = (T-thermo$opt$Theta)^-2,        # 1/(T-Theta)^2
    "invPPsi" = (P+thermo$opt$Psi)^-1,             # 1/(P-Psi)
    "invPPsiTTheta" = (P+thermo$opt$Psi)^-1 * (T-thermo$opt$Theta)^-1,  # 1/[(P-Psi)(T-Theta)]
    "V" = water(var,T=T,P=P)[,1],
    "E" = water(var,T=T,P=P)[,1],
    "kT" = water(var,T=T,P=P)[,1],
    "alpha" = water(var,T=T,P=P)[,1],
    "beta" = water(var,T=T,P=P)[,1],
    "X" = water(var,T=T,P=P)[,1],
    "Q" = water(var,T=T,P=P)[,1],
    "TX" = T*water("X",T=T,P=P)[,1],
    "drho.dT" = -water("rho",T=T,P=P)[,1]*water("E",T=T,P=P)[,1],
    "V.kT" = water("V",T=T,P=P)[,1]*water("kT",T=T,P=P)[,1],
    NA
  )
  return(out)
}

EOSlab <- function(var,coeff="") {
  # make pretty labels for the variables
  lab <- switch(EXPR = var,
    "(Intercept)" = substitute(YYY*" ",list(YYY=coeff)),
    "TTheta" = substitute(YYY%*%(italic(T)-Theta),list(YYY=coeff)),
    "invTTheta" = substitute(YYY/(italic(T)-Theta),list(YYY=coeff)),
    "TTheta2" = substitute(YYY%*%(italic(T)-Theta)^2,list(YYY=coeff)),
    "invTTheta2" = substitute(YYY/(italic(T)-Theta)^2,list(YYY=coeff)),
    "T" = substitute(YYY%*%italic(XXX),list(XXX=var,YYY=coeff)),
    "P" = substitute(YYY%*%italic(XXX),list(XXX=var,YYY=coeff)),
    "V" = substitute(YYY%*%italic(XXX),list(XXX=var,YYY=coeff)),
    "E" = substitute(YYY%*%italic(XXX),list(XXX=var,YYY=coeff)),
    "X" = substitute(YYY%*%italic(XXX),list(XXX=var,YYY=coeff)),
    "Q" = substitute(YYY%*%italic(XXX),list(XXX=var,YYY=coeff)),
    "TX" = substitute(YYY%*%italic(XXX),list(XXX=var,YYY=coeff)),
    "kT" = substitute(YYY%*%kappa[italic(T)],list(YYY=coeff)),
    "alpha" = substitute(YYY%*%alpha,list(YYY=coeff)),
    "beta" = substitute(YYY%*%beta,list(YYY=coeff)),
    "drho.dT" = substitute(YYY%*%(d~rho/dT),list(YYY=coeff)),
    "V.kT" = substitute(YYY%*%V~kappa[italic(T)],list(YYY=coeff)),
    NA
  )
  return(lab)
}

EOSregress <- function(exptdata,var="",T.max=9999) {
  # regress exptdata using terms listed in fun 
  # which values to use
  iT <- which(exptdata$T <= T.max)
  exptdata <- exptdata[iT,]
  # temperature and pressure
  T <- exptdata$T
  P <- exptdata$P
  # the third column is the property of interest: Cp or V
  X <- exptdata[,3]
  # now build a regression formula 
  if(length(var) == 0) stop("var is missing")
  fmla <- as.formula(paste("X ~ ",paste(var,collapse="+")))
  # retrieve the values of the variables
  for(i in seq_along(var)) assign(var[i],EOSvar(var[i],T=T,P=P))
  # now regress away!
  EOSlm <- lm(fmla)
  return(EOSlm)
}

EOScalc <- function(coefficients,T,P) {
  # calculate values of volume
  # or heat capacity from regression fit
  X <- 0
  for(i in 1:length(coefficients)) {
    coeff.i <- coefficients[[i]]
    fun.i <- EOSvar(names(coefficients)[i],T,P)
    X <- X + coeff.i * fun.i
  }
  return(X)
}

EOSplot <- function(exptdata,var=NULL,T.max=9999,T.plot=NULL,
  P=NULL,fun.legend="topleft",coefficients=NULL) {
  # plot experimental and modelled volumes and heat capacities
  # first figure out the property (Cp or V) from the exptdata
  prop <- colnames(exptdata)[3]
  # if var is NULL use HKF equations
  if(is.null(var)) {
    if(prop=="Cp") var <- c("invTTheta2","TX")
    if(prop=="V") var <- c("invTTheta","Q")
  }
  expt <- exptdata
  # perform the regression, only using temperatures up to T.max
  if(is.null(coefficients)) {
    EOSlm <- EOSregress(expt,var,T.max)
    coefficients <- EOSlm$coefficients
  }
  # only plot points below a certain temperature
  iexpt <- 1:nrow(expt)
  if(!is.null(T.plot)) iexpt <- which(expt$T < T.plot)
  iX <- match(prop,colnames(expt))
  ylim <- extendrange(expt[iexpt,iX],f=0.1)
  xlim <- extendrange(expt$T[iexpt],f=0.1)
  # start plot
  thermo.plot.new(xlim=xlim,ylim=ylim,xlab=axis.label("T", units="K"),
    ylab=axis.label(paste(prop,"0",sep="")),yline=2,mar=NULL)
  # we group the data by pressure ranges;
  # assume increasing temperatures are in the
  # same pressure range but a decrease in temperature
  # signals the next pressure range
  idrop <- c(1,which(diff(expt$T)<0)+1,length(expt$T)+1)
  Plab <- character()
  pch.open <- c(1,0,2)
  pch.filled <- c(16,15,17)
  for(i in 1:(length(idrop)-1)) {
    ip <- idrop[i]:(idrop[i+1]-1)
    # find the calculated values at these conditions
    myT <- expt$T[ip]
    myP <- expt$P[ip]
    calc.X <- EOScalc(coefficients,myT,myP)
    expt.X <- expt[ip,iX]
    # are we within 10% of the values
    in10 <- which(abs((calc.X-expt.X)/expt.X) < 0.1)
    pch <- rep(pch.open[i],length(myT))
    pch[in10] <- pch.filled[i]
    points(myT,expt[ip,iX],pch=pch)
    # if we calculate lines at a constant P, do that
    xs <- seq(xlim[1],xlim[2],length.out=200)
    if(!is.null(P)) {
      myT <- xs
      myP <- P
      calc.X <- EOScalc(coefficients,myT,myP)
    } 
    # take out NAs and Infinite values
    iNA <- is.na(calc.X) | is.infinite(calc.X)
    xs <- xs[!iNA]
    calc.X <- calc.X[!iNA]
    myT <- myT[!iNA]
    # plot regression line
    lines(xs,splinefun(myT,calc.X,method="monoH.FC")(xs))
    Plim <- range(expt$P[ip])
    Plab <- c(Plab,paste(Plim[1],"-",Plim[2],"bar"))
  }
  # make legend
  if(!is.null(fun.legend)) {
    coeffs <- as.character(round(as.numeric(coefficients),4))
    # so that positive ones appear with a plus sign
    ipos <- which(coeffs >= 0)
    coeffs[ipos] <- paste("+",coeffs[ipos],sep="")
    # make labels for the functions
    fun.lab <- as.expression(lapply(1:length(coeffs),
      function(x) {EOSlab(names(coefficients)[x],coeffs[x])} ))
    #fun.lab <- paste(names(coeffs),round(as.numeric(coeffs),4))
    legend(fun.legend,legend=fun.lab,pt.cex=0.1)
  }
  return(invisible(list(xlim=range(expt$T[iexpt]))))
}

EOScoeffs <- function(species,property) {
  # get the HKF coefficients for species in the database
  iis <- info(info(species,"aq"))
  if(property=="Cp") {
    out <- iis[,c("c1","c2","omega")]
    names(out) <- c("(Intercept)","invTTheta2","TX")
  } else if(property=="V") {
    iis <- iis[,c("a1","a2","a3","a4","omega")]
    sigma <- ( iis$a1 + iis$a2 / (2600 + 1) ) * 41.84
    xi <- ( iis$a3 + iis$a4 / (2600 + 1) ) * 41.84
    # watch for the negative sign on omega here!
    out <- data.frame(sigma,xi,-iis$omega)
    names(out) <- c("(Intercept)","invTTheta","Q")
  }
  return(out)
}

