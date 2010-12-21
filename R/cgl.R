# CHNOSZ/cgl.R
# calculate standard thermodynamic properties of non-aqueous species
# 20060729 jmd

cgl <- function(property=NULL,T=298.15,P=1,ghs=NULL,eos=NULL) {
  # calculate properties of crystalline, liquid (except H2O) and gas species
  # argument handling
  Tr <- thermo$opt$Tr; Pr <- thermo$opt$Pr
  eargs <- eos.args('mk',property=property)
  prop <- eargs$prop
  props <- eargs$props
  Prop <- eargs$Prop

  # a list - the result
  x <- list()
  for(k in 1:nrow(ghs)) {
    # loop over each species
    GHS <- ghs[k,]
    EOS <- eos[k,]
    w <- NULL
    for(i in 1:length(prop)) {
      property <- prop[i]
      # equations for lambda adapted from HOK+98
      if(property=='cp')
        p <- EOS$a + EOS$b*T + EOS$c*T^-2 + EOS$d*T^-0.5 + EOS$e*T^2 + EOS$f*T^EOS$lambda
      if(property=='v')
        p <- rep(EOS$V,length(T))
      if(property %in% c('e','kt')) {
        p <- rep(NA,length(T))
        warning('cgl: E and/or kT of cr, gas and/or liq species are NA.')
      }
      if(property=='g') {
        p <-   0 - GHS$S*(T-Tr) + EOS$a*(T-Tr-T*log(T/Tr)) - 
               EOS$b*(T-Tr)^2/2 - EOS$c*(1/T + T/Tr^2 - 2/Tr)/2 -
               EOS$d*(T^0.5-0.5*T*Tr^-0.5-0.5*Tr^0.5)/-0.25 -
               EOS$e*(T^3-3*T*Tr^2+2*Tr^3)/6 +
               convert(EOS$V*(P-Pr),'calories')
        p[is.na(p)] <- 0
        p <- GHS$G + p
        if(!is.na(EOS$f) & !is.na(EOS$lambda)) if(EOS$f!=0) {
          if(EOS$lambda== -1) p <- p + EOS$f*(log(T/Tr)-T*(1/Tr-1/T))
          else p <- p + EOS$f * ( T^(EOS$lambda+1) - (EOS$lambda+1)*T*Tr^EOS$lambda + 
            EOS$lambda*Tr^(EOS$lambda+1) ) / ( EOS$lambda*(EOS$lambda+1) ) 
        }
      }
      if(property=='h') { 
        p <- EOS$a*(T-Tr) + EOS$b*(T^2-Tr^2)/2 +
             EOS$c*(1/T-1/Tr)/-1 + EOS$d*(T^0.5-Tr^0.5)/0.5 + 
             EOS$e*(T^3-Tr^3)/3 
             # SUPCRT seems to ignore this term? ... 20070802
             # + convert(EOS$V*(P-Pr),'calories')
        p[is.na(p)] <- 0
        p <- GHS$H + p
        if(!is.na(EOS$f) & !is.na(EOS$lambda)) if(EOS$f!=0) {
           if(EOS$lambda == -1) p <- p + EOS$f*log(T/Tr) 
           else p <- p - EOS$f * ( T^(EOS$lambda+1) - Tr^(EOS$lambda+1) ) / (EOS$lambda+1)
        }
      }
      if(property=='s') {
        p <- EOS$a*log(T/Tr) + EOS$b*(T-Tr) + 
             EOS$c*(T^(-2)-Tr^(-2))/(-2) + EOS$e*(T^2-Tr^2)/2 + 
             EOS$d*(T^-0.5-Tr^-0.5)/-0.5
        p[is.na(p)] <- 0
        p <- GHS$S + p
        if(!is.na(EOS$f) & !is.na(EOS$lambda)) if(EOS$f!=0) {
          p <- p + EOS$f*(T^EOS$lambda-Tr^EOS$lambda)/EOS$lambda
        }
      }
      wnew <- data.frame(p)
      if(i>1) w <- cbind(w,wnew) else w <- wnew
    }
  colnames(w) <- Prop
  x[[k]] <- w
 }
 return(x)
}

