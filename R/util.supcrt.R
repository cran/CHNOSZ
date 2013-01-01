# CHNOSZ/util.supcrt.R
# import and export thermodynamic data in SUPCRT format (e.g., slop98.dat)
# 20051105 jmd


read.supcrt <- function(file) {
#cat(paste('read.supcrt: reading',file,'\n'))

  # the file must not have any comments in the data blocks
  # (the reference and data must be on a line by themselves)
  # 999999 will be converted to NA

  # read the entire file
  tab <- scan(file,what='character')
  # scanning the file fooled by an entry starting with
  # RIBOSE-5-PHOSPHATE-2C5H9O8P2-
  # which really means "RIBOSE-5-PHOSPHATE-2" "C5H9O8P2-"


  # function to identify separators
  # (lines that begin with ******)
  is.separator <- function(s) {
    if( nchar(s) > 5 & identical(substr(s,1,1),'*') ) return(TRUE)
    return(FALSE)
  }

  # function to find the first data position
  # for our desired state
  ifirst <- function(tab,state,transitions=0) {
    isep <- which(substr(tab,1,5)=="*****")
    isep <- isep[!isep==1]
    if( state=='aq' ) return(isep[match("species",tab[isep-1])]+1)
    if( state=='gas' ) return(isep[match("gases",tab[isep-1])]+1)
    if( state=="cr" ) {
      # what word we search for
      txt <- (c("undergo","one","two","three"))[transitions+1]
      return(isep[match(txt,tab[isep-3])]+1)
    }
  }

  # a function that strips trailing decimal points
  # (so e.g. -66.310. typo Gf for LYSINE- in SLOP07.dat can be processed)
  striptrailingdot <- function(x) {
    istart <- rep(1,length(x))
    iend <- nchar(x)
    hastrailingdot <- grep("\\.$",x)
    if(length(hastrailingdot) > 0) iend[hastrailingdot] <- iend[hastrailingdot] - 1
    x <- substr(x,istart,iend)
    return(x)
  }

  name.lower <- function(name,abbrv,formula,state) {
    name <- as.character(name)
    abbrv <- as.character(abbrv)
    formula <- as.character(formula)
    for(i in 1:length(name)) {
      if( state=='aq' & (identical(name[i],formula[i]) 
        | identical(name[i],abbrv[i]) | identical(abbrv[i],formula[i]) ) ) next
      if( state=='cr' & (identical(name[i],formula[i]) 
        | identical(name[i],abbrv[i])) ) next
      name[i] <- tolower(name[i])
    }
    return(name)
  }
  
  strip <- function(s) {
    for(i in 1:length(s)) {
      si <- s[i]
      if(nchar(si)>3) {
        ss <- substr(si,nchar(si)-2,nchar(si))
        if(ss %in% c('(g)',',aq',',AQ','(s)','(-)','(+)','(0)','(1)','(2)') ) {
          s[i] <- substr(si,1,nchar(si)-3)
          if(identical(ss,'(-)')) s[i] <- paste(s[i],'-',sep='')
          if(identical(ss,'(+)')) s[i] <- paste(s[i],'+',sep='')
          if(identical(ss,'(1)')) s[i] <- paste(s[i],'+1',sep='')
          if(identical(ss,'(2)')) s[i] <- paste(s[i],'+2',sep='')
        }
      }
      if(nchar(si)>4) {
        ss <- substr(si,nchar(si)-3,nchar(si)) 
        if(ss %in% c('(+0)','(+1)','(+2)','(+3)','(+4)',
          '(-1)','(-2)','(-3)','(aq)','.dat') ) {
          s[i] <- substr(si,1,nchar(si)-4)
          if(identical(ss,'(-3)')) s[i] <- paste(s[i],'-3',sep='')
          if(identical(ss,'(-2)')) s[i] <- paste(s[i],'-2',sep='')
          if(identical(ss,'(-1)')) s[i] <- paste(s[i],'-1',sep='')
          if(identical(ss,'(+1)')) s[i] <- paste(s[i],'+1',sep='')
          if(identical(ss,'(+2)')) s[i] <- paste(s[i],'+2',sep='')
          if(identical(ss,'(+3)')) s[i] <- paste(s[i],'+3',sep='')
          if(identical(ss,'(+4)')) s[i] <- paste(s[i],'+4',sep='')
        }
      }
    }
    return(s)
  }

  source <- toupper(strip(file))
  
  read.cr.transitions <- function(tab,istart,transitions) {
    ntail <- c(8,15,22,29)[transitions+1]
    gsd <- get.state.data(tab,istart,nhead=6,ntail=ntail)
    cr <- gsd$gas
    ghs <- data.frame( name=name.lower(strip(cr[[1]]),strip(cr[[3]]),
      strip(cr[[2]]),'cr'), abbrv=cr[[3]], formula=strip(cr[[2]]),
      state=rep('cr',length(cr[[1]])), source=strip(cr[[5]]),
      date=strip(cr[[6]]), Gf=cr[[7]], Hf=cr[[8]], S=cr[[9]] )
    eos <- data.frame( name=ghs$name, source=strip(cr[[5]]),
      a=cr[[11]], b=cr[[12]], c=cr[[13]], V=cr[[10]], T=cr[[14]] )
    trans <- data.frame()
    for(j in 0:transitions) {
      if(identical(transitions,0)) break
      it <- 1:7 + 7 + 7 * j
      t <- data.frame( name=ghs$name, Htr=cr[[it[1]]], Vtr=cr[[it[2]]], dP.dT=cr[[it[3]]],
        a=cr[[it[4]]], b=cr[[it[5]]], c=cr[[it[6]]], T=cr[[it[7]]] )
      trans <- rbind(trans,t)
    }
    return(list(ghs=ghs,eos=eos,transitions=trans,ilast=gsd$ilast))
  }

  read.cr <- function(tab) {
    icr0 <- ifirst(tab,'cr')
    cr0 <- read.cr.transitions(tab,icr0,0)
    icr1 <- ifirst(tab,'cr',1)
    cat(paste('read.supcrt:',nrow(cr0$ghs),
      'minerals that do not undergo phase transition\n'))
    cr1 <- read.cr.transitions(tab,icr1,1)
    icr2 <- ifirst(tab,'cr',2)
    cat(paste('read.supcrt:',nrow(cr1$ghs),
      'minerals that undergo one phase transition\n'))
    cr2 <- read.cr.transitions(tab,icr2,2)
    icr3 <- ifirst(tab,'cr',3)
    cat(paste('read.supcrt:',nrow(cr2$ghs),
      'minerals that undergo two phase transitions\n'))
    cr3 <- read.cr.transitions(tab,icr3,3)
    cat(paste('read.supcrt:',nrow(cr3$ghs),
      'minerals that undergo three phase transitions\n'))
    ghs <- rbind(cr0$ghs,cr1$ghs,cr2$ghs,cr3$ghs)
    eos <- rbind(cr0$eos,cr1$eos,cr2$eos,cr3$eos)
    transitions <- rbind(cr0$transition,cr1$transitions,cr2$transitions,cr3$transitions)
    return(list(ghs=ghs,eos=eos,transitions=transitions,ilast=cr3$ilast))
  }

  get.state.data <- function(tab,istart,nhead=6,ntail=8) {
    # get all the data for particular section
    gas <- list()
    gas[[nhead+ntail]] <- numeric()
    # populating the list
    inext <- istart
    for(i in istart:length(tab)) {
      if(is.separator(tab[i])) break
      # don't go into the footer of the file
      if(tab[i+1]=="minerals" & tab[i+2]=="that") break
      if(i==inext) {
        # entries start with nhead identifying fields
        # followed by ntail numeric (thermodynamic data) fields
        for(j in 1:nhead) gas[[j]] <- c(gas[[j]],tab[i+j-1])
      } else if((i > (inext+nhead-1)) & can.be.numeric(striptrailingdot(tab[i]))) {
        for(j in 1:ntail) gas[[j+nhead]] <- c(gas[[j+nhead]],tab[i+j-1])
        inext <- i+j
      }
    }
    return(list(gas=gas,ilast=i))
  }

  read.gas <- function(tab,istart) {
    gsd <- get.state.data(tab,istart,nhead=6,ntail=8)
    gas <- gsd$gas
    # processing the list
    ghs <- data.frame( name=tolower(gas[[2]]), 
      abbrv=strip(gas[[3]]), formula=strip(gas[[3]]),
      state=rep('gas',length(gas[[1]])), source=strip(gas[[5]]),
      date=strip(gas[[6]]),Gf=gas[[7]], Hf=gas[[8]], S=gas[[9]] )
    eos <- data.frame( name=ghs$name, source=strip(gas[[5]]),
      a=gas[[11]], b=gas[[12]], c=gas[[13]], V=gas[[10]], T=gas[[14]] )
    return(list(ghs=ghs,eos=eos,ilast=gsd$ilast))
  }

  read.aq <- function(tab,istart) {
    gsd <- get.state.data(tab,istart,6,11)
    aq <- gsd$gas
    ghs <- data.frame( name=strip(aq[[1]]),
      abbrv=strip(aq[[3]]), formula=strip(aq[[2]]),
      state=rep('aq',length(aq[[1]])), source=strip(aq[[5]]),
      date=strip(aq[[6]]), Gf=aq[[7]], Hf=aq[[8]], S=aq[[9]] )
    eos <- data.frame( name=ghs$name, source=strip(aq[[5]]),
      a1=aq[[10]], a2=aq[[11]], a3=aq[[12]], a4=aq[[13]], 
      c1=aq[[14]], c2=aq[[15]], omega=aq[[16]] )
    return(list(ghs=ghs,eos=eos,ilast=gsd$ilast))
  }

  cr <- read.cr(tab)
  igas <- ifirst(tab,'gas')
  gas <- read.gas(tab,igas)
  cat(paste('read.supcrt:',nrow(gas$ghs),'gases\n'))
  iaq <- ifirst(tab,'aq',gas$ilast)
  aq <- read.aq(tab,iaq)
  cat(paste('read.supcrt:',nrow(aq$ghs),'aqueous species\n'))

  supcrt <- list(ghs=rbind(cr$ghs,gas$ghs,aq$ghs),
    eos.cr=cr$eos,eos.gas=gas$eos,eos.aq=aq$eos,transitions=cr$transitions)
  return(supcrt)

  # some numerical conversions
  as.speciate <- function(supcrt) {
    for(i in 1:ncol(supcrt)) {
      supcrt[,i] <- as.numeric(as.character(supcrt[,i]))
      supcrt[supcrt[,i]==999999,i] <- NA
      cname <- colnames(supcrt)[i]
      if(cname=='b') supcrt[,i] <- supcrt[,i] * 10^(-3)
      if(cname=='c') supcrt[,i] <- supcrt[,i] * 10^(5)
      if(cname=='a1') supcrt[,i] <- supcrt[,i] * 10^(-1)
      if(cname=='a2') supcrt[,i] <- supcrt[,i] * 10^(2)
      if(cname %in% c('a4','c2')) supcrt[,i] <- supcrt[,i] * 10^(4)
      if(cname=='omega') supcrt[,i] <- supcrt[,i] * 10^(5)
    }
    return(supcrt)
  }

  supcrt$ghs[,6:8] <- as.speciate(supcrt$ghs[,6:8])
  supcrt$eos.cr[,3:7] <- as.speciate(supcrt$eos.cr[,3:7])
  supcrt$transitions[,2:8] <- as.speciate(supcrt$transitions[,2:8])
  supcrt$eos.gas[,3:7] <- as.speciate(supcrt$eos.gas[,3:7])
  supcrt$eos.aq[,3:9] <- as.speciate(supcrt$eos.aq[,3:9])

  return(supcrt)
}


write.supcrt <- function(file='supcrt.dat',obigt=thermo$obigt) {
  # create a SUPCRT-like data file that can be imported
  # into the OBIGT or SUPCRT programs. originally written
  # for OBIGT import; made to work with SUPCRT/CPRONS May-June 2007
  # set units to a known value
  T.units('K')
  separator <- c2s(rep('*',72),sep='')
  spaces <- function(n) c2s(rep(' ',n),sep='')
  min0text <- 'minerals that do not undergo phase transitions'
  min0textalt <- 'minerals that do not undergo phase transition'
  min1text <- 'minerals that undergo one phase transition'
  min2text <- 'minerals that undergo two phase transitions'
  min3text <- 'minerals that undergo three phase transitions'
  gastext <- 'gases'
  aqtext <- 'aqueous species'
  commtext <- 'comment lines precede top of min1 block'

  Z <- function(obigt) {
    mu <- makeup(as.numeric(rownames(obigt)))
    if('Z' %in% rownames(mu)) Z <- mu$count[rownames(mu)=='Z']
    else Z <- 0
    return(Z)
  }
  NAME <- function(o) {
    n <- as.character(o$name)
    f <- as.character(o$formula)
    # spaces and colons in names make OBIGT unhappy
    if(length(grep(' ',n))>0) {
      n <- c2s(s2c(n,sep=' ',keep.sep=FALSE),sep='-')
    }
    if(length(grep(':',n))>0) {
      n <- c2s(s2c(n,sep=':',keep.sep=FALSE),sep='-')
    }
    # names can be as long as 20 characters for SUPCRT
    # (and run into the formula field) but OBIGT wants a space
    # so cut them off at 19 characters
    if(as.character(o$state)=='liq') {
      if(nchar(n)>17) n <- substr(n,1,17)
      n <- paste(n,',l',sep='')
    } else if(as.character(o$state)=='gas') {
      if(nchar(n)>17) n <- substr(n,1,17)
      n <- paste(n,',g',sep='')
    } else if(nchar(n)>19) n <- substr(n,1,19)
    .space <- c(20,51)
    .length <- c(nchar(n),nchar(f))
    .pad <- .space - .length
    n.pad <- c2s(rep(' ',.pad[1]),sep='')
    f.pad <- c2s(rep(' ',.pad[2]),sep='')
    c2s(c(n,n.pad,f,f.pad),sep='')
  }
  ABBRV <- function(o) {
    f <- as.character(o$formula)
    # rewrite the formula of minerals (asterisks confuse OBIGT)
    if(as.character(o$state)%in% c('cr','cr1','cr2','cr3','cr4')) {
      f <- makeup(f,'')
      # any 'Z's here (as in Pd-oxyannite) will confuse OBIGT
      f <- f[,!colnames(f)%in%'Z',drop=FALSE]
      f <- makeup(f,'')
    }
    abbrv <- o$abbrv
    if(is.na(abbrv)) abbrv <- o$name
    abbrv <- as.character(abbrv)
    # spaces in abbreviations make OBIGT unhappy
    # (it thinks text after a space should be a formula)
    if(length(grep(' ',abbrv))>0) {
      abbrv <- c2s(s2c(abbrv,sep=' ',keep.sep=FALSE),sep='-')
    }
    if(nchar(abbrv)>19) abbrv <- substr(abbrv,1,19)
    .space <- c(20,51)
    .length <- c(nchar(abbrv),nchar(f))
    .pad <- .space - .length
    abbrv.pad <- c2s(rep(' ',.pad[1]),sep='')
    f.pad <- c2s(rep(' ',.pad[2]),sep='')
    c2s(c(abbrv,abbrv.pad,f,f.pad),sep='')
  }
  DATE <- function(o) {
    source <- as.character(o$ref1)
    if(!is.na(o$ref2)) if(!o$ref2=='CHNOSZ') 
      source <- paste(source,',',as.character(o$ref2),sep='')
    date <- as.character(o$date)
    .space <- c(20,51)
    .length <- c(nchar(source),nchar(date))
    .pad <- .space - .length
    source.pad <- c2s(rep(' ',.pad[1]),sep='')
    date.pad <- c2s(rep(' ',.pad[2]),sep='')
    c2s(c(source,source.pad,date,date.pad),sep='')
  }
  newround <- function(x,n) {
    # to round and insert decimal point
    if(is.na(x)) return(x)
    x <- round(x,n)
    if(trunc(x)==x & length(grep('e+',x))==0) x <- paste(x,'.',sep='')
    return(x)
  }
  GHSV <- function(o) {
    G <- newround(o$G,1)
    H <- newround(o$H,1)
    S <- newround(o$S,3)
    V <- newround(o$V,3)
    if(is.na(G)) G <- 999999
    if(is.na(H)) H <- 999999
    if(is.na(S)) S <- 999999
    if(is.na(V)) V <- 999999
    .space <- c(17,14,10,10)
    .length <- c(nchar(G),nchar(H),nchar(S),nchar(V))
    .pad <- .space - .length
    G.pad <- c2s(rep(' ',.pad[1]),sep='')
    H.pad <- c2s(rep(' ',.pad[2]),sep='')
    S.pad <- c2s(rep(' ',.pad[3]),sep='')
    V.pad <- c2s(rep(' ',.pad[4]),sep='')
    c2s(c(G.pad,G,H.pad,H,S.pad,S,V.pad,V),sep='')
  }
  GHS <- function(o) {
    G <- newround(o$G,1)
    H <- newround(o$H,1)
    S <- newround(o$S,3)
    if(is.na(G)) G <- 999999
    if(is.na(H)) H <- 999999
    if(is.na(S)) S <- 999999
    .space <- c(15,12,12)
    .length <- c(nchar(G),nchar(H),nchar(S))
    .pad <- .space - .length
    G.pad <- c2s(rep(' ',.pad[1]),sep='')
    H.pad <- c2s(rep(' ',.pad[2]),sep='')
    S.pad <- c2s(rep(' ',.pad[3]),sep='')
    c2s(c(G.pad,G,H.pad,H,S.pad,S),sep='')
  }
  HKF1 <- function(o) {
    a1 <- o$a1.a; a2 <- o$a2.b; a3 <- o$a3.c; a4 <- o$a4.d
    if(is.na(a1)) a1 <- 999999 else a1 <- newround(a1,4)
    if(is.na(a2)) a2 <- 999999 else a2 <- newround(a2,4)
    if(is.na(a3)) a3 <- 999999 else a3 <- newround(a3,4)
    if(is.na(a4)) a4 <- 999999 else a4 <- newround(a4,4)
    .space <- c(13,12,12,12)
    .length <- c(nchar(a1),nchar(a2),nchar(a3),nchar(a4))
    .pad <- .space - .length
    a1.pad <- c2s(rep(' ',.pad[1]),sep='')
    a2.pad <- c2s(rep(' ',.pad[2]),sep='')
    a3.pad <- c2s(rep(' ',.pad[3]),sep='')
    a4.pad <- c2s(rep(' ',.pad[4]),sep='')
    c2s(c(a1.pad,a1,a2.pad,a2,a3.pad,a3,a4.pad,a4),sep='')
  }
  HKF2 <- function(o) {
    c1 <- o$c1.e; c2 <- o$c2.f; omega <- o$omega.lambda; z <- o$z.T
    if(is.na(c1)) c1 <- 999999 else c1 <- newround(c1,4)
    if(is.na(c2)) c2 <- 999999 else c2 <- newround(c2,4)
    if(is.na(omega)) omega <- 999999 else omega <- newround(omega,4)
    if(is.na(z)) z <- 999999 else z <- newround(z,4)
    .space <- c(13,12,12,14)
    .length <- c(nchar(c1),nchar(c2),nchar(omega),nchar(z))
    .pad <- .space - .length
    c1.pad <- c2s(rep(' ',.pad[1]),sep='')
    c2.pad <- c2s(rep(' ',.pad[2]),sep='')
    omega.pad <- c2s(rep(' ',.pad[3]),sep='')
    z.pad <- c2s(rep(' ',.pad[4]),sep='')
    c2s(c(c1.pad,c1,c2.pad,c2,omega.pad,omega,z.pad,z),sep='')
  }
  MK <- function(o) {
    a <- newround(o$a1.a,6)
    b <- newround(o$a2.b,6)
    c <- newround(o$a3.c,6)
    if(is.na(a)) a <- 999999
    if(is.na(b)) b <- 999999
    if(is.na(c)) c <- 999999
    .space <- c(17,14,14)
    .length <- c(nchar(a),nchar(b),nchar(c))
    .pad <- .space - .length
    a.pad <- c2s(rep(' ',.pad[1]),sep='')
    b.pad <- c2s(rep(' ',.pad[2]),sep='')
    c.pad <- c2s(rep(' ',.pad[3]),sep='')
    c2s(c(a.pad,a,b.pad,b,c.pad,c),sep='')
  }
  T <- function(o) {
    t <- thermo$obigt$z.T[(as.numeric(rownames(o)))]
    if(is.na(t)) t <- 1000
    t <- newround(t,2)
    .space <- c(9)
    .length <- c(nchar(t))
    .pad <- .space - .length
    t.pad <- c2s(rep(' ',.pad[1]),sep='')
    c2s(c(t.pad,t),sep='')
  }
  # to get the phase transition properties
  HVD <- function(o) {
    T <- thermo$obigt$z.T[(as.numeric(rownames(o))-1)]
    s1 <- subcrt(as.character(o$name),as.character(obigt$state[as.numeric(rownames(o))-1]),T=T,P=1,convert=FALSE)
    s2 <- subcrt(as.character(o$name),as.character(o$state),T=T,P=1,convert=FALSE)
    if(!is.na(o$H)) {
      H <- newround(s2[[2]][[1]]$H - s1[[2]][[1]]$H,1)
    } else H <- 999999
    VV <- s2[[2]][[1]]$V - s1[[2]][[1]]$V
    V <- newround(VV,3)
    if(VV==0) V <- 999999
    if(!V==999999) {
      S <- s2[[2]][[1]]$S - s1[[2]][[1]]$S
      #D <- newround(S/VV,1)
      print(rownames(o))
      D <- newround(dPdTtr(as.numeric(rownames(o))-1),1)
      # for e.g. ca-phillipsite dP.dT is small enough to
      # be essentially zero but cause problems for formatting
      if(D < 0.01) D <- 0
    }
    else D <- 999999
    .space <- c(10,12,12)
    .length <- c(nchar(H),nchar(V),nchar(D))
    .pad <- .space - .length
    H.pad <- c2s(rep(' ',.pad[1]),sep='')
    V.pad <- c2s(rep(' ',.pad[2]),sep='')
    D.pad <- c2s(rep(' ',.pad[3]),sep='')
    c2s(c(H.pad,H,V.pad,V,D.pad,D),sep='')
  }
  

  # find out what we're dealing with
  i <- 1:nrow(obigt)
  iwater <- i[obigt$state[i]=='liq' & obigt$name[i]=='water']
  if(length(iwater)>0) {
    cat('write.supcrt:excluding water\n')
    i <- i[!i %in% iwater]
  }
  ielec <- i[obigt$state[i]=='aq' & obigt$name[i]=='e-']
  if(length(iwater)>0) {
    cat('write.supcrt:excluding electron\n')
    i <- i[!i %in% ielec]
  }
  # we need to drop species with Am or Cm
  # for a successful import to OBIGT
  notelement <- grep('Am|Cm',as.character(obigt$formula[i]))
  if(length(notelement)>0) {
    cat(paste('write.supcrt: excluding',c2s(as.character(obigt$name[i[notelement]])),'\n'))
    i[notelement] <- NA
    i <- i[!is.na(i)]
  }
  # drop species with parameters beyond the maier-kelley
  d <- obigt$a4.d[i]; e <- obigt$c1.e[i]; f <- obigt$c2.f[i]
  d[is.na(d)] <- 0; e[is.na(e)] <- 0; f[is.na(f)] <- 0
  hasd <- obigt$state[i]!='aq' & d != 0
  hase <- obigt$state[i]!='aq' & e != 0
  hasf <- obigt$state[i]!='aq' & f != 0
  hasstuff <- hasd | hase | hasf
  if(TRUE %in% hasstuff) {
    cat(paste('write.supcrt: excluding',c2s(obigt$name[i[hasstuff]]),'\n'))
    i <- i[!hasstuff]
  }
  imin3 <- i[obigt$state[i]=='cr4']
  if(length(imin3)>0) {
    iminin3 <- c(imin3,imin3-1,imin3-2,imin3-3)
    i <- i[!i %in% iminin3]
  }
  imin2 <- i[obigt$state[i]=='cr3']
  if(length(imin2)>0) {
    iminin2 <- c(imin2,imin2-1,imin2-2)
    i <- i[!i %in% iminin2]
  }
  imin1 <- i[obigt$state[i]=='cr2']
  if(length(imin1)>0) {
    iminin1 <- c(imin1,imin1-1)
    i <- i[!i %in% iminin1]
  }
  imin0 <- i[obigt$state[i] %in% c('cr','liq')]
  if(length(imin0)>0) {
    i <- i[!i %in% imin0]
  }
  igas <- i[obigt$state[i]=='gas']
  if(length(igas)>0) {
    i <- i[!i %in% igas]
  }
  iaq <- i[obigt$state[i]=='aq']
  if(length(iaq)>0) {
    i <- i[!i %in% iaq]
  }

  # exclude non-published data
  dropsource <- c("CHNOSZ")
  # drop aqueous species from particular sources
  # (used to make the result SUPCRT92-friendly
  # presumably by making the file shorter)
  #dropsource <- c(dropsource,'SK93','SK95','PSK99','HSS95','SS98','SSW+97')
  idrop <- iaq[obigt$ref1[iaq] %in% dropsource]
  if(length(idrop)>0) {
    iaq <- iaq[!iaq %in% idrop]
    print(paste('dropping these species from',dropsource))
    print(idrop)
  }

  if(length(i)>0) warning(paste('unused states for species',c2s(i)))
  print(data.frame(type=c('imin0','imin1','imin2','imin3','igas','iaq'),count=c(length(imin0),length(imin1),length(imin2),length(imin3),length(igas),length(iaq))))

  # usefull for testing
  #imin0 <- imin0[1]
  #imin1 <- imin1[1]
  #imin2 <- imin2[1]
  #imin3 <- imin3[1]
  #igas <- igas[1]
  #iaq <- iaq[1]

  # header: identify the file (change this as appropriate)
  write(paste('***',file,'of',date()),file)
  # header: the sources
  s1 <- as.character(obigt$ref1[c(iaq,igas,imin0,imin1,imin2,imin3)])
  s2 <- as.character(obigt$ref2[c(iaq,igas,imin0,imin1,imin2,imin3)])
  s2 <- s2[!is.na(s2)]
  sources <- sort(unique(c(s1,s2)))
  sources <- sources[!sources=='CHNOSZ']
  write('*** references follow:',file,append=TRUE)
  for(i in 1:length(sources)) {
    write(paste('*',sources[i],thermo$source$reference[match(sources[i],thermo$source$source)]),file,append=TRUE)
  }
  write('*** journal abbreviations follow:',file,append=TRUE)
  iabbrv <- grep('_',thermo$source$source)
  for(i in 1:length(iabbrv)) {
    abbrv <- thermo$source$source[iabbrv[i]]
    abbrv <- substr(abbrv,2,nchar(abbrv))
    write(paste('*',abbrv,thermo$source$reference[i]),file,append=TRUE)
  }
  write('*** end of comments, start thermodynamic data blocks',file,append=TRUE)
  # the number of comment lines
  icomm <- 5 + length(sources) + length(iabbrv)

  # minerals that do not undergo phase transitions
  write(separator,file,append=TRUE)
  write(paste(spaces(10),min0text,sep=''),file,append=TRUE)
  write(separator,file,append=TRUE)
  if(length(imin0)>0)
  for(i in 1:length(imin0)) {
    print(paste('imin0',i,obigt$name[imin0[i]]))
    write(paste('',NAME(obigt[imin0[i],])),file,append=TRUE) 
    write(paste('',ABBRV(obigt[imin0[i],])),file,append=TRUE)
    write(paste('',DATE(obigt[imin0[i],])),file,append=TRUE)
    write(paste('',GHSV(obigt[imin0[i],])),file,append=TRUE)
    write(paste('',MK(obigt[imin0[i],])),file,append=TRUE)
    write(paste('     ',T(obigt[imin0[i],])),file,append=TRUE)
  } else imin0 <- 0
  # minerals that undergo one phase transition
  write(separator,file,append=TRUE)
  write(paste(spaces(10),min1text,sep=''),file,append=TRUE)
  write(separator,file,append=TRUE)
  if(length(imin1)>0)
  for(i in 1:length(imin1)) {
    print(paste('imin1',i,obigt$name[imin1[i]-1]))
    #print('NAME')
    write(paste('',NAME(obigt[(imin1[i]-1),])),file,append=TRUE) 
    #print('ABBRV')
    write(paste('',ABBRV(obigt[(imin1[i]-1),])),file,append=TRUE)
    #print('DATE')
    write(paste('',DATE(obigt[(imin1[i]-1),])),file,append=TRUE)
    #print('GHSV')
    write(paste('',GHSV(obigt[(imin1[i]-1),])),file,append=TRUE)
    #print('MK')
    write(paste(' ',MK(obigt[(imin1[i]-1),]),T(obigt[(imin1[i]-1),]),HVD(obigt[imin1[i],]),sep=''),file,append=TRUE)
    #print('MK2')
    write(paste('',MK(obigt[imin1[i],])),file,append=TRUE)
    #print('T')
    write(paste('     ',T(obigt[imin1[i],])),file,append=TRUE)
  } else imin1 <- 0
  # minerals that undergo two phase transitions
  write(separator,file,append=TRUE)
  write(paste(spaces(10),min2text,sep=''),file,append=TRUE)
  write(separator,file,append=TRUE)
  if(length(imin2)>0)
  for(i in 1:length(imin2)) {
    print(paste('imin2',i,obigt$name[imin2[i]-2]))
    write(paste('',NAME(obigt[(imin2[i]-2),])),file,append=TRUE) 
    write(paste('',ABBRV(obigt[(imin2[i]-2),])),file,append=TRUE)
    write(paste('',DATE(obigt[(imin2[i]-2),])),file,append=TRUE)
    write(paste('',GHSV(obigt[(imin2[i]-2),])),file,append=TRUE)
    write(paste(' ',MK(obigt[(imin2[i]-2),]),T(obigt[(imin2[i]-2),]),HVD(obigt[(imin2[i]-1),]),sep=''),file,append=TRUE)
    write(paste(' ',MK(obigt[(imin2[i]-1),]),T(obigt[(imin2[i]-1),]),HVD(obigt[imin2[i],]),sep=''),file,append=TRUE)
    write(paste('',MK(obigt[imin2[i],])),file,append=TRUE)
    write(paste('     ',T(obigt[imin2[i],])),file,append=TRUE)
  } else imin2 <- 0
  # minerals that undergo three phase transitions
  write(separator,file,append=TRUE)
  write(paste(spaces(10),min3text,sep=''),file,append=TRUE)
  write(separator,file,append=TRUE)
  if(length(imin3)>0)
  for(i in 1:length(imin3)) {
    print(paste('imin3',i,obigt$name[imin3[i]-3]))
    write(paste('',NAME(obigt[(imin3[i]-3),])),file,append=TRUE) 
    write(paste('',ABBRV(obigt[(imin3[i]-3),])),file,append=TRUE)
    write(paste('',DATE(obigt[(imin3[i]-3),])),file,append=TRUE)
    write(paste('',GHSV(obigt[(imin3[i]-3),])),file,append=TRUE)
    write(paste(' ',MK(obigt[(imin3[i]-3),]),T(obigt[(imin3[i]-3),]),HVD(obigt[(imin3[i]-2),]),sep=''),file,append=TRUE)
    write(paste(' ',MK(obigt[(imin3[i]-2),]),T(obigt[(imin3[i]-2),]),HVD(obigt[(imin3[i]-1),]),sep=''),file,append=TRUE)
    write(paste(' ',MK(obigt[(imin3[i]-1),]),T(obigt[(imin3[i]-1),]),HVD(obigt[imin3[i],]),sep=''),file,append=TRUE)
    write(paste('',MK(obigt[imin3[i],])),file,append=TRUE)
    write(paste('     ',T(obigt[imin3[i],])),file,append=TRUE)
  } else imin3 <- 0
  # gases
  write(separator,file,append=TRUE)
  write(paste(spaces(10),gastext,sep=''),file,append=TRUE)
  write(separator,file,append=TRUE)
  if(length(igas)>0)
  for(i in 1:length(igas)) {
    print(paste('igas',i,obigt$name[igas[i]]))
    write(paste('',NAME(obigt[igas[i],])),file,append=TRUE) 
    write(paste('',ABBRV(obigt[igas[i],])),file,append=TRUE)
    write(paste('',DATE(obigt[igas[i],])),file,append=TRUE)
    write(paste('',GHSV(obigt[igas[i],])),file,append=TRUE)
    write(paste('',MK(obigt[igas[i],])),file,append=TRUE)
    write(paste('     ',T(obigt[igas[i],])),file,append=TRUE)
  } else igas <- 0
  # aqueous species
  write(separator,file,append=TRUE)
  write(paste(spaces(10),aqtext,sep=''),file,append=TRUE)
  write(separator,file,append=TRUE)
  if(length(iaq)>0)
  for(i in 1:length(iaq)) {
    print(paste('iaq',i,obigt$name[iaq[i]]))
    write(paste('',NAME(obigt[iaq[i],])),file,append=TRUE) 
    write(paste('',ABBRV(obigt[iaq[i],])),file,append=TRUE)
    write(paste('',DATE(obigt[iaq[i],])),file,append=TRUE)
    write(paste('',GHS(obigt[iaq[i],])),file,append=TRUE)
    write(paste('',HKF1(obigt[iaq[i],])),file,append=TRUE)
    write(paste('',HKF2(obigt[iaq[i],])),file,append=TRUE)
  } else iaq <- 0

  # footer
  write('\n',file,append=TRUE)
  write(c2s(c(spaces(9-nchar(length(imin0))),length(imin0),'  ',min0textalt,'  (nmin1)'),sep=''),file,append=TRUE)
  write(c2s(c(spaces(9-nchar(length(imin1))),length(imin1),'  ',min1text,'     (nmin2)'),sep=''),file,append=TRUE)
  write(c2s(c(spaces(9-nchar(length(imin2))),length(imin2),'  ',min2text,'    (nmin3)'),sep=''),file,append=TRUE)
  write(c2s(c(spaces(9-nchar(length(imin3))),length(imin3),'  ',min3text,'  (nmin4)'),sep=''),file,append=TRUE)
  write(c2s(c(spaces(9-nchar(length(igas))),length(igas),'  ',gastext,'                                          (ngas)'),sep=''),file,append=TRUE)
  write(c2s(c(spaces(9-nchar(length(iaq))),length(iaq),'  ',aqtext,'                                (naqs)'),sep=''),file,append=TRUE)
  write('     ----  ----------------------------',file,append=TRUE)
  total <- sum(length(imin0),length(imin1),length(imin2),length(imin3),length(igas),length(iaq))
  write(c2s(c(spaces(9-nchar(total)),total,'  ','total'),sep=''),file,append=TRUE)
  write('',file,append=TRUE)
  write(c2s(c(spaces(9-nchar(icomm)),icomm,'  ',commtext),sep=''),file,append=TRUE)
  # convert to dos line endings
  #if(grep('linux',R.version$platform)) system(paste('unix2dos',file))

  return()
}

