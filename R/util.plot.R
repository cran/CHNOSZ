# CHNOSZ/util.plot.R
# Copyright (C) 2006-2009 Jeffrey M. Dick
# Functions to create and modify plots

thermo.plot.new <- function(xlim,ylim,xlab,ylab,cex=par('cex'),mar=NULL,lwd=par('lwd'),side=c(1,2,3,4),
  mgp=c(1.2,0.3,0),cex.axis=par('cex'),col=par('col'),yline=NULL,axs='i',do.box=TRUE,ticks=NULL,
  las=1,xline=NULL) {
  # start a new plot with some customized settings
  # 20091108 changed argument name from 'ticks' to 'side' but
  # keep 'ticks' for backward compatibility
  if(!is.null(ticks)) side <- ticks 
  # 20090324 mar handling: NULL - a default setting; NA - par's setting
  # 20090413 changed mar of top side from 2 to 2.5
  if(is.null(mar)) mar <- c(3,3.5,2.5,1) else if(is.na(mar[1])) mar <- par('mar')
  par(mar=mar,mgp=mgp,tcl=0.3,las=las,xaxs=axs,yaxs=axs,cex=cex,lwd=lwd)
  plot.new()
  plot.window(xlim=xlim,ylim=ylim)
  if(do.box) box()
  # labels
  if(is.null(xline)) xline <- mgp[1]
  thermo.axis(xlab,side=1,line=xline,cex=cex.axis,lwd=NULL)
  if(is.null(yline)) yline <- mgp[1]
  thermo.axis(ylab,side=2,line=yline,cex=cex.axis,lwd=NULL)
  # (optional) tick marks
  if(1 %in% side) thermo.axis(NULL,side=1,lwd=lwd,col=par('col'))
  if(2 %in% side) thermo.axis(NULL,side=2,lwd=lwd,col=par('col'))
  if(3 %in% side) thermo.axis(NULL,side=3,lwd=lwd,col=par('col'))
  if(4 %in% side) thermo.axis(NULL,side=4,lwd=lwd,col=par('col'))
}

thermo.postscript <- function(file,family='Helvetica',width=8,height=6,horizontal=FALSE) {
  postscript(onefile=FALSE,horizontal=horizontal,paper='special',width=width,height=height,file=file,family=family)
}

thermo.axis <- function(lab='x-axis',side=1,line=1.5,cex=par('cex'),lwd=par('lwd'),T=NULL,col=par('col')) {
  # if T isn't NULL, looks like we want make a second
  # oxidation scale corresponding to one already plotted.
  # e.g.,  Eh-pe, Eh-logfO2, or logfO2-Eh
  if(!is.null(T)) {
    usr <- par('usr')
    if(side %in% c(1,3)) lim <- usr[1:2] else lim <- usr[3:4]
    if(length(grep('pe',lab)) > 0) {
      lim <- convert(lim,'pe',T=T)
    } else if(length(grep('O2',lab)) > 0) {
      lim <- convert(lim,'logfO2',T=T)
    } else if(length(grep('Eh',lab)) > 0) {
      lim <- convert(lim,'E0',T=T)
    }
    if(side %in% c(1,3)) usr[1:2] <- lim else usr[3:4] <- lim
    opar <- par(usr=usr)
  }
  if(!is.null(lwd)) {
    ## plot major tick marks and numeric labels
    do.label <- TRUE
    if(missing(cex) & side %in% c(3,4) & is.null(T)) do.label <- FALSE
    at <- axis(side,labels=do.label,tick=TRUE,lwd=lwd,col=col,col.axis=col) 
    ## plot minor tick marks
    # the distance between major tick marks
    da <- abs(diff(at[1:2]))
    # distance between minor tick marks
    di <- da / 4
    if(da %% 2 | !(da %% 10)) di <- da / 5
    # number of minor tick marks
    if(side %in% c(1,3)) {
      ii <- c(1,2) 
      myasp <- par('xaxp')
    } else {
      ii <- c(3,4)
      myasp <- par('yaxp')
    }
    myusr <- par('usr')[ii]
    daxis <- abs(diff(myusr))
    nt <- daxis / di + 1
    ## if nt isn't an integer, it probably
    ## means the axis limits don't correspond
    ## to major tick marks (expect problems)
    ##at <- seq(myusr[1],myusr[2],length.out=nt)
    # start from (bottom/left) of axis?
    bl <- 1
    #if(myasp[2]==myusr[2]) bl <- 2
    # is forward direction (top/right)?
    tr <- 1
    if(xor(myusr[2] < myusr[1] , bl==2)) tr <- -1
    #at <- myusr[bl] + tr * di * seq(0:(nt-1))
    # well all of that doesn't work in a lot of cases,
    # where none of the axis limits correspond to
    # major tick marks. perhaps the following will work
    at <- myusr[1] + tr * di * (0:(nt-1))
    # apply an offset
    axt <- axTicks(side)[1]
    daxt <- (axt - myusr[1])/di
    daxt <- (daxt-round(daxt))*di
    at <- at + daxt
    tcl <- par('tcl') * 0.5
    axis(side,labels=FALSE,tick=TRUE,lwd=lwd,col=col,col.axis=col,at=at,tcl=tcl)
  }

  # rotate labels on side axes
  if(side %in% c(2,4)) las <- 0 else las <- 1
  if(!is.null(lab)) mtext(lab,side=side,line=line,cex=cex,las=las)
  # reset limits if we were plotting a second axis
  if(!is.null(T)) par(opar)
}

label.plot <- function(x,xfrac=0.95,yfrac=0.9,cex=1,paren=TRUE,adj=1) {
  # make a text label e.g., "(a)" in the corner of a plot
  # xfrac, yfrac: fraction of axis where to put label (default top right)
  # paren: put a parenthesis around the text, and italicize it?
  if(paren) x <- as.expression(substitute(group('(',italic(a),')'),list(a=x)))
  pu <- par('usr')
  text( pu[1]+xfrac*(pu[2]-pu[1]), pu[3]+yfrac*(pu[4]-pu[3]), labels=x, cex=cex , adj=adj)
}

axis.label <- function(lab,opt=NULL,do.state=TRUE,oldstyle=FALSE,do.upper=FALSE,mol='mol') {
  # make axis labels
  # 20090826: just return the argument if a comma is already present
  if(length(grep(",",lab)) > 0) return(lab)
  if(missing(opt)) do.opt <- TRUE else do.opt <- FALSE
  if(!is.null(opt)) if(is.na(opt)) do.opt <- TRUE
  if(lab %in% c('T','P','Eh','pH','pe','logK','IS')) {
    # the label is one of these properties
    if(lab=='Eh') lab <- paste(lab,'(volt)')
    else if(lab=='T') {
      if(do.opt) T.units <- nuts('T') else T.units <- opt
      if(T.units=='C') lab <- as.expression(quote(list(italic(T),degree*C)))
      else lab <- as.expression(quote(list(italic(T),K)))
    } 
    else if(lab=='P') {
      if(do.opt) P.units <- nuts('P') else P.units <- opt
      if(P.units=='bar') lab <- as.expression(quote(list(italic(P),bar)))
      else lab <- as.expression(quote(list(italic(P),MPa)))
    } 
    else if(lab=='logK') lab <- as.expression(quote(log~italic(K)))
    else if(lab=='IS') lab <- as.expression(quote(list(IS,mol~~kg^{-1})))
    return(lab)
  } else {
    # the label is a chemical activity or fugacity
    if(is.null(thermo$basis)) rn <- '' else rn <- rownames(basis())
    if(lab %in% rn) {
      # 20090215: the state this basis species is in
      if(do.opt) opt <- as.character(thermo$basis$state)[match(lab,rn)]
      state <- opt
      if(oldstyle) {
        # append (log a) or (log f)
        if(state %in% c('gas')) llab <- '(log f)' else llab <- '(log a)'
        newlab <- paste(lab,llab,sep=' ')
        return(newlab)
      } else {
        return(species.label(lab,do.state=do.state,state=state,do.log=TRUE))
      }
    } else {
      # a way to make expressions for various properties
      # e.g. axis.label('DG0r','k') for standard molal Gibbs energy
      # change of reaction in kcal/mol
      clab <- s2c(lab)
      wlab <- ''
      doital <- TRUE; dosub <- FALSE
      for(i in 1:length(clab)) {
        blab <- clab[i]
        #if(dosub) blab <- substitute(phantom0[a],list(a=blab))
        # D for Delta
        if(i==1 & blab=='D') wlab <- quote(Delta)
        else if(i > 1 & blab=='0') {
          wlab <- substitute(a*degree,list(a=wlab))
          dosub <- TRUE
        }
        else if(i > 1 & (can.be.numeric(blab) | blab=='P' | dosub)) 
          wlab <- substitute(a[italic(b)],list(a=wlab,b=blab))
        else if(i > 1 & blab==',') {
          wlab <- substitute(a*b,list(a=wlab,b=blab))
          doital <- FALSE
        }
        else {
          if(blab=='A') blab <- substitute(bold(a),list(a=blab))
          if(i > 1) {
            if(blab=='p') wlab <- substitute(a[italic(b)],list(a=wlab,b=toupper(blab)))
            else if(doital) wlab <- substitute(a*italic(b),list(a=wlab,b=blab))
            else wlab <- substitute(a*b,list(a=wlab,b=blab))
          } else wlab <- substitute(italic(b),list(b=blab))
        }
      }
      # now to the nuts
      if(do.state) {
        mytoupper <- function(lab) {
          if(do.upper & lab!='g') return(toupper(lab))
          else return(lab)
        }
        if(clab[1]=='D') clab <- clab[-1]
        ulab <- mytoupper(nuts('E'))
        if(clab[1] %in% c('C','S')) mylab <- substitute(a~~K^-1,list(a=ulab))
        else if(clab[1] == c('V')) mylab <- substitute(a^3,list(a=mytoupper('cm')))
        else if(clab[1] == 'E') mylab <- substitute(a^3~~K^-1,list(a=mytoupper('cm')))
        else mylab <- ulab
        if(!is.null(opt)) {
          if(can.be.numeric(opt)) mylab <- substitute(10^a~~b,list(a=opt,b=mylab))
          else {
            opt <- mytoupper(opt)
            mylab <- substitute(a*b,list(a=opt,b=mylab))
          }
        }
        mylab <- substitute(a~~b^-1,list(a=mylab,b=mytoupper(mol)))
        wlab <- substitute(list(a,b),list(a=wlab,b=mylab))
      }
      return(as.expression(wlab))
    }
  }
}

species.label <- function(formula,do.state=FALSE,state="",do.log=FALSE) {
  # make plotting expressions for chemical formulas
  # that include subscripts, superscripts (if charged)
  # and optionally designations of states +/- loga or logf prefix
  labform <- makeup(formula)
  newlab <- ''
  for(i in 1:nrow(labform)) {
      if(rownames(labform)[i] != 'Z') {
        newlab <- substitute(paste(a,b),list(a=newlab,b=rownames(labform)[i]))
        if(labform$count[i]!=1) {
          # subscripts within subscripts (do.log) are too small
          if(do.log) newlab <- substitute(a*b,list(a=newlab,b=labform$count[i]))
          else newlab <- substitute(a[b],list(a=newlab,b=labform$count[i]))
        }
      } else {
        # for charged species, don't show "Z" but do show e.g. "+2"
        lc <- labform$count[i]
        if(lc==-1) lc <- "-"
        else if(lc==1) lc <- "+"
        else if(lc > 0) lc <- paste("+",as.character(lc),sep="")
        if(do.log) newlab <- substitute(paste(a,b),list(a=newlab,b=lc))
        else newlab <- substitute(a^b,list(a=newlab,b=lc))
      }
  }
  if(do.state) {
    llab <- 'f'; slab <- 'g'
    if(state %in% c('aq','cr','liq')) {llab <- 'a'; slab <- state}
    #newlab <- substitute(a[group('(',italic(b),')')],list(a=newlab,b=slab))
    newlab <- substitute(a*group('(',italic(b),')'),list(a=newlab,b=slab))
    if(do.log) {
      llab <- substitute(log*italic(a),list(a=llab))
      newlab <- substitute(a[b],list(a=llab,b=newlab))
    }
  }
  return(as.expression(newlab))
}


water.lines <- function(xaxis='pH',yaxis='Eh',T=298.15,P='Psat',which=c('oxidation','reduction'),logaH2O=0,lty=2,col=par('fg'),xpoints=NULL) {
  # draw water stability limits
  # if we're on an Eh-pH diagram, or logfO2-pH diagram,
  # or logfO2-T or Eh-T
  # calculate them exactly (nicer looking lines), otherwise 
  # (TODO) add them using affinity() and diagram()
  
  # get the x and y limits from the current plot
  pu <- par('usr')
  xlim <- pu[1:2]
  ylim <- pu[3:4]
  # exact lines
  # warning: Eh calculations are reliable only at a single T
  if(xaxis=='pH' & (yaxis=='Eh' | yaxis=='O2' | yaxis=="pe")) {
    if('reduction' %in% which) {
      logfH2 <- 0
      logK <- subcrt(c('H2O','oxygen','hydrogen'),c(-1,0.5,1),T=T,P=P,convert=FALSE)$out$logK 
      logfO2 <- 2 * logK - logfH2 + 2 * logaH2O
      if(yaxis=='O2') abline(h=logfO2,lty=lty,col=col) 
      else if(yaxis=="Eh") lines(xlim,convert(logfO2,'E0',T=T,P=P,pH=xlim),lty=lty,col=col)
      else if(yaxis=="pe") lines(xlim,convert(convert(logfO2,'E0',T=T,P=P,pH=xlim),"pe",T=T),lty=lty,col=col)
    }
    if('oxidation' %in% which) {
      logfO2 <- 0
      if(yaxis=='O2') abline(h=logfO2,lty=lty,col=col) 
      else if(yaxis=="Eh") lines(xlim,convert(logfO2,'E0',T=T,P=P,pH=xlim),lty=lty,col=col)
      else if(yaxis=="pe") lines(xlim,convert(convert(logfO2,'E0',T=T,P=P,pH=xlim),"pe",T=T),lty=lty,col=col)
    }
  } else if(xaxis %in% c('T','P') & yaxis %in% c('Eh','O2') ) {
    #if(xaxis=='T') if(is.null(xpoints)) xpoints <- T
    # 20090212 get T values from plot limits
    # TODO: make this work for T on y-axis too
    if(xaxis=='T') {
      if(missing(T)) {
        xpoints <- seq(xlim[1],xlim[2],length.out=100)
        T <- envert(xpoints,"K")
      }
    }
    if(xaxis=='P') if(is.null(xpoints)) xpoints <- P
    if('oxidation' %in% which) {
      logfO2 <- rep(0,length(xpoints))
      if(yaxis=='Eh') lines(xpoints,convert(logfO2,'E0',T=T,P=P,pH=xlim),lty=lty,col=col)
      else lines(xpoints,logfO2,lty=lty,col=col)
    }
    if('reduction' %in% which) {
      logfH2 <- 0
      logK <- subcrt(c('H2O','oxygen','hydrogen'),c(-1,0.5,1),T=T,P=P,convert=FALSE)$out$logK 
      logfO2 <- 2 * logK - logfH2 + 2 * logaH2O
      if(yaxis=='Eh') lines(xpoints,convert(logfO2,'E0',T=T,P=P,pH=xlim),lty=lty,col=col)
      else lines(xpoints,logfO2,lty=lty,col=col)
    }
  } else {
    # inexact lines
    #
  }
}

mtitle <- function(main,...) {
  # make a possibly multi-line plot title 
  # useful for including expressions on multiple lines 
  l <- length(main)
  for(i in 1:l) mtext(main[i],line=l-i,...)
}

