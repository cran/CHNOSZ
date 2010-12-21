# CHNOSZ/affinity.R
# calculate affinities of formation reactions

energy <- function(what,vars,vals,lims,T=thermo$opt$Tr,P="Psat",IS=0,sout=NULL,do.phases=FALSE,transect=FALSE) {
  # 20090329 extracted from affinity() and made to
  # deal with >2 dimensions (variables)

  # calculate "what" property
  # logK - logK of formation reactions
  # logact.basis - logacts of basis species in reactions
  # logact.species - logacts of species in reactions
  # logQ - logQ of formation reactions
  # A - A/2.303RT of formation reactions
  # "vars": vector of variables (T,P,basisnames,IS)
  # "vals": list of values for each variable
  # "lims": used to get the dimensions (resolution for each variable)
  # (TODO: remove G.basis from documentation)

  ### some argument checking
  if(length(unique(vars)) != length(vars)) stop("please supply unique variable names")
  ## basis definition / number of basis species 
  mybasis <- thermo$basis  # FIXME: use basis() here
  nbasis <- nrow(mybasis)
  ## species definition / number of species
  myspecies <- species(quiet=TRUE)
  if(is.character(what)) {
    if(is.null(myspecies)) stop('species properties requested, but species have not been defined')
    nspecies <- nrow(myspecies)
    if(!identical(what,"logact.basis")) ispecies <- 1:nspecies
  }
  ## the dimensions of our arrays
  resfun <- function(lim) lim[3]
  mydim <- sapply(lims,resfun)
  if(transect) {
    if(!isTRUE(all.equal(min(mydim),max(mydim)))) stop("variables define a transect but their lengths are not all equal")
    mydim <- mydim[[1]]
  }
  ## the number of dimensions we have
  if(transect) nd <- 1 else nd <- length(vars) # == length(mydim)
  ## basis names / which vars denote basis species
  basisnames <- rownames(mybasis)
  ibasisvar <- match(vars,basisnames)
  varisbasis <- !is.na(ibasisvar)
  ibasisvar <- ibasisvar[!is.na(ibasisvar)]
  ## which vars are in P,T,IS
  varissubcrt <- vars %in% c("P","T","IS")
  if(length(which(varissubcrt)) > 2) stop("sorry, currently only up to 2 of P,T,IS are supported")
  ## categorize the basis species:
  # 0 - not in the vars; 1 - one of the vars
  ibasis <- 1:nbasis
  ibasis0 <- ibasis[!ibasis %in% ibasisvar]
  ibasis1 <- ibasis[ibasis %in% ibasisvar]
  if(identical(what,"logact.basis")) ispecies <- ibasis
  ## what subcrt variable is used to make a 2-D grid?
  if(length(which(varissubcrt)) > 1 & !transect) {
    if("IS" %in% vars) grid <- "IS"
    else grid <- vars[varissubcrt][1]
  } else grid <- NULL
  ### done argument processing

  ### function to index the variables in a permutated order
  # by swapping ivar for the first
  # e.g. for nd=5, ivars(4)=c(4,2,3,1,5)
  ivars <- function(ivar,iv=NULL) {
    if(nd==0) return(1)
    if(is.null(iv)) iv <- 1:nd
    iv.1 <- iv[1]
    iv[1] <- ivar
    iv[ivar] <- iv.1
    return(iv)
  }

  ### functions for logact / logQ
  # a basis species not in var
  logact.basis0.fun <- function(ibasis) {
    logact <- mybasis$logact[ibasis]
    # for the case where this basis species is buffered
    if(!can.be.numeric(logact)) logact <- 0
    else logact <- as.numeric(logact)
    return(array(logact,mydim))
  }
  # a basis species in var
  logact.basis1.fun <- function(ivar) {
    dim.in <- dim(vals[[ivar]])
    if(is.null(dim.in)) dim.in <- 0
    # why doesn't this work?
    # if(identical(dim.in,mydim))
    if(all(dim.in==mydim) | transect) return(vals[[ivar]])
    else return(aperm(array(vals[[ivar]],mydim[ivars(ivar)]),ivars(ivar)))
  }
  # any basis species
  logact.basis.fun <- function(ibasis) {
    if(ibasis %in% ibasis0) return(logact.basis0.fun(ibasis))
    else return(logact.basis1.fun(match(basisnames[ibasis],vars)))
  }
  # all basis species
  logact.basis <- function() lapply(ibasis,logact.basis.fun)
  # logact of a single species
  logact.species.fun <- function(ispecies) array(myspecies$logact[ispecies],mydim)
  ## contributions to logQ
  # from a single basis species 
  logQ.basis.fun <- function(ibasis,coeff) - coeff * logact.basis.fun(ibasis)
  # from all basis species in a single formation reaction
  logQ.basis.species <- function(ispecies) 
    psum(mapply(logQ.basis.fun,ibasis,myspecies[ispecies,1:nbasis],SIMPLIFY=FALSE))
  # all basis species in all reactions
  logQ.basis <- function() mapply(logQ.basis.species,1:nspecies,SIMPLIFY=FALSE)
  # by a single species
  logQ.species.fun <- function(ispecies,coeff) coeff * logact.species.fun(ispecies)
  # by all species
  logQ.species <- function() 
    mapply(logQ.species.fun,1:nspecies,1,SIMPLIFY=FALSE)
  # total logQ of all reactions
  logQ <- function() lsum(logQ.basis(),logQ.species())

  ### function for calling subcrt
  sout.fun <- function(property="logK") {
    if(!is.null(sout)) return(sout) else {
      ## subcrt arguments
      species <- c(mybasis$ispecies,myspecies$ispecies)
      if("T" %in% vars) T <- vals[[which(vars=="T")]]
      if("P" %in% vars) P <- vals[[which(vars=="P")]]
      s.args <- list(species=species,property=property,T=T,P=P,grid=grid,convert=FALSE,do.phases=do.phases)
      if("IS" %in% vars) {
        IS <- vals[[which(vars=="IS")]]
        # do the calculation in parts: basis species (IS=0)
        s.args$species <- mybasis$ispecies
        s.out.basis <- do.call("subcrt",s.args)$out
        # and species of interest (IS=IS)
        s.args$species <- myspecies$ispecies
        s.args <- c(s.args,list(IS=IS))
        s.out.species <- do.call("subcrt",s.args)$out
        # put them together
        return(c(s.out.basis,s.out.species))
      } else return(do.call("subcrt",s.args)$out)
    }
  }

  ### functions for logK/subcrt props
  # the logK contribution by any species or basis species
  X.species <- function(ispecies,coeff,X) coeff * sout[[ispecies]][,names(sout[[ispecies]])==X]
  # the logK contribution by all basis species in a reaction
  X.basis <- function(ispecies,X) psum(mapply(X.species,ibasis,-myspecies[ispecies,ibasis],X,SIMPLIFY=FALSE))
  # the logK of any reaction
  X.reaction <- function(ispecies,X) X.species((ispecies+nbasis),1,X) + X.basis(ispecies,X)
  # to get logK or subcrt props or other values into the dimensions we are using
  dim.fun <- function(x,idim=NULL) {
    if(is.null(idim)) {
      if(transect) idim <- 1
      else if(is.null(grid)) {
        # one of T,P,IS
        ivar <- which(vars %in% c("T","P","IS"))
        if(length(ivar)==0) ivar <- 1
        idim <- ivars(ivar)
      } else {
        # two of T,P,IS
        ivar1 <- which(varissubcrt)[1]
        ivar2 <- which(varissubcrt)[2]
        idim <- ivars(ivar2,ivars(ivar1))
      }
    }
    return(aperm(array(x,mydim[idim]),idim))
  }
  # properties of all species
  X.fun <- function(X) lapply(lapply(ispecies,X.reaction,X),dim.fun)
  logK <- function() lapply(ispecies,X.reaction,"logK")
  # A/2.303RT
  A <- function() lsub(X.fun("logK"),logQ())

  ### call the necessary functions
  if(!is.character(what)) {
    # expand numeric values into our dimensions
    # (used by energy.args() for calculating pe=f(Eh,T) )
    # TODO: document that sout here denotes the dimension
    # we're expanding into
    return(dim.fun(what,ivars(sout)))
  } else if(what %in% c('G','H','S','Cp','V','E','kT','logK')) {
    # get subcrt properties for reactions
    sout <- sout.fun(what)
    a <- X.fun(what)
  } else if(length(agrep("species",what)) > 0) {
    # get subcrt properties for species
    # e.g. what=G.species, Cp.species etc
    # NOTE: change from previous calling style
    # (particularly affects Cp example in ionize.Rd
    # and G.Z in protein.info function)
    mywhat <- s2c(what,sep=".",keep.sep=FALSE)[1]
    sout <- sout.fun(mywhat)
    a <- lapply(mapply(X.species,ispecies+nbasis,1,mywhat,SIMPLIFY=FALSE),dim.fun)
  } else {
    # get affinities or logact.basis, logact.species or logQ
    if(what=="A") sout <- sout.fun("logK")
    a <- do.call(what,list())
  }

  ### use species indices as the names 
  if(what=="logact.basis") names(a) <- mybasis$ispecies
  else names(a) <- myspecies$ispecies
  ### return the results
  return(list(sout=sout,a=a))
}

energy.args <- function(args,quiet=FALSE) {
  ## extracted from affinity() and modified 20090329 jmd
  # takes a list of arguments which define the dimensions
  # over which to calculate logQ, logK and affinity
  # the names should be T, P, IS and names of basis species
  # (or pH, pe, Eh)
  ## inputs are like c(T1,T2,res)
  # and outputs are like seq(T1,T2,length.out=res)
  # unless transect: do the variables specify a transect? 20090627
  transect <- any(sapply(args,length) > 3)
  ## convert T, P args and take care of 
  # single values for T, P or IS
  T <- thermo$opt$Tr
  P <- "Psat"
  IS <- 0
  T.is.var <- P.is.var <- IS.is.var <- FALSE
  arg.is.T <- names(args)=="T"
  arg.is.P <- names(args)=="P"
  arg.is.IS <- names(args)=="IS"
  if(any(arg.is.T)) {
    T <- args[[which(arg.is.T)]]
    if(length(T) > 1) T.is.var <- TRUE
    if(transect) args[[which(arg.is.T)]] <- T <- envert(T,'K')
    else args[[which(arg.is.T)]][1:2] <- T[1:2] <- envert(T[1:2],'K')
    T <- T[!is.na(T)]
  }
  if(any(arg.is.P)) {
    P <- args[[which(arg.is.P)]]
    if(length(P) > 1) P.is.var <- TRUE
    if(!identical(P,"Psat")) {
      if(transect) args[[which(arg.is.P)]] <- P <- envert(P,'bar')
      else args[[which(arg.is.P)]][1:2] <- P[1:2] <- envert(P[1:2],'bar')
    }
    P <- P[!is.na(P)]
  }
  if(any(arg.is.IS)) {
    IS <- args[[which(arg.is.IS)]]
    if(length(IS) > 1) IS.is.var <- TRUE
  }
  # report non-variables to user
  if(!T.is.var & !quiet)
    cat(paste('affinity: temperature is ',outvert(T,'K'),' ',nuts('t'),'\n',sep=''))
  if(!P.is.var & !quiet) {
    if(identical(P,"Psat")) cat("energy.args: pressure is Psat\n")
    else cat(paste('affinity: pressure is ',outvert(P,'bar'),' ',nuts('p'),'\n',sep=''))
  }
  if(!IS.is.var & !identical(IS,0) & !quiet) cat(paste('affinity: ionic strength is ',IS,'\n',sep=''))
  # default values for resolution
  res <- 128
  # where we store the output
  what <- "A"
  vars <- character()
  vals <- list(0)
  lims <- list(c(0,0,1))
  # clean out non-variables
  if(any(arg.is.T) & !T.is.var) args <- args[names(args)!="T"]
  if(any(arg.is.P) & !P.is.var) args <- args[names(args)!="P"]
  if(any(arg.is.IS) & !IS.is.var) args <- args[names(args)!="IS"]
  # the property we're interested in
  if("what" %in% names(args)) {
    what <- args[[names(args)=="what"]]
    args <- args[names(args)!="what"]
  }
  # assemble the variables
  if(length(args) > 0) {
    for(i in 1:length(args)) {
      names.orig <- names(args)[i]
      if(transect) lims.orig <- c(min(args[[i]]),max(args[[i]]))
      else lims.orig <- args[[i]][1:2]
      if(names(args)[i]=="pH") {
        names(args)[i] <- "H+"
        if(transect) args[[i]] <- -args[[i]]
        else args[[i]][1:2] <- -args[[i]][1:2]
        if(!'H+' %in% rownames(thermo$basis)) 
          warning('energy.args: pH requested, but no H+ in the basis',immediate.=TRUE)
      } 
      if(names(args)[i]=="pe") {
        names(args)[i] <- "e-"
        if(!'e-' %in% rownames(thermo$basis)) 
          warning('energy.args: pe requested, but no e- in the basis',immediate.=TRUE)
        if(transect) args[[i]] <- -args[[i]]
        else args[[i]][1:2] <- -args[[i]][1:2]
      }
      if(length(args[[i]]) < 3 & !transect) args[[i]] <- c(args[[i]],res)
      vars[length(vars)+1] <- names(args)[i]
      if(transect) {
        vals[[length(vars)]] <- args[[i]]
        lims[[length(vars)]] <- c(lims.orig,length(vals[[i]]))
      } else {
        vals[[length(vars)]] <- seq(args[[i]][1],args[[i]][2],length.out=args[[i]][3])
        lims[[length(vars)]] <- args[[i]]
      }
      names(lims)[length(vars)] <- names(args)[i]
      # say something
      if(transect) n <- length(args[[i]]) else n <- args[[i]][3]
      if(!quiet) cat(paste("energy.args: variable",length(vars),"is",names.orig,"at",
        n,"increments from",lims.orig[1],"to",lims.orig[2],"\n"))
    }
  }
  args <- list(what=what,vars=vars,vals=vals,lims=lims,T=T,P=P,IS=IS,transect=transect)

  # convert Eh to pe
  if("Eh" %in% args$vars) {
    # get Eh into our dimensions
    Eh.args <- args
    # what variable is Eh
    Eh.var <- which(args$vars=="Eh")
    Eh.args$what <- args$vals[[Eh.var]]
    Eh.args$sout <- Eh.var
    Eh <- do.call("energy",Eh.args)
    # get temperature into our dimensions
    T.args <- args  
    if("T" %in% args$vars) {
      T.var <- which(args$vars=="T")
      T.args$what <- args$vals[[T.var]]
    } else {
      T.var <- 1
      T.args$what <- T
    }
    T.args$sout <- T.var
    T <- do.call("energy",T.args)
    # do the conversion on vectors
    mydim <- dim(Eh)
    Eh <- as.vector(Eh)
    T <- as.vector(T)
    pe <- convert(Eh,"pe",T=T)
    dim(pe) <- mydim
    # update the arguments list
    args$vars[Eh.var] <- "e-"
    args$vals[[Eh.var]] <- -pe
  }
  return(args)
}

affinity <- function(...,property=NULL,sout=NULL,do.phases=FALSE,
  return.buffer=FALSE,balance="PBB",quiet=FALSE,iprotein=NULL,loga.protein=-3) {
  # ...: variables over which to calculate
  # property: what type of energy
  #   (G.basis,G.species,logact.basis,logK,logQ,A)
  # return.buffer: return buffered activities
  # balance: balance protein buffers on PBB
  # do.phases: take 999999 for the Gibbs energies
  #   of minerals beyond their T-limits?
  # sout: provide a previously calculated output from subcrt
  # quiet: make less noise
  # iprotein: build these proteins from residues (speed optimization)

  # history: 20061027 jmd version 1
  # this is where energy.args() used to sit
  # this is where energy() used to sit

  # the argument list
  args <- energy.args(list(...),quiet=quiet)
  args <- c(args,list(sout=sout,do.phases=do.phases))

  # the species we're given
  mybasis <- thermo$basis
  myspecies <- thermo$species

  if(!is.null(property)) {
    # the user just wants an energy property
    buffer <- FALSE
    args$what <- property
    out <- do.call("energy",args)
    a <- out$a
    sout <- out$sout
  } else {

    # affinity calculations
    property <- args$what

    # iprotein stuff
    if(!is.null(iprotein)) {
      # add protein residues to the species list
      resnames <- c("H2O",aminoacids(nchar=3))
      # residue activities set to zero;
      # account for protein activities later
      ires <- species(paste(resnames,"RESIDUE",sep="_"),0,quiet=TRUE)
    }

    # buffer stuff
    buffer <- FALSE
    ibufbasis <- which(!can.be.numeric(mybasis$logact))
    if(!is.null(mybasis) & length(ibufbasis) > 0) {
      buffer <- TRUE
      if(!quiet) cat('affinity: loading buffer species\n')
      if(!is.null(thermo$species)) is.species <- 1:nrow(thermo$species) else is.species <- numeric()
      is.buffer <- buffer(logK=NULL)
      is.buff <- numeric()
      for(i in 1:length(is.buffer)) is.buff <- c(is.buff,as.numeric(is.buffer[[i]]))
      is.only.buffer <- is.buff[!is.buff %in% is.species]
      buffers <- names(is.buffer)
      # reorder the buffers according to thermo$buffer
      buffers <- buffers[order(match(buffers,thermo$buffer$name))]
    }

    # ionization stuff
    ionize <- FALSE
    if( (!is.null(iprotein) | length(grep('_',as.character(thermo$species$name))) > 0) & 
        'H+' %in% rownames(mybasis) & thermo$opt$ionize) {
      ionize <- TRUE
      if(!quiet) cat('affinity: loading ionizable protein groups\n')
      is.species <- 1:nrow(thermo$species)
      is.ion <- ionize(affinity=NULL)
      is.only.ion <- is.ion[!is.ion %in% is.species]
      # keep the value of H+ if we're buffering it
      H.act <- as.character(mybasis$logact[rownames(mybasis)=='H+'])
    }

    # here we call energy
    aa <- do.call("energy",args)
    a <- aa$a
    sout <- aa$sout

    # more ionization stuff
    if(ionize) {
      #charge <- ionize(a,a)
      if(length(is.only.ion)>0) a <- ionize(a)[-is.only.ion]
      species(is.only.ion,delete=TRUE,quiet=TRUE)
      #if(!can.be.numeric(H.act)) thermo$basis$logact[rownames(thermo$basis)=='H+'] <<- H.act
    }

    # more buffer stuff
    if(buffer) {
      args$what <- "logact.basis"
      args$sout <- sout
      logact.basis.new <- logact.basis <- do.call("energy",args)$a
      ibasis.new <- numeric()
      for(k in 1:length(buffers)) {
        ibasis <- which(as.character(mybasis$logact)==buffers[k])
        # calculate the logKs from the affinities
        logK <- a
        for(i in 1:length(logK)) {
          logK[[i]] <- logK[[i]] + thermo$species$logact[i]
          for(j in 1:length(logact.basis.new)) {
            logK[[i]] <- logK[[i]] - logact.basis.new[[j]] * thermo$species[i,j]
            # add ionization correction to proteins
            #if(i %in% is.buffer & length(grep('_',as.character(thermo$species$name[i])))>0 & 
            #  thermo$opt$ionize & rownames(mybasis)[j]=='H+') {
            #  logK[[i]] <- logK[[i]] - logact.basis[[j]] * 
            #    as.data.frame(charge[[match(thermo$species$ispecies[i],names(charge))]]) 
            #}
          }
        }
        lbn <- buffer(logK=logK,ibasis=ibasis,logact.basis=logact.basis.new,
          is.buffer=as.numeric(is.buffer[[which(names(is.buffer)==buffers[k])]]),balance=balance)
        for(j in 1:length(logact.basis.new)) if(j %in% ibasis) logact.basis.new[[j]] <- lbn[[2]][[j]]
        # calculation of the buffered activities' effect on chemical affinities
        is.only.buffer.new <- is.only.buffer[is.only.buffer %in% is.buffer[[k]]]
        for(i in 1:length(a)) {
          if(i %in% is.only.buffer.new) next
          for(j in 1:nrow(mybasis)) {
            # let's only do this for the basis species specified by the user
            # even if others could be buffered
            if(!j %in% ibufbasis) next
            if(!j %in% ibasis) next
            aa <- a[[i]]
            a[[i]] <- aa + (logact.basis.new[[j]] - logact.basis[[j]]) * thermo$species[i,j]
            #if(!identical(a[[i]],aa)) print(paste(i,j))
          }
        }
        if(k==length(buffers) & return.buffer) {
          logact.basis.new <- lbn[[2]]
          ibasis.new <- c(ibasis.new,lbn[[1]])
        } else ibasis.new <- c(ibasis.new,ibasis)
      }
      species(is.only.buffer,delete=TRUE,quiet=TRUE)
      if(length(is.only.buffer) > 0) a <- a[-is.only.buffer]
      # to return the activities of buffered basis species
      tb <- logact.basis.new[unique(ibasis.new)]
      if(!is.null(ncol(tb[[1]]))) {
        nd <- length(which(dim(tb[[1]]) > 1))
        # TODO: apply names for more than two dimensions
        if(nd < 3) {
          for(i in 1:length(tb)) {
            #tb[[i]] <- as.data.frame(tb[[i]])
            if(nd > 0) colnames(tb[[i]]) <- 
              seq(args$lims[[1]][1],args$lims[[1]][2],length.out=args$lims[[1]][3])
            if(nd > 1) rownames(tb[[i]]) <- 
              seq(args$lims[[2]][1],args$lims[[2]][2],length.out=args$lims[[2]][3])
          }
        }
      }
      if(return.buffer) return(tb)
    }

    # more iprotein stuff
    if(!is.null(iprotein)) {
      # 20090331 fast protein calculations
      # function to calculate affinity of formation reactions
      # from those of residues
      loga.protein <- rep(loga.protein,length.out=length(iprotein))
      protein.fun <- function(ip) {
        if(ip %% 50 == 0) cat(paste(ip,"..",sep=""))
        psum(pprod(a[ires],as.numeric(thermo$protein[iprotein[ip],5:25])))-loga.protein[ip]
      }
      # use another level of indexing to let the function
      # report on its progress
      jprotein <- 1:length(iprotein)
      protein.affinity <- mylapply(jprotein,protein.fun)
      if(length(iprotein) > 49) cat("\n")
      ## update the species list
      # we use negative values for ispecies to denote that
      # they index thermo$protein and not thermo$species
      ispecies <- -iprotein
      # the current species list, containing the residues
      resspecies <- thermo$species
      # now we can delete the residues from the species list
      species(ires,delete=TRUE)
      # state and protein names
      state <- resspecies$state[1]
      name <- paste(thermo$protein$protein[iprotein],thermo$protein$organism[iprotein],sep="_")
      # the numbers of basis species in formation reactions of the proteins
      protbasis <- t(t((resspecies[ires,1:nrow(mybasis)])) %*% t((thermo$protein[iprotein,5:25])))
      # put them together
      protspecies <- cbind(protbasis,data.frame(ispecies=ispecies,logact=loga.protein,state=state,name=name))
      myspecies <- rbind(myspecies,protspecies)
      rownames(myspecies) <- 1:nrow(myspecies)
      ## update the affinity values
      names(protein.affinity) <- ispecies
      a <- c(a,protein.affinity)
      a <- a[-ires]
    }

  }

  # put together return values
  T <- args$T
  P <- args$P
  xname <- yname <- ""
  xlim <- ylim <- ""
  xvals <- yvals <- ""
  if(length(args$vars) > 0) {
    xname <- names(args$lims)[1]
    xlim <- args$lims[[1]]
    xvals <- args$vals[[1]]
    if(xname=="T") xvals <- outvert(xvals,"K")
  }
  if(length(args$vars) > 1 & !args$transect) {
    #yname <- args$vars[2]
    yname <- names(args$lims)[2]
    ylim <- args$lims[[2]]
    yvals <- args$vals[[2]]
    if(yname=="T") yvals <- outvert(yvals,"K")
  }

  a <- list(sout=sout,property=property,basis=mybasis,species=myspecies,T=T,P=P,
    xname=xname,xlim=xlim,xvals=xvals,yname=yname,ylim=ylim,yvals=yvals,values=a)
  if(buffer) a <- c(a,list(buffer=tb))
  return(a)

}

ionize <- function(affinity=NULL,other=NULL) {
  # load ionizable groups (if affinity=NULL)
  # or calculate net charges (if return.charge=TRUE) in
  # formation affinities of proteins
  # or unload ionizable groups (if affinity=FALSE)
  neutral <- c('[Asp]','[Glu]','[His]','[Lys]','[Arg]','[Cys]','[Tyr]','[AABB]')
  charged <- c('[Asp-]','[Glu-]','[His+]','[Lys+]','[Arg+]','[Cys-]','[Tyr-]','[AABB+]','[AABB-]')
  if(identical(affinity,FALSE)) {
    s <- species(quiet=TRUE)
    is <- match(c(neutral,charged),s$name)
    if(any(is.na(is))) is <- is[!is.na(is)]
    if(length(is) > 0) species(is,delete=TRUE)
    return(invisible())
  }
  if(is.null(affinity)) {
    is <- species(c(neutral,charged),quiet=TRUE)
    return(is)
  }
  ineutral <- match(neutral,as.character(thermo$species$name))
  ineutral <- c(ineutral,ineutral[length(ineutral)])
  icharged <- match(charged,as.character(thermo$species$name))
  # this is about reacting neutral to charged
  z <- c(1,1,1,1,1,1,1,1,1)
  alpha <- list()
  for(i in 1:length(ineutral)) {
    pK.minus.pH <- affinity[[ineutral[i]]] - affinity[[icharged[i]]]
    alpha[[i]] <- 1 / (1 + 10 ^ (z[i] * (pK.minus.pH)) )
  }
  iprotein <- grep('_',as.character(thermo$species$name))
  # 20090331 added chains here
  iaa <- match(c('Asp','Glu','His','Lys','Arg','Cys','Tyr','chains'),colnames(thermo$protein))
  return.charge <- FALSE
  if(!is.null(other)) {
    if(identical(affinity,other)) return.charge <- TRUE
  }
  for(i in 1:length(iprotein)) {
    po <- s2c(as.character(thermo$species$name[iprotein[i]]),sep='_',keep.sep=FALSE)
    iiprotein <- which(as.character(thermo$protein$protein)==po[1] & 
      as.character(thermo$protein$organism)==po[2])
    b <- function(j) return(thermo$protein[iiprotein,iaa[j]][1])
    if(return.charge) {
      a <- 0
      z <- function(j) return(c(-1,-1,1,1,1,-1,-1,1,-1)[j])
    } else {
      if(is.null(other)) {
        a <- affinity[[iprotein[i]]]
        z <- function(j) return(affinity[[icharged[j]]] - affinity[[ineutral[j]]])
      } else {
        a <- other[[iprotein[i]]]
        z <- function(j) return(other[[icharged[j]]] - other[[ineutral[j]]])
      }
    }
    # sidechain group ionizations
    for(j in 1:(length(iaa)-1)) a <- a + alpha[[j]] * b(j) * z(j) 
    # 20090331 we count numbers of chains now (b(j+1))
    # to determine backbone group ionizations 
    # (these are its cation and anion formation)
    a <- a + alpha[[j+1]] * b(j+1) * z(j+1)
    a <- a + alpha[[j+2]] * b(j+1) * z(j+2)
    affinity[[iprotein[i]]] <- a
  }
  #affinity <- affinity[-c(ineutral,icharged)]
  return(affinity)
}


