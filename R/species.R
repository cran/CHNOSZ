# CHNOSZ/species.R
# Copyright (C) 2006-2008 Jeffrey M. Dick
# define basis species and species of interest 
# 20051106 renamed from components.R (2005-02-25) jmd

basis <- function(basis=NULL,values=NULL,values2=NULL,delete=FALSE) {

  # define, modify or delete the basis species of a thermodynamic system
  if(delete) {
    thermo$basis <<- NULL
    if(is.null(basis)) return(basis)
  }

  # depending on type, values2 become the logacts of
  # the basis species, or their states
  # called like basis(c('H2O','CO2','H2'),rep('aq',3),c(0,-3,-3))
  # or basis(c('H2O','CO2','H2'),c(0,-3,-3),rep('aq',3))
  # results should be same
  if(!is.null(values2[1])) {
    # first do it with the states (character) so
    # we get the right species
    if(all(is.character(values))) { myv1 <- values; myv2 <- values2 }
    else { myv1 <- values; myv2 <- values }
    basis(basis,myv1)
    return(basis(basis,myv2))
  }

  if(is.null(basis)) {
    if(is.null(thermo$basis)) return(NULL)
    # optionally print the basis definition
    if(sys.nframe()<2) print(thermo$basis)
    # return the basis matrix
    return(invisible(thermo$basis[,1:nrow(thermo$basis),drop=FALSE]))
  }

  # a filter to replace 'aqueous' H2O with 'liq' (avoids warning output)
  if(is.character(values[1])) {
    bw <- match('H2O',basis)
    vw <- match('aq',values)
    t <- bw==vw
    if(!is.na(t)) if(t) values[bw] <- 'liq'
  }

  # keywords, a function and a catch for preset basis definitions
  use.basis <- function(basis,values=NULL) {
    basis.use <- c('CHNOS','CHNOS+','CHNOSe','CHNOPS+','MgCHNOPS+','FeCHNOS','FeCHNOS+')
    if(is.null(basis)) return(basis.use)
    #if(!is.null(thermo$species)) species(delete=TRUE)
    if(!is.null(thermo$basis)) thermo$basis <<- NULL
    basis.logact <- function(b) {
      if(length(b) > 1) return(as.numeric(sapply(b,basis.logact)))
      bases <- c('H2O','CO2','NH3','H2S','O2','H+','e-','Fe2O3')
      logact <- c(0,-3,-4,-7,-80,-7,-7,0)
      ibase <- match(b,bases)
      if(!is.na(ibase)) return(logact[ibase]) else return(-3)
    }
    ibase <- match(basis,basis.use)
    if(ibase==1) b <- c('CO2','H2O','NH3','H2S','O2')
    if(ibase==2) b <- c('CO2','H2O','NH3','H2S','O2','H+')
    if(ibase==3) b <- c('CO2','H2O','NH3','H2S','e-','H+')
    if(ibase==4) b <- c('CO2','H2O','NH3','H3PO4','H2S','e-','H+')
    if(ibase==5) b <- c('Mg+2','CO2','H2O','NH3','H3PO4','H2S','e-','H+')
    if(ibase==6) b <- c('Fe2O3','CO2','H2O','NH3','H2S','O2')
    if(ibase==7) b <- c('Fe2O3','CO2','H2O','NH3','H2S','O2','H+')
    basis(b,basis.logact(b))  
    return(thermo$basis)
  }
  if(basis[1] %in% use.basis(NULL)) return(use.basis(basis[1]))

  # does the supplied dataframe qualify as a basis matrix?
  if(is.data.frame(basis)) {
    # the first test: matrix is square
    if( nrow(basis) != ncol(basis) ) {
    #  warning('basis: non-square matrix.',call.=FALSE)
      return(FALSE)
    }
    # the second test: matrix is invertible
    # to avoid misleading errors, only do this test if we seem to be dealing with basis
    # so far (i.e., the matrix is square and has no zero rows)
    if(class(try(solve(basis),silent=TRUE))=='try-error') { 
    #  warning('basis: non-invertible matrix.',call.=FALSE)
      return(FALSE)
    }
    # they're basis unless proven otherwise above
    return(TRUE)
  }

  # update logact values in thermo$basis
  basis.values <- function(basis,values) {
    values <- rep(values,length.out=length(basis))
    # assign new values for numeric
    # 20090209 these don't seem to be used anymore
    #pHchanged <- FALSE; pechanged <- FALSE
    for(j in 1:length(values)) {
      jold <- which(basis[j]==rownames(thermo$basis) | basis[j]==thermo$obigt$name[thermo$basis$ispecies])
      if(can.be.numeric(values[j])) {
        # fiddling with the name of the electron
        basis[basis=='Z-1'] <- 'e-'
        # check for rownames after the pH,pe transformations
        if(length(jold)==0) {
          warning(paste('basis:',basis[j],
            'is not in the basis. skipping new value of',values[j],'.'),call.=FALSE)
          next()
        }
        thermo$basis$logact[jold] <<- as.numeric(values[j])
      } else {
        # for character, we're looking at states or buffers
          i <- match(basis[j], rownames(thermo$basis))
          if(values[j] %in% as.character(unique(thermo$obigt$state))) {
            # grr... escape this action if different phases are found
            if(values[j]=='cr' & thermo$basis$state[i] %in% c('cr1','cr2','cr3',
              'cr4','cr5','cr6','cr7','cr8','cr9')) values[j] <- as.character(thermo$basis$state[i])
            if(values[j]==thermo$basis$state[i]) is <- thermo$basis$ispecies[i]
            else is <- info(rownames(thermo$basis)[i],values[j],quiet=TRUE)
            if(is.na(is)) {
              cat(paste('basis: state',values[j],'not found for',rownames(thermo$basis)[i],'\n'))
              next
            }
            thermo$basis$ispecies[i] <<- is
            stemp <- as.character(thermo$basis$state); stemp[i] <- values[j]
            thermo$basis$state <<- stemp
            #cat(paste('basis: changed state of',rownames(thermo$basis)[i],'to',values[j],'\n'))
          } else {
            # test if this is the name of a buffer
            if(values[j] %in% as.character(thermo$buffer$name)) {
              ibuff <- which(as.character(thermo$buffer$name)==values[j])
              for(k in 1:length(ibuff)) {
                ispecies <- info(as.character(thermo$buffer$species)[ibuff[k]],as.character(thermo$buffer$state)[ibuff[k]],quiet=TRUE)
                bufmakeup <- makeup(ispecies)
                inbasis <- rownames(bufmakeup) %in% colnames(basis()) 
                if(FALSE %in% inbasis) {
                  cat(paste('basis: the elements',c2s(rownames(bufmakeup)[!inbasis]),'of',thermo$buffer$species[ibuff[k]],'in buffer',values[j],'are not in the basis. skipping.\n'))
                  next
                }
                stemp <- as.character(thermo$basis$logact); stemp[i] <- values[j]
                thermo$basis$logact <<- stemp
              }
            } else {
              cat(paste('basis:',values[j],'is not the name of a state or a buffer.\n'))
            }
          }
        }
    }
    #if(can.be.numeric(values[1])) cat(paste('basis: set logarithm of activity of ',c2s(basis),'.\n',sep=''))
    return(thermo$basis)
  } 

  # main: proceed if basis argument is present
  #if(!is.null(basis) & !(NA %in% basis)) {

  put.basis <- function(basis,states) {
    # indices of the species in obigt
    if(all(is.character(states))) ispecies <- info(basis,states,quiet=TRUE)
    else ispecies <- info(basis,quiet=TRUE)
    if(NA %in% ispecies | is.list(ispecies)) {
      if(is.list(ispecies)) for(i in 1:length(ispecies))
        if(length(ispecies[[i]]) > 1) ispecies[[i]] <- NA
      stop(paste('basis: species',c2s(basis[which(is.na(ispecies))]),'is not available.'),call.=FALSE)
    }
    state <- as.character(thermo$obigt$state)[ispecies]
    # figure out what elements are needed
    elements <- rownames(makeup(ispecies))
    comp <- expand.formula(elements)
    # each row is a basis species
    for(i in 1:length(basis)) {
      # the species makeup
      comp <- rbind(comp,expand.formula(elements,makeup(as.character(thermo$obigt$formula[ispecies[i]]))))
    }
    rownames(comp) <- as.character(thermo$obigt$formula[ispecies])
    colnames(comp) <- elements
    comp <- as.data.frame(comp)
    # order elements alphabetically ('CHNOSZ')
    comp <- comp[,sort(colnames(comp),index.return=TRUE)$ix,drop=FALSE]
    # the formula of the electron is special XXX
    if('Z-1' %in% rownames(comp)) rownames(comp)[rownames(comp)=='Z-1'] <- 'e-'
    # now check it for validity of basis species
    if(!basis(comp)) {
      # values weren't given, and this isn't a basis matrix. stop.
      warning(paste('basis:',nrow(comp),
        'compounds (',c2s(thermo$obigt$formula[ispecies]),')'),call.=FALSE)
      warning(paste('basis:',ncol(comp),
        'elements (',c2s(colnames(comp)),')'),call.=FALSE)
      stop('this is not a valid stoichiometric matrix.')
    }
    comp <- cbind(comp,ispecies,logact,state)
    # assign to the global thermo object
    thermo$basis <<- comp
    return(thermo$basis)
  }

  # get species stoichiometries if numeric first argument
  if(is.numeric(basis[1])) return(basis.comp(basis))

  # character argument: create/modify the definition of basis species 
  if(length(basis)==0) stop("basis: argument was empty.")
  # pH and pe transformations
  if('pH' %in% basis) {
    if(!is.null(values)) values[basis=='pH'] <- -values[basis=='pH']
    basis[basis=='pH'] <- 'H+'
  }
  if(TRUE %in% (c('pe','Eh') %in% basis)) {
    if(!is.null(values)) {
      values[basis=='pe'] <- -values[basis=='pe']
      # 20090209 should be careful with this conversion as it only works at 25 deg C
      # to be sure, just don't call basis('Eh') to be careful
      values[basis=='Eh'] <- -convert(values[basis=='Eh'],'pe')
    }
    basis[basis %in% c('pe','Eh')] <- 'e-'
  }
  if(!length(unique(basis))==length(basis)) stop('basis: please give unique species names.')
  state <- character(0); logact <- numeric(0)
  # state, activity and formulas
  sdefault <- thermo$opt$state; adefault <- 0
  formula <- character()
  for(i in 1:length(basis)) {
    state <- c(state,sdefault); logact <- c(logact,adefault)
    if(!basis[i] %in% thermo$obigt$formula) {
      iname <- match(basis[i],thermo$obigt$name)
      if(length(iname)>0) {
        if(!is.null(thermo$basis)) {
          iold <- match(basis[i],thermo$obigt$name[thermo$basis$ispecies])
          if(!is.na(iold)) formula <- c(formula,rownames(thermo$basis)[iold])
          else formula <- c(formula,as.character(thermo$obigt$formula[iname]))
        } else formula <- c(formula,as.character(thermo$obigt$formula[iname]))
      } else formula <- c(formula,basis[i])
    } else formula <- c(formula,basis[i])
  }
  if('Z-1' %in% formula) formula[formula=='Z-1'] <- 'e-'
  # if these basis species are already present, update their values
  # or if no values are present and all the requested species are 
  # in the current basis, just return the current basis
  if(!is.null(thermo$basis)) {
    # replace names of basis species with their formulas
    #myformula <- formula
    #iinobigt <- match(myformula,thermo$obigt$name[thermo$basis$ispecies])
    #myformula[!is.na(iinobigt)] <- rownames(thermo$basis)[iinobigt[!is.na(iinobigt)]]
    #print(myformula)
    if(!FALSE %in% (formula %in% rownames(thermo$basis))) {
      if(!is.null(values)) return(basis.values(formula,values)) 
      #if(length(formula)==nrow(thermo$basis)) return(thermo$basis)
    }
  }
  # remove old species object
  species.remove <- function() {
    ts <- NULL
    if(!is.null(thermo$species)) {
      ts <- thermo$species
      thermo$species <<- NULL
      #cat('basis: removed species definitions.\n')
    }
    return(ts)
  }
  # keep track of the old basis species
  basis.old <- thermo$basis
  mystates <- rep("",length(basis))
  for(i in 1:length(mystates)) {
    doit <- TRUE
    #if(is.null(values)) doit <- TRUE
    #else if(can.be.numeric(values[i])) doit <- TRUE
    if(doit) {
      iold <- which(basis[i]==rownames(basis.old) | basis[i]==thermo$obigt$name[basis.old$ispecies])
      if(!is.null(values)) {
        if(values[i] %in% unique(thermo$obigt$state)) mystates[i] <- values[i]
        else if(length(iold) > 0) mystates[i] <- as.character(basis.old$state)[iold]
      } else if(length(iold) > 0) mystates[i] <- as.character(basis.old$state)[iold]
    }
  }
  # provisionally load new basis species
  basis.new <- put.basis(basis,mystates)
  # if (not state) values are specified use them
  isstate <- FALSE
  if(!is.null(values)) isstate <- all(values %in% thermo$obigt$state)
  if(!is.null(values) & !isstate) {
    species.remove()
    return(basis.values(formula,values))
  } else {
    # otherwise, if old and new basis species have same elements,
    # maintain elemental chemical potentials (at Tr, Pr)
    # but first, check that the requested basis is actually different
    if(identical(basis.old$ispecies,basis.new$ispecies)) return(thermo$basis)
    if(!is.null(basis.old)) {
      cat(paste('basis: changed basis to ',c2s(rownames(basis.new)),'.\n',sep=''))
    }
    if(!is.null(basis.old)) {
      # isolate stoichiometric matrices
      basis.old.mat <- basis.old[,1:nrow(basis.old),drop=FALSE]
      basis.new.mat <- thermo$basis[,1:nrow(thermo$basis),drop=FALSE]
      # order columns (elements)
      basis.old.mat <- basis.old.mat[,sort(colnames(basis.old.mat),index.return=TRUE)$ix,drop=FALSE]
      basis.new.mat <- basis.new.mat[,sort(colnames(basis.new.mat),index.return=TRUE)$ix,drop=FALSE]
      # proceed if the elements are the same
      if(identical(colnames(basis.old.mat),colnames(basis.new.mat))) {
        cat('basis: elements unchanged. maintaining elemental chemical potentials.\n')
        # chemical potentials of the old basis species
        basis.old.mu <- thermo$obigt$G[basis.old$ispecies] - convert(basis.old$logact,'g')
        # chemical potentials of the elements
        element.mu <- solve(basis.old.mat,basis.old.mu)
        # chemical potentials and activities of the new basis species
        basis.new.mu <- numeric()
        for(i in 1:nrow(basis.new.mat)) basis.new.mu <- c(basis.new.mu,
          sum(element.mu*basis.new.mat[i,]) - thermo$obigt$G[basis.new$ispecies[i]])
        basis.new.logact <- - convert(basis.new.mu,'logk')
        # eliminate small numbers
        for(j in 1:length(basis.new.logact)) 
          if(abs(basis.new.logact[j]) < thermo$opt$cutoff) basis.new.logact[j] <- 0
        # update basis and restore species
        ts <- species.remove()
        basis.values(formula,basis.new.logact)
        if(!is.null(ts)) {
          cat('basis: restoring species using new basis definition.\n')
          species(ts$ispecies,quiet=TRUE)
          species(1:nrow(thermo$species),ts$logact,quiet=TRUE)
        }
        return(thermo$basis)
      }
    }
    species.remove()
    return(thermo$basis)
  }
}


species <- function(species=NULL,state=NULL,delete=FALSE,quiet=TRUE) {
# 20080925 changed default to quiet=TRUE 
# (TODO: update CHNOSZ via Rpad accordingly)
  
  missingstate <- missing(state)
  state <- state.args(state)

  # if no species or states are given, just return the species list
  if(is.null(species) & is.null(state) & !delete) return(thermo$species)
  # if no species are given use all of them if available
  if(is.null(species) & !is.null(thermo$species)) species <- 1:nrow(thermo$species)

  if(length(species) > length(state)) state <- rep(state,length.out=length(species)) else 
  if(length(state) > length(species)) species <- rep(species,length.out=length(state))

  # if they don't look like states (aq,gas,cr) or activities (numeric), 
  # use them as a suffix for species name (e.g., a protein-organism)
  if( length(which(state %in% unique(as.character(thermo$obigt$state)))) < 
    length(state) & !can.be.numeric(state[1]) & !can.be.numeric(species[1]) ) {
      for(i in 1:length(state)) species[i] <- paste(species[i],'_',state[i],sep='')
      state <- rep(thermo$opt$state,length.out=length(state))
  }

  # delete a specie entry or warn that it won't be duplicated
  if(delete) {
    if(is.null(species)) {
      if(!is.null(thermo$species)) {
        cat('species: removing species definitions.\n')
        thermo$species <<- NULL
      }
      return(thermo$species)
    }
    if(is.character(species[1])) {
      is <- match(species,thermo$species$name)
      if(NA %in% is) 
        cat(paste('species:',c2s(species[is.na(is)]),'not deleted because it is not present.\n'))
      is <- is[!is.na(is)]
    } else is <- species
    isspecies <- is %in% 1:nrow(thermo$species)
    if(FALSE %in% isspecies) 
      cat(paste('species:',c2s(is[!is%in%nrow(thermo$species)]),'not deleted because it is not present.\n'))
    is <- is[isspecies]
    if(length(is)==0) return(species())
    thermo$species <<- thermo$species[-is,]
    if(nrow(thermo$species)==0) thermo$species <<- NULL
    else rownames(thermo$species) <<- 1:nrow(thermo$species)
  } else {
      # append/change species entries
      if(is.null(thermo$basis)) {
        stop('first you must define a basis.')
      }
      if(is.character(species[1])) {
        # character first argument, species in thermo$obigt
        # but only give states if they are numeric
        is <- NULL; if(!can.be.numeric(state[[1]])) is <- state
        ispecies <- info(species,is,quiet=TRUE)
        ispecies <- ispecies[!is.na(ispecies)]
        if(length(ispecies)==0) return(species())
        was.character <- TRUE
      } else if(is.numeric(species[1])) {
        ispecies <- species
        ispecies <- ispecies[!is.na(ispecies)]
        if(length(ispecies)==0) return(species())
        was.character <- FALSE
      }
      jspecies <- ispecies
      myspecies <- NULL
      if(!is.null(thermo$species)) 
        if(TRUE %in% (ispecies %in% thermo$species$ispecies)) {
          myspecies <- match(ispecies,thermo$species$ispecies)
          myspecies <- thermo$species$ispecies[myspecies[!is.na(myspecies)]]
          ispecies <- ispecies[!ispecies%in%myspecies]
        }
      # only add to an existing dataframe if the indices can't
      # all possibly refer to the rows of the dataframe
      doit <- TRUE
      ## might (not) work well when you want to add e.g. H2O after some others
      if(!is.null(thermo$species)) if(all(jspecies %in% 1:nrow(thermo$species))) 
        if(!was.character) doit <- FALSE
      if(length(ispecies) > 0 & !(is.numeric(ispecies[1]) & is.numeric(state[1]) ) & doit) {
        f <- (basis(ispecies))
        # the default states and activities
        state <- as.character(thermo$obigt$state[ispecies])
        logact <- numeric()
        for(i in 1:length(state)) {
          if(length(agrep('missing',rownames(f)[i]))>0) la <- NA
          else { if(state[i]=='aq') la <- -3 else la <- 0 }
          logact <- c(logact,la)
        }
        # yes, the species names too
        name <- as.character(thermo$obigt$name[ispecies])
        # add the ispecies values
        t <- data.frame(f,ispecies=ispecies,logact=logact,state=state,name=name,stringsAsFactors=FALSE)
        # nasty for R, but "H2PO4-" looks better than "H2PO4."
        colnames(t)[1:nrow(thermo$basis)] <- rownames(thermo$basis)
        if(is.null(thermo$species)) thermo$species <<- t else 
          thermo$species <<- rbind(thermo$species,t)
        rownames(thermo$species) <<- seq(1:nrow(thermo$species))
      }
      #  update activities or states
        if(!is.null(state)) {
          state <- rep(state,length.out=length(jspecies))
          # if is looks like species aren't set yet, try to do so
          if(is.null(thermo$species)) { 
            species(jspecies)
          } else {
            mj <- jspecies[!jspecies %in% thermo$species$ispecies]
            if(!can.be.numeric(species[1])) species(mj)
          }
          # we bet that the number of rows is smaller
          # than the indices of whatever species we have
          if(can.be.numeric(species[1]) & max(jspecies) <= nrow(thermo$species))
            jspecies <- thermo$species$ispecies[jspecies]
          mj <- match(jspecies,thermo$species$ispecies)
          if(can.be.numeric(state[1])) {
            if(NA %in% mj[1]) warning(paste('can\'t update activity of species',
              c2s(which(is.na(mj))),' requested'),call.=FALSE)
            thermo$species$logact[mj] <<- state
          } else {
            mj <- match(jspecies,thermo$species$ispecies)
            state <- rep(state,length.out=length(mj))
            name <- thermo$species$name[mj]
            # try to check that the states actually exist
            for(k in 1:length(mj)) {
              doit <- TRUE
              if(NA %in% mj[k]) doit <- FALSE
              myform <- thermo$obigt$formula[thermo$species$ispecies[mj[k]]]
              # 20080925 don't match formula -- two proteins might have the
              # same formula (e.g. YLR367W and YJL190C)
              #iobigt <- which(thermo$obigt$name==thermo$species$name[mj[k]] | thermo$obigt$formula==myform)
              iobigt <- which(thermo$obigt$name==thermo$species$name[mj[k]])
              if(!state[k] %in% thermo$obigt$state[iobigt]) 
                doit <- FALSE
              if(!doit) warning(paste('can\'t update state of species ',
                mj[k],' to ',state[k],'.\n',sep=''),call.=FALSE)
              else {
                ii <- match(state[k],thermo$obigt$state[iobigt])
                thermo$species$state[mj[k]] <<- state[k]
                thermo$species$name[mj[k]] <<- thermo$obigt$name[iobigt[ii]]
              }
              #if(TRUE %in% can.be.numeric(s2c(thermo$species$state[mj]))) {
              #  warning(paste("forbidden: i won't change",
              #    thermo$species$state[mj],'of',thermo$species$name[mj]))
              #  next
              #}
            }
          }
        } else {
          # this message turns out to be kinda distracting. what should be here?
          #if(!is.null(myspecies)) 
            #cat(paste('species: keeping ',c2s(thermo$obigt$name[myspecies],sep=', '),'.\n',sep='')) 
        }
      if(quiet) { if(sys.nframe() < 2) print(species()) }
      else print(species())
      return(invisible(match(jspecies,thermo$species$ispecies)))
  }
  if(quiet) { if(sys.nframe() < 2) print(species()) }
  else print(species())
  return(invisible(nrow(thermo$species)))
}

