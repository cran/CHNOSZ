# CHNOSZ/species.R
# define species of interest 

species <- function(species=NULL,state=NULL,delete=FALSE,quiet=FALSE) {
# 20080925 changed default to quiet=TRUE 
# 20101003 changed default to quiet=FALSE
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
        ispecies <- numeric()
        for(i in 1:length(species)) {
          ispec <- info(species[i],is[i],quiet=TRUE)
          #ispecies <- c(ispecies,ispec[min(thermo$opt$level,length(ispec))])
          ispecies <- c(ispecies,ispec[min(1,length(ispec))])
        }
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
            species(jspecies,quiet=TRUE)
          } else {
            mj <- jspecies[!jspecies %in% thermo$species$ispecies]
            if(!can.be.numeric(species[1])) species(mj,quiet=TRUE)
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
              #iobigt <- which(thermo$obigt$name==thermo$species$name[mj[k]] | thermo$obigt$formula==myform)
              # 20080925 don't match formula -- two proteins might have the
              # same formula (e.g. YLR367W and YJL190C)
              #iobigt <- which(thermo$obigt$name==thermo$species$name[mj[k]])
              # 20091112 do match formula if it's not a protein -- be able to 
              # change "carbon dioxide(g)" to "CO2(aq)"
              if(length(grep("_",thermo$species$name[mj[k]])) > 0)  
                iobigt <- which(thermo$obigt$name==thermo$species$name[mj[k]])
              else {
                iobigt <- which(thermo$obigt$name==thermo$species$name[mj[k]] & thermo$obigt$state==state[k])
                if(length(iobigt)==0)
                  iobigt <- which(thermo$obigt$name==thermo$species$name[mj[k]] | thermo$obigt$formula==myform)
              }
              if(!state[k] %in% thermo$obigt$state[iobigt]) 
                doit <- FALSE
              if(!doit) warning(paste('can\'t update state of species ',
                mj[k],' to ',state[k],'.\n',sep=''),call.=FALSE)
              else {
                ii <- match(state[k],thermo$obigt$state[iobigt])
                thermo$species$state[mj[k]] <<- state[k]
                thermo$species$name[mj[k]] <<- thermo$obigt$name[iobigt[ii]]
                thermo$species$ispecies[mj[k]] <<- as.numeric(rownames(thermo$obigt)[iobigt[ii]])
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
      if(!quiet) print(species(quiet=TRUE))
      return(invisible(match(jspecies,thermo$species$ispecies)))
  }
  #if(quiet) { if(sys.nframe() < 2) print(species()) }
  #else print(species())
  if(!quiet) print(species(quiet=TRUE))
  return(invisible(nrow(thermo$species)))
}

