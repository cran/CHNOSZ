# CHNOSZ/info.R
# Copyright (C) 2006-2009 Jeffrey M. Dick
# search information and thermodynamic properties of species
# 20061024 extraced from species.R jmd

info <- function(species=NULL,states=NULL,quiet=FALSE,file='',return.approx=TRUE) {
  # file argument: where to cat the messages about inconsistencies among
  # GHS, and between eos parameters and Cp, V ('' for console)

  # quiet=TRUE disables checks of self-consistency of ghs
  # and eos/Cp/V values and also makes things run faster
  #if(missing(quiet)) if(sys.nframe()>1) quiet <- TRUE

  # if species is missing, run a consistency check for each species,
  # returning a nice table at the end ...
  missing.species <- FALSE
  if(missing(species)) {
    missing.species <- TRUE
    species <- 1:nrow(thermo$obigt)
    cat(paste('info: checking consistency of parameters for',length(species),'species.\n'))
    #for(i in 1:nrow(thermo$obigt)) {
    #  if(i %% 10 == 0) cat(paste(i,'\n'),file=file,append=TRUE)
    #  t <- info(i)
    #}
  }

  # show license or release notes
  if(identical(species,'license') | identical(species,'licence') | identical(species,'copying')) {
    licensefile <- file.path(system.file(package="CHNOSZ"), "COPYING")
    file.show(licensefile)
    return()
  }
  if(identical(species,'CHNOSZ') | identical(species,'news')) {
    newsfile <- file.path(system.file(package="CHNOSZ"), "NEWS")
    file.show(newsfile)
    return()
  }

  species.na <- states.na <- FALSE
  if(!is.null(species)) if(is.na(species[1])) species.na <- TRUE 
  if(!is.null(states)) if(is.na(states[1])) states.na <- TRUE 
  if(species.na | states.na) stop('info: species and/or states arguments are NA.')

  # argument handling
  if(missing(states) | is.null(states)) missing.states <- TRUE else missing.states <- FALSE
  if(!is.null(states) & length(species)!=length(states))
    states <- rep(states,length.out=length(species))  
  # fill-in missing state with default setting
  if(is.null(states)) states <- rep(thermo$opt$state,length(species))


  # if arguments are character, return the matching rownumbers in thermo$obigt
  if(is.character(species[1])) {
    # first make sure the rownames are their numbers
    rownames(thermo$obigt) <<- 1:nrow(thermo$obigt)
    # the lists of species addresses
    ighs <- numeric(0); ieos <- numeric(0)
    ighs.list <- list()
    # dataframe for telling the user what species are found
    thisghs <- thermo$obigt[0,]
    for(i in 1:length(species)) {
        # if there's more than one match use the first by default
        iighs <- 1

        use.other.states <- FALSE
        if(length(grep('_',species[i]))>0) {
          # we're dealing with a protein
          tghs <- thermo$obigt[(thermo$obigt$name %in% species[i]) & thermo$obigt$state %in% states[i],]
          # try to add up protein
          if(nrow(tghs)==0) {
          # split the name at the underscore
            us <- match(TRUE, s2c(species[i]) == '_')
            protein <- substr(species[i],1,us-1)
            organism <- substr(species[i],us+1,nchar(species[i]))
            ip <- protein(protein,organism)
            # did we find a protein? add its properties to obigt
            if(length(ip) > 0) {
              newrow <- protein(ip,states[i])
              colnames(newrow) <- colnames(thermo$obigt)
              thermo$obigt <<- rbind(thermo$obigt,newrow)
              rownames(thermo$obigt) <<- 1:nrow(thermo$obigt)
              tghs <- thermo$obigt[nrow(thermo$obigt),]
            }
          }
        } else {
          # selection criteria: search names and abbreviations and formulas
          # and state, if the states argument is supplied
          if(missing.states) use.other.states <- TRUE
          else if(states[i]=="") use.other.states <- TRUE 
          if(!use.other.states) {
            if(states[i]=='cr') searchstates <- c('cr','cr1','cr2','cr3',
              'cr4','cr5','cr6','cr7','cr8','cr9') else searchstates <- states[i]
            tghs <- thermo$obigt[(thermo$obigt$name %in% species[i] | 
              thermo$obigt$abbrv %in% species[i] | thermo$obigt$formula %in% species[i]) & 
              thermo$obigt$state %in% searchstates,]
            if(states[i]=='cr') iighs <- 1:nrow(tghs)
          }
          # allow other states if the argument is missing
          else {
            tghs <- thermo$obigt[(thermo$obigt$name %in% species[i] | 
              thermo$obigt$abbrv %in% species[i] | thermo$obigt$formula %in% species[i]),]
            # 20090416 if the name matched put that first
            # (e.g. so that "oxygen" comes up as a gas)
            iname <- match(species[i],tghs$name)[1]
            if(!is.na(iname)) tghs <- tghs[c(iname,(1:nrow(tghs))[-iname]),]
            # look for the species with the preferred state
            # changed to NA 20090416 -- only do H2 and O2
            desiredstate <- NA
            # special handling for O2 and H2
            # FIXME: remove this
            if(species[i] %in% c('O2','H2')) desiredstate <- 'gas'
            istate <- match(desiredstate,tghs$state)
            if(!is.na(istate)) iighs <- istate else iighs <- 1
            states[i] <- as.character(tghs$state[iighs])
          }
        }

        # notify the user of multiple matches
        if(nrow(tghs)>1 & states[i]!='cr') { 
          if(!quiet) {
            tg <- tghs[,1:5]
            if(length(unique(tg$formula))==1 & length(unique(tg$name))==1) {
              otext <- ''
              if(species[i] %in% tg$formula) otext <- tg$name[1]
              if(species[i] %in% tg$name) otext <- tg$formula[1]
              if(otext==species[i]) otext <- '' else otext <- paste(' (',otext,')',sep='')
              stext <- ''
              for(j in 1:nrow(tg)) {
                if(j==nrow(tg)) ntext <- '.\n' else ntext <- ', '
                stext <- paste(stext,tg$state[j],ntext,sep='')
              }
              cat(paste('info: ',species[i],otext,' available in ',stext,sep=''))
            } else {
              cat(paste('info:',species[i],'matches these species:\n'))
              print(tghs[,1:5])
            }
              
          } 
          #if(length(unique(tghs$state))>1) {
          #    ostate <- unique(as.character(tghs$state))
          #    ostate <- ostate[-match(states[i],ostate)]
          #    cat(paste('info: matched ',tghs$formula[iighs],' ',states[i],' (',tghs$name[iighs],') ... other states are ',c2s(ostate),'.\n',sep=''))
          #}
        }

        # notify the user of no matches
        if(nrow(tghs)==0) {
          findnames <- function(species,state=NA,max) {
            if(!is.na(state)) altghs <- thermo$obigt[thermo$obigt$state==state,]
            else altghs <- thermo$obigt
            anames <- unique(c(agrep(species,as.character(altghs$name),value=TRUE,max=max),
                               agrep(species,as.character(altghs$abbrv),value=TRUE,max=max),
                               agrep(species,as.character(altghs$formula),value=TRUE,max=max)))
            if(length(anames)>0) {
              # finish off the message and print the approximate matches
              # 20090301 observe quiet
              if(!quiet) cat('.\n')
              # which species are approximately matching
              iga <- which(altghs$name %in% anames | 
                altghs$abbrv %in% anames | altghs$formula %in% anames,1:4)
              # print species names
              # 20090301 only if return.approx is TRUE
              if(return.approx) {
                if(length(anames) > 20) {
                  cat('info: similar species names, abbreviations, or formulas are:\n')
                  nt <- 200
                  if(length(anames) > nt) {
                    cat('info: (truncated at) ',nt,'.\n',sep='')
                    anames <- anames[1:200]
                  }
                  print(anames)
                # print species info
                } else {
                  cat('info: approximately matching species are:\n')
                  print(altghs[iga,1:4])
                }
              }
              return(iga)  
            } else return()
          }
          if(use.other.states) {
            # 20090301 observe quiet
            if(!quiet) cat(paste('info: no match for ',species[i],sep=''))
            t <- findnames(species[i],max=0.1)
            if(is.null(t)) t <- findnames(species[i],max=0.3)
          } else {
            if(!quiet) cat(paste('info: no match for ',species[i],' ',states[i],sep=''))
            t <- findnames(species[i],states[i],max=0.1)
            if(is.null(t)) t <- findnames(species[i],states[i],max=0.3)
          }
          if(is.null(t)) {
            if(!quiet) cat(', and no approximate matches.\n')
            if(return.approx) ighs.list[length(ighs.list)+1] <- NA
          } else {
            #  if we want to return the indices of approximat matches
            if(return.approx) ighs.list[length(ighs.list)+1] <- list(t)
          }
        } 
        tighs <- as.numeric(rownames(tghs)[iighs])
        ighs <- c(ighs,tighs)
        if(!is.na(tighs[1])) for(i in 1:length(tighs))
          ighs.list[[length(ighs.list)+1]] <- tighs[i]
        # build the table of species to show the user
        if(!is.na(tighs[[1]])) thisghs <- rbind(thisghs,tghs[iighs,])
    }
    if(nrow(thisghs)>0) {
      if(quiet) {
        #cat(paste('info: matching species indices are ',c2s(rownames(thermo$obigt)[ighs],sep=', '),'.\n',sep=''))
      } else {
        #if(nrow(thisghs)==1) {iword <- 'index'; word <- 'this'}
        #else {iword <- 'indices'; word <- 'these'}
        #cat('info: returning',iword,'for',word,'species:\n')
        #print(thisghs[,1:7])
        for(j in 1:nrow(thisghs)) {
          tf <- thisghs$formula[j]
          if(tf %in% species) tf <- '' else tf <- paste(', ',tf,sep='')
          thissource <- thisghs$source1[j]
          ts2 <- thisghs$source2[j]
          if(!is.na(ts2)) thissource <- paste(thissource,', ',ts2,sep='')
          td <- thisghs$date[j]
          if(!is.na(td)) thissource <- paste(thissource,', ',td,sep='')
          ii <- ighs[j]
          if(is.na(ii)) ii <- rownames(thisghs)[j]
          cat('info: ',ii,' refers to ',thisghs$name[j],tf,
            ' ',thisghs$state[j],' (',thissource,').\n',sep='')
        }
      }
    }
    # try to return a vector not a list
    t <- try(as.numeric(ighs.list),silent=TRUE)
    if(class(t)=='try-error') return(invisible(ighs.list))
    else return(invisible(as.numeric(ighs.list)))
  }

  # if numeric arguments are given, return ghs and eos
  if(is.numeric(species[1])) {
    nnspecies <- species[species > 0]
    myghs <- thermo$obigt[nnspecies,]
    # to keep track of the results of consistency checks
    DCp <- DV <- DCp.supcrt <- DV.supcrt <- DG <- Di <- numeric()
    Dname <- Dstate <- character()
    #for(i in 1:nrow(myghs)) {
    for(i in 1:length(species)) {
      if(species[i] > nrow(thermo$obigt) | species[i] < 1) {
        cat(paste("info: there aren't",species[i],"species.\n"))
        next
      }
      dCp <- dV <- dCp.supcrt <- dV.supcrt <- dG <- NA
      # convert the eos parameters depending on state
      # and check them for NAs and consistency with Cp, V values
      if(myghs$state[i]=='aq') {
        myghs[i,13:20] <- myghs[i,13:20] * 10^c(-1,2,0,4,0,4,5,0)
        colnames(myghs)[13:20] <- c('a1','a2','a3','a4','c1','c2','omega','Z') 
        # "quiet" also means quick ... to skip some of these checks
        #if(!quiet) {
        naEOS <- which(is.na(myghs[i,13:20]))
        #if(length(naEOS)>0) {
        #  cat(paste('info:',c2s(colnames(myghs)[13:20][naEOS],sep=','),'of',
        #    myghs$name[i],myghs$state[i],'are NA; set to 0.\n'))
        #  myghs[i,naEOS+12] <- 0
        #}
        #}
        if(!quiet) {
          # value of X consistent with CHNOSZ
          X <- -2.773788E-7
          # value of X consistent with SUPCRT
          X <- -3.055586E-7
          Cp <- myghs$c1[i] + myghs$c2[i]/(298.15-228)^2 + myghs$omega[i]*298.15*X
          if(!is.na(myghs$Cp[i])) {
            if(abs(Cp-myghs$Cp[i])>1) {
              cat(paste('info: Cp (from EOS) of',
                myghs$name[i],myghs$state[i],'differs by',round(Cp-myghs$Cp[i],2),
                'from tabulated value.\n'),file=file,append=TRUE)
              dCp.supcrt <- round(Cp-myghs$Cp[i],2)
            }
          } else {
            if(!is.na(Cp)) {
              cat(paste('info: Cp of ',myghs$name[i],' ',myghs$state[i],
                ' is NA; set by EOS parameters to ',round(Cp,2),'.\n',sep=''))
              myghs$Cp[i] <- as.numeric(Cp)
            }
          }
        }
        if(!quiet) {
          # value of Q consistent with IAPWS95/AW90
          Q <- 0.00002483137
          # value of Q consistent with SUPCRT92
          Q <- 0.00002775729
          V <- convert(myghs$a1[i],'cm3bar') + convert(myghs$a2[i],'cm3bar')/(2601) + 
            (convert(myghs$a3[i],'cm3bar') + convert(myghs$a4[i],'cm3bar')/(2601))/(298.15-228) - 
              Q * myghs$omega[i]
          if(!is.na(myghs$V[i]) & !is.na(V)) {
            if(abs(V-myghs$V[i])>1) {
              cat(paste('info: V (from EOS) of',
                myghs$name[i],myghs$state[i],'differs by',round(V-myghs$V[i],2),
                'from tabulated value.\n'),file=file,append=TRUE)
              dV.supcrt <- round(V-myghs$V[i],2)
            }
          } else {
            if(!is.na(V)) {
              cat(paste('info: V of ',myghs$name[i],' ',myghs$state[i],
                ' is NA; set by EOS parameters to ',round(V,2),'.\n',sep=''))
              myghs$V[i] <- as.numeric(V)
            }
          }
        }
      } else {
        myghs[i,13:20] <- myghs[i,13:20] * 10^c(0,-3,5,0,-5,0,0,0)
        colnames(myghs)[13:20] <- c('a','b','c','d','e','f','lambda','T')
        #if(!quiet) {
        naEOS <- which(is.na(myghs[i,12:19]))
        if(length(naEOS)>0) {
            cat(paste('info:',c2s(colnames(myghs)[12:19][naEOS],sep=','),'of',
              myghs$name[i],myghs$state[i],'are NA; set to 0.\n'))
          myghs[i,naEOS+11] <- 0
        }
        #}
        if(!quiet) {
          Cp <- cgl('Cp',eos=myghs[i,],ghs=myghs[i,],T=thermo$opt$Tr,P=thermo$opt$Pr)[[1]]
          if(!is.na(myghs$Cp[i]) & !is.na(Cp)) {
            if(abs(Cp-myghs$Cp[i])>1) {
              cat(paste('info: Cp (from EOS) of',
                myghs$name[i],myghs$state[i],'differs by',round(Cp-myghs$Cp[i],2),
                'from tabulated value.\n'),file=file,append=TRUE)
              dCp <- round(Cp-myghs$Cp[i],2)
            }
          } else {
            if(!is.na(Cp)) {
              cat(paste('info: Cp of ',myghs$name[i],' ',myghs$state[i],
                ' is NA; set by EOS parameters to ',round(Cp,2),'.\n',sep=''))
              myghs$Cp[i] <- as.numeric(Cp)
            }
          }
        }
      }
      # check the GHS values
      #if(!quiet) {
        naGHS <- which(is.na(myghs[i,8:10]))
        j <- as.numeric(rownames(myghs)[i])
        if(length(naGHS)>1) {
          # don't do this. keep NAs in there
          ##cat(paste('info:',c2s(colnames(myghs)[8:10][naGHS],sep=','),'of',
          ##  myghs$name[i],myghs$state[i],'are NA; set to 0.\n'))
          ##myghs[i,naGHS+7] <- 0
        } else if(length(naGHS)==1) {
          GHS <- GHS(as.character(myghs$formula[i]),DG=myghs[i,8],DH=myghs[i,9],S=myghs[i,10])
          if(!quiet) cat(paste('info: ',c2s(colnames(myghs)[8:10][naGHS]),' of ',
            myghs$name[i],' ',myghs$state[i],' is NA; set to ',round(GHS,2),'.\n',sep=''))
          myghs[i,naGHS+7] <- GHS
        } else {
          G <- GHS(as.character(myghs$formula[i]),DG=myghs[i,8],DH=myghs[i,9],S=myghs[i,10])
          warn <- FALSE
          if(!is.na(G)) {
            if(abs(G-myghs[i,8])>1000) warn <- TRUE
            #if(myghs[i,8]!=0) if(abs((G-myghs[i,8])/myghs[i,8])>0.05) warn <- TRUE
            if(warn) {
              if(!quiet) cat(paste('info: G (from H and S) of',myghs$name[i],myghs$state[i],'differs by',
                round(G-myghs[i,8]),'from tabulated value.\n'),file=file,append=TRUE)
              dG <- round(G-myghs[i,8])
            }
          } else {
            if(!quiet) cat(paste('info: G of',myghs$name[i],myghs$state[i],'is NA!!! (maybe a missing element?)\n'),file=file,append=TRUE)
          }
        }
      #}
      if(missing.species) {
        if(!(is.na(dCp)&is.na(dV)&is.na(dCp.supcrt)&is.na(dV.supcrt)&is.na(dG)) ) {
          Dname <- c(Dname,as.character(myghs$name[i]))
          Dstate <- c(Dstate,as.character(myghs$state[i]))
          DCp <- c(DCp,dCp)
          DV <- c(DV,dV)
          DCp.supcrt <- c(DCp.supcrt,dCp.supcrt)
          DV.supcrt <- c(DV.supcrt,dV.supcrt)
          DG <- c(DG,dG)
          Di <- c(Di,species[i])
        }
      }
    } # end loop over i species

    if(length(unique(myghs$state))!=1 & 'aq' %in% myghs$state) {
      cat('info: the species are in aqueous and other states\n')
      colnames(myghs)[13:20] <- colnames(thermo$obigt)[13:20]
    }
    if(missing.species | all(is.na(myghs))) {
      # for some reason, DCp likes to become list. make it numeric
      DCp <- as.numeric(DCp)
      #t <- data.frame(name=Dname,state=Dstate,DCp=DCp,DCp.supcrt=DCp.supcrt,DV=DV,DV.supcrt=DV.supcrt,DG=DG)
      t <- data.frame(ispecies=Di,name=Dname,state=Dstate,DCp=DCp,DV=DV.supcrt,DG=DG)
      return(t)
    } else return(myghs)
  }
}

