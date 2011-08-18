# CHNOSZ/info.R
# search information and thermodynamic properties of species
# 20061024 extraced from species.R jmd

info <- function(species=NULL,states=NULL,quiet=FALSE,return.approx=TRUE) {
  # quiet=TRUE disables checks of self-consistency of ghs
  # and eos/Cp/V values and also makes things run faster

  if(missing(species)) {
    # a friendly summary of thermodynamic information? 20101129
    cat("info: species is NULL; summarizing information about thermodynamic data...\n")
    cat(paste("thermo$obigt has",nrow(thermo$obigt[thermo$obigt$state=="aq",]),"aqueous,",
      nrow(thermo$obigt),"total species\n"))
    cat(paste("number of literature sources: ",nrow(thermo$refs),", elements: ",
      nrow(thermo$element),", buffers: ",length(unique(thermo$buffers$name)),"\n",sep=""))
    cat(paste("number of proteins in thermo$protein is",nrow(thermo$protein),"from",
      length(unique(thermo$protein$organism)),"organisms\n"))
    # print information about SGD.csv, ECO.csv, HUM.csv
    get.protein(organism="SGD")
    get.protein(organism="ECO")
    #get.protein(organism="HUM")
    # print information about yeastgfp.csv
    yeastgfp()
  }

  # argument handling
  species.na <- states.na <- FALSE
  if(!is.null(species)) if(is.na(species[1])) species.na <- TRUE 
  if(!is.null(states)) if(is.na(states[1])) states.na <- TRUE 
  if(species.na | states.na) stop('info: species and/or states arguments are NA.')
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
          iighs <- nrow(tghs)
          if(iighs==0) iighs <- 1
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
            iighs <- 1:nrow(tghs)
          }
          # allow other states if the argument is missing
          else {
            tghs <- thermo$obigt[(thermo$obigt$name %in% species[i] | 
              thermo$obigt$abbrv %in% species[i] | thermo$obigt$formula %in% species[i]),]
            # 20090416 if the name matched put that first
            iname <- match(species[i],tghs$name)[1]
            if(!is.na(iname)) tghs <- tghs[c(iname,(1:nrow(tghs))[-iname]),]
            # except we like O2(g)
            if(species[i]=="O2") {
              iname <- match("oxygen",tghs$name)
              if(!is.na(iname)) tghs <- tghs[c(iname,(1:nrow(tghs))[-iname]),]
            }
            # here we go to the last (duplicated) entry that thermo$opt$level allows
            if(nrow(tghs)==0) iiighs <- 99
            else iiighs <- which(tghs$state==tghs$state[1])
            #iighs <- max(min(thermo$opt$level,max(iiighs)),1)
            iighs <- max(min(1,max(iiighs)),1)
            states[i] <- as.character(tghs$state[iighs])
          }
        }

        # notify the user of multiple matches
        if(nrow(tghs)>1) { 
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
        }

        # notify the user of no matches
        if(nrow(tghs)==0) {
          findnames <- function(species,state=NA,max) {
            if(!is.na(state)) altghs <- thermo$obigt[thermo$obigt$state==state,]
            else altghs <- thermo$obigt
            anames <- unique(c(agrep(species,as.character(altghs$name),value=TRUE,max.distance=max),
                               agrep(species,as.character(altghs$abbrv),value=TRUE,max.distance=max),
                               agrep(species,as.character(altghs$formula),value=TRUE,max.distance=max)))
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
            #  if we want to return the indices of approximate matches
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
          thissource <- thisghs$ref1[j]
          ts2 <- thisghs$ref2[j]
          if(!is.na(ts2)) thissource <- paste(thissource,', ',ts2,sep='')
          td <- thisghs$date[j]
          if(!is.na(td)) thissource <- paste(thissource,', ',td,sep='')
          ii <- ighs[j]
          if(is.na(ii)) ii <- rownames(thisghs)[j]
          cat('info: ',ii,' refers to ',thisghs$name[j],tf,
            ' ',thisghs$state[j],' (',thissource,')\n',sep='')
        }
      }
    }
    # try to return a vector not a list
    t <- try(as.numeric(ighs.list),silent=TRUE)
    if(class(t)=='try-error') return(invisible(ighs.list))
    else return(invisible(as.numeric(ighs.list)))
  }

  # if numeric arguments are given, return ghs and eos parameters
  if(is.numeric(species[1])) {
    nnspecies <- species[species > 0]
    myghs <- thermo$obigt[nnspecies,]
    # species indices don't exceed this value
    ispeciesmax <- nrow(thermo$obigt)
    for(i in 1:length(species)) {
      if(species[i] > ispeciesmax | species[i] < 1) {
        cat(paste("info: there aren't",species[i],"species.\n"))
        next
      }
      # remove scaling factors on EOS parameters depending on state
      # and check them for NAs and consistency with Cp, V values
      if(myghs$state[i]=='aq') {
        # use new obigt2eos function here
        myghs.conv <- obigt2eos(myghs[i,],"aq")
        myghs[i,] <- myghs.conv
        colnames(myghs)[13:20] <- colnames(myghs.conv)[13:20]
        if(!quiet) {
          # check heat capacities
          calcCp <- checkEOS(myghs[i,],"aq","Cp")
          if(!is.na(calcCp) & is.na(myghs$Cp[i])) {
            cat(paste('info: Cp of',myghs$name[i],myghs$state[i],
              'is NA; set by EOS parameters to',round(calcCp,2),'\n'))
            myghs$Cp[i] <- as.numeric(calcCp)
          }
          # check volumes
          calcV <- checkEOS(myghs[i,],"aq","V")
          if(!is.na(calcV) & is.na(myghs$V[i])) {
            cat(paste('info: V of',myghs$name[i],myghs$state[i],
              'is NA; set by EOS parameters to',round(calcV,2),'\n'))
            myghs$V[i] <- as.numeric(calcV)
          }
        }
      } else {
        # for states other than aq
        myghs.conv <- obigt2eos(myghs[i,],"cr")
        myghs[i,] <- myghs.conv
        colnames(myghs)[13:20] <- colnames(myghs.conv)[13:20]
        # set some EOS parameters to zero if they are NA
        # why do we need this? 2010808
#        naEOS <- which(is.na(myghs[i,12:19]))
#        if(length(naEOS)>0) {
#          if(!quiet) cat(paste('info:',c2s(colnames(myghs)[12:19][naEOS],sep=','),'of',
#              myghs$name[i],myghs$state[i],'are NA; set to 0.\n'))
#          myghs[i,naEOS+11] <- 0
#        }
        if(!quiet) {
          calcCp <- checkEOS(myghs[i,],"notaq","Cp")
          if(!is.na(calcCp) & is.na(myghs$Cp[i])) {
            cat(paste('info: Cp of',myghs$name[i],myghs$state[i],
              'is NA; set by EOS parameters to',round(calcCp,2),'\n'))
            myghs$Cp[i] <- as.numeric(calcCp)
          }
        }
      }

      # check the GHS values
      naGHS <- which(is.na(myghs[i,8:10]))
      if(length(naGHS)==1) {
        # calculate a single missing one of G, H, or S
        # we do this even if quiet=TRUE,
        # because NA for one of G, H or S will hamper calculations at high T
        GHS <- GHS(as.character(myghs$formula[i]),DG=myghs[i,8],DH=myghs[i,9],S=myghs[i,10])
        if(!quiet) cat(paste('info: ',c2s(colnames(myghs)[8:10][naGHS]),' of ',
          myghs$name[i],' ',myghs$state[i],' is NA; set to ',round(GHS,2),'.\n',sep=''))
        myghs[i,naGHS+7] <- GHS
      } else if(length(naGHS)==0 & !(quiet)) {
        # check the values
        calcG <- checkGHS(myghs[i,])
        if(is.null(calcG)) {
          cat(paste('info: calculated G of',myghs$name[i],myghs$state[i],
            'is NA!!! (maybe a missing element?)\n'))
        }
      }

    } # end loop over i species

    if(length(unique(myghs$state))!=1 & 'aq' %in% myghs$state) {
      if(!quiet) cat('info: the species are in aqueous and other states\n')
      colnames(myghs)[13:20] <- colnames(thermo$obigt)[13:20]
    }
    return(myghs)
  }
}

