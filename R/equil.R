# CHNOSZ/equil.R
# functions to calculation logarithm of activity
# of species in (metastable) equilibrium

equil.boltz <- function(Astar,nbalance,thisloga) {
  # 20090217 new "abundance" function
  # return logactivity of species
  # works using Maxwell-Boltzmann distribution
  # A/At = e^(Astar/nbalance) / sum(e^(Astar/nbalance))
  # A/At = e^(Astar/nbalance) / sum(e^(Astar/nbalance))
  # where A is activity of the ith residue and
  # At is total activity of residues
  # advantages over abundance.old
  # 1) works on vectors (also matrices - 2D plots now possible)
  # 2) loops over species only - way faster
  # 3) always works (no root finding games)
  # disadvantage:
  # 1) only works for residue reactions
  # 2) can create NaN logacts if the Astars are huge/small

  # initialize output object
  A <- Astar
  # remember the dimensions of elements of Astar (could be vector,matrix)
  Astardim <- dim(Astar[[1]])
  Anames <- names(Astar)
  # first loop: make vectors
  A <- mylapply(1:length(A),function(i) as.vector(A[[i]]))
  # second loop: get the exponentiated Astars (numerators)
  # need to convert /2.303RT to /RT
  #A[[i]] <- exp(log(10)*Astar[[i]]/nbalance[i])/nbalance[i]
  A <- mylapply(1:length(A),function(i) exp(log(10)*Astar[[i]]/nbalance[i]))
  # third loop: accumulate the denominator
  # initialize variable to hold the sum
  At <- A[[1]]; At[] <- 0
  for(i in 1:length(A)) At <- At + A[[i]]*nbalance[i]
  # fourth loop: calculate log abundances and replace the dimensions
  A <- mylapply(1:length(A),function(i) thisloga + log10(A[[i]]/At))
  # fifth loop: replace dimensions
  for(i in 1:length(A)) dim(A[[i]]) <- Astardim
  # add names and we're done!
  names(A) <- Anames
  return(A)
}

equil.react <- function(Astar,av,nbalance,thisloga) {
  # 20090217 extracted from diagram and renamed to abundance.old
  # to turn the affinities/RT (A) of formation reactions into 
  # logactivities of species (logact(things)) at metastable equilibrium
  # 20080217 idea: for any reaction stuff = thing,
  # logQ = logact(thing) - logact(stuff),
  # A = logK - logQ = logK - logact(thing) + logact(stuff),
  # logact(thing) = Astar - A
  # where Astar = logK + logact(stuff)
  # ( Astar = A + logact(thing) )
  # and Abar = ( Astar - logact(thing) ) / n
  # ( or logact(thing) = Astar + Abar * n )
  # where thing has n of the balanced quantity
  # below, j indexes species and i indexes conditions
  # remember the dimensions (could be vector,matrix)
  Adim <- dim(Astar[[1]])
  # first loop: make vectors
  for(i in 1:length(Astar)) {
    Astar[[i]] <- as.vector(Astar[[i]])
    av[[i]] <- as.vector(av[[i]])
  }
  A <- Astar2 <- av

  # some function definitions
  # calculate logact(thing) from Abar and Astar
  activityfun <- function(Abar,j,i) (Astar[[j]][i] - Abar * nbalance[j])
  # difference between total activity of balanced quantity
  # computed from affinities and the mass-balanced quantity (thisloga)
  activityfun2 <- function(Abar,i) {
    act <- 0
    for(j in 1:length(A)) act <- act + (10^activityfun(Abar,j,i))*nbalance[j]
    if(act < 0) {
      act <- -act
      logact <- log10(act)
      diff <- thisloga - logact
    } else {
      logact <- log10(act)
      diff <- logact - thisloga
    }
    return(diff)
  }
  for(i in 1:length(A[[1]])) {
    # gather the min/max values of original affinities
    # to constrain our search interval
    Abar.max <- Abar.min <- NULL
    for(j in 1:length(A)) {
      thisAbar <- (A[[j]][i])/nbalance[j]
      thisAbarstar <- (Astar[[j]][i])/nbalance[j]
      thisAbarstar2 <- (Astar2[[j]][i])/nbalance[j]
      if(!is.infinite(activityfun2(thisAbar,i))) {
        if(is.null(Abar.max)) Abar.max <- thisAbar
          else if(thisAbar > Abar.max) Abar.max <-thisAbar
        if(is.null(Abar.min)) Abar.min <- thisAbar
          else if(thisAbar < Abar.min) Abar.min <-thisAbar
      }
      if(!is.infinite(activityfun2(thisAbarstar,i))) {
        if(is.null(Abar.max)) Abar.max <- thisAbarstar
          else if(thisAbarstar > Abar.max) Abar.max <-thisAbarstar
        if(is.null(Abar.min)) Abar.min <- thisAbarstar
          else if(thisAbarstar < Abar.min) Abar.min <-thisAbarstar
      }
    }
    # make sure A.min < A.max
    #if(Abar.min >= Abar.max) Abar.min <- Abar.max - 1
    # avoid the complication when using uniroot,
    # "f() values at end points not of opposite sign"
    fmin <- function(Abar.min) activityfun2(Abar.min,i)
    fmax <- function(Abar.max) activityfun2(Abar.max,i)
    maxiter <- 1000
    myiter <- 0
    while(sign(fmin(Abar.min))==sign(fmax(Abar.max))) {
      # spread out the values until we have opposite signs 
      diff <- Abar.max - Abar.min
      Abar.min <- Abar.min - diff/2
      Abar.max <- Abar.max + diff/2
      myiter <- myiter + 1
      if(myiter==maxiter) stop(paste('diagram: i tried it',maxiter,
        'times but can\'t make it work :< ... giving up on Abar'))
    }
    # how badly we want the right answer, might
    # have to be adjusted in certain cases
    Atol <- 1e-5
    # find the affinity that gives us the right amount of stuff
    Abar.new <- uniroot(activityfun2,interval=c(Abar.min,Abar.max),i=i,tol=Atol)$root
    # test: did we converge? this shouldn't be required,
    # as uniroot would spit out warnings or errors ... but it doesn't, 
    # even when the tolerance isn't reached by a factor of 100?!
    shouldbezero <- activityfun2(A=Abar.new,i=i)
    if(abs(shouldbezero) > Atol*100) 
      cat(paste('diagram: poor convergence in step ',i,' (remainder in logact of ',
        shouldbezero,')\n',sep=''))
    # and save the activities of the species
    for(j in 1:length(A)) A[[j]][i] <- activityfun(Abar.new,j,i)
  }
  # replace dimensions
  for(i in 1:length(A)) {
    dim(A[[i]]) <- Adim
    dim(av[[i]]) <- Adim
  }
  return(A)
}

