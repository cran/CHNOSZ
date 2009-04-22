# CHNOSZ/makeup.R
# Copyright (C) 2006-2008 Jeffrey M. Dick
# manipulate objects describing chemical compositions of species
# 20060806 jmd

makeup <- function(compound='',component=NULL) {
  # formula: character, such as C2H3O89
  # makeup: dataframe equiv of formula (1 column, rownames=elements)
  # basis: dataframe equiv of formula (1 row, colnames=basis species)
  
  # if both arguments are dataframe, return their sum
  if(is.data.frame(compound) & is.data.frame(component)) {
    elements <- unique(c(rownames(compound),rownames(component)))
    nelement <- rep(NA,length(elements))
    for(i in 1:length(elements)) {
      n <- 0
      n1 <- match(elements[i],rownames(compound))
      n2 <- match(elements[i],rownames(component))
      if(!is.na(n1)) n <- n + compound[n1,]
      if(!is.na(n2)) n <- n + component[n2,]
      nelement[i] <- n
    }
    f <- data.frame(count=nelement)
    rownames(f) <- elements
    return(f)
  } 
  
  # if first argument is numeric, get the corresponding species formulas
  if(is.numeric(compound[1])) compound <- as.character(thermo$obigt$formula[compound])

  # if first argument is of length > 1, return the makeup of
  # that (those) species, multiplied by the coefficients
  # specified in the second argument (default=1)
  if(length(compound) > 1 & !is.character(component[1]) & !is.logical(component[1])) {
    if(missing(component)) component <- rep(1,length(compound))
    if(length(compound)==1) {
      return(component * makeup(compound))
    } else {
      m <- data.frame()
      for(i in 1:length(compound)) {
        m <- makeup(m,component[i] * makeup(compound[i]))
      }
      return(m)
    }
  }

  #if(is.character(compound) & is.numeric(component)) {
  formula.word <- function(compound,component) {
  # return the word beginning at the specified (2nd arg)
  # position of the formula (1st word)
    # a word in the chemical formula of a compound can be:
    # - the abbreviation of an element 
    # (an uppercase letter followed by zero or more lowercase letters)
    # - the coefficient on the element
    # (a positive, possibly non-integer value which may be absent)
    # - the net electronic (not protonic) charge
    # ('+' or '-', possibly followed by a possibly non-integer value)
    ###  
    # words are juxtaposed, may appear multiple times, and
    # can have any arrangement within these conditions:
    # - the first word must be an element-word
    # - up to one coefficient-word may follow a given element-word
    # - the net electronic charge, if it is present, is the last word
    ###
    # this function isolates a given word by juxtaposing the letter 
    # identified by pos with the following ones, up to the one preceding 
    # the beginning of the next word (or the end of the formula)
    ###
    # consequently, the value of pos given in the function call 
    # should correspond to the start of a word

    # what type is the current position?
    pos <- component
    if(can.be.numeric(substr(compound,pos,pos))) 
      pos.is.value <- TRUE else pos.is.value <- FALSE
    
    # what is last position of this word?
    pos.end <- nchar(compound)
    flag <- FALSE
    for(i in pos:nchar(compound)) {
      if(flag) next()
      letter <- substr(compound,i,i)
      if(can.be.numeric(letter) != pos.is.value & letter!='.') {pos.end <- i-1; flag <- TRUE}
      if(!can.be.numeric(letter)) {
        if(letter==toupper(letter) & i != pos) {pos.end <- i-1; flag <- TRUE}
      }
      if(letter %in% c('+','-') & i != pos) {pos.end <- i-1; flag <- TRUE}
    }

    # get and return our word, correctly typed
    word <- substr(compound,pos,pos.end)
    if(word=='+') return(list(word=1,length=1))
    if(word=='-') return(list(word=-1,length=1))
    if(pos.is.value) return(list(word=as.numeric(word),length=nchar(word)))
    else return(list(word=word,length=nchar(word)))
  }

  # if the first argument is character and the second is missing
  # return the makeup -- a dataframe with the elemental coefficients
  # of the formula (with names of elements on rows)
  # call the formula string fcomp
  if(is.character(compound) & missing(component)) {
    fcomp <- compound
    # split formula at asterisk or colon
    # as in K2(Al2Si5)O14*5H2O, or C6H14N2O2:HCl
    makeup.dot <- data.frame()
    fcomp.c <- s2c(fcomp)
    is.dot <- fcomp.c == '*' | fcomp.c == ':'
    if(TRUE %in% is.dot) {
      idot <- match(TRUE,is.dot)
      notdotcomp <- substr(fcomp,1,idot-1)
      dotcomp <- substr(fcomp,idot+1,nchar(fcomp))
      fcomp <- notdotcomp
      # does the dot stuff have a coefficient
      dotcoeff <- 1
      dotfirstword <- formula.word(dotcomp,1)
      if(can.be.numeric(dotfirstword$word)) {
        dotcoeff <- dotfirstword$word
        dotcomp <- substr(dotcomp,1+dotfirstword$length,nchar(dotcomp))
      }
      makeup.dot <- makeup(dotcomp) * dotcoeff
    }
    # remove parentheses, if necessary multiplying the stuff
    # inside by some coefficient
    # look for a closing paren (grep doesn't like an opening one)
    makeup.par <- data.frame()
    if(length(grep(')',fcomp))>0) {
      par.ends <- which(fcomp.c == ')')
      par.starts <- which(fcomp.c == '(')
      # we can deal with multiple pairs (but not nested) parenthesis
      # unpaired parens will give wacky results
      for(i in 1:length(par.ends)) { 
        fcomp.c <- s2c(fcomp)
        # extract the parenthetical and (maybe) coefficient
        par.end <- par.ends[i]
        par.start <- par.starts[i]
        fcomp.par <- substr(fcomp,par.start+1,par.end-1)
        coeff <- formula.word(fcomp,par.end+1)
        par.coeff <- 1
        if(can.be.numeric(coeff$word) & !( substr(fcomp,par.end,par.end) %in% c('-','+') ) ) {
          if(!substr(fcomp,par.end+1,par.end+1) %in% c('+','-')) {
            par.end <- par.end + coeff$length
            par.coeff <- coeff$word
          }
        } 
        makeup.par <- makeup(makeup.par,(makeup(fcomp.par)*par.coeff))
        # rewrite the formula, without parenthetical
        fcomp.start <- fcomp.end <- ''
        if( (par.start) > 1 ) fcomp.start <- substr(fcomp,1,par.start-1)
        if( (par.end) < nchar(fcomp) ) fcomp.end <- substr(fcomp,par.end+1,nchar(fcomp))
        # a workaround for compounds with only charge after substituting parens
        if(can.be.numeric(fcomp.end) & fcomp.start=='') fcomp.start <- 'C0'
        fcomp <- paste(fcomp.start,fcomp.end,sep='')
        # recalculate our positions
        par.starts <- par.starts - (par.end - par.start + 1)
        par.ends <- par.ends - (par.end - par.start + 1)
      }
    }
    # quick calc: NULL value
    if(identical(fcomp,'')) {
      e <- data.frame()
    } else {
      # count the elements
      elements <- character(0); count <- numeric(0)
      ie <- 0; j <- 0; Z <- 0
      for(i in 1:nchar(fcomp)) {
        if(j > 0) {j <- j - 1; next()}
        letter <- substr(fcomp,i,i)
        # get the word
        cw <- formula.word(fcomp,i)
        if(!can.be.numeric(letter)) {
          ie <- ie + 1; elements[ie] <- cw$word
          j <- cw$length - 1; count[ie] <- 1
        }
        # 20070923 charge only if it's the last word
        if(i+cw$length-1==nchar(fcomp)) can.be.charge <- TRUE else can.be.charge <- FALSE
        # 20071128 and if it's not preceded by 'Z'
        if(can.be.charge) if(formula.word(fcomp,i-1)$word=='Z') can.be.charge <- FALSE
        if(letter %in% c('+','-') & can.be.charge) {
          Z <- cw$word; j <- cw$length - 1
        } else {
          if(can.be.numeric(letter)) {
            count[ie] <- cw$word; j <- cw$length - 1
          }
        }
      }
      # build the formula by successively summing it
      e <- data.frame()
      for(i in 1:ie) {
        f <- data.frame(count=count[i])
        rownames(f) <- elements[i]
        # it seems like these should be counted
        #if(rownames(f)=='Z') f$count <- 0 
        e <- makeup(e,f)
      }
      # charge as an element
      if(Z!=0) {
        z <- data.frame(count=Z)
        rownames(z) <- 'Z'
        e <- makeup(e,z)
      }
    }
    # if there were parenthesis or dots, add those compositions
    if(length(makeup.dot)>0) e <- makeup(e,makeup.dot)
    if(length(makeup.par)>0) e <- makeup(e,makeup.par)
    # complain if there are any lower-case element symbols
    ilc <- !(rownames(e) %in% thermo$element$element)
    if(length(which(ilc))>1) p <- 's' else p <- ''
    if(length(which(ilc))>1) p2 <- 'are not recognized elements.' 
      else p2 <- 'is not a recognized element.'
    if(any(ilc)) warning(paste('\'',
      c2s(rownames(e)[ilc],sep='\' \''),'\' ',p2,sep=''))
    return(e)
  }

  # if component is logical, return the basis,
  # or eliminate small numbers in the given basis
  if(is.logical(component)) {
    # FALSE: cutoff to eliminate small numbers
    if(!component) {
      for(j in 1:ncol(compound)) if(abs(compound[,j]) < thermo$opt$cutoff) compound[,j] <- 0
      # this might also be the best place to make sure we keep
      # pluses and minuses in the colnames (anti-R behavior)
      #colnames(compound) <- rownames(thermo$basis)[1:nrow(thermo$basis)]
      return(compound)
    # TRUE: return the basis (coefficients of the formation reaction)
    } else {
      if(is.null(thermo$basis)) stop("makeup: missing basis species; load with basis() function")
      basis.elements <- colnames(basis())
      if(!is.data.frame(compound)) compound.makeup <- makeup(compound) else compound.makeup <- compound
      compound.elements <- makeup(compound.makeup,'')
      okayelements <- TRUE
      iselement <- colnames(compound.elements) %in% basis.elements
      for(i in 1:length(iselement)) {
        if(!iselement[i]) {
          cat(paste('basis:',colnames(compound.elements)[i],'of',compound,'is not contained by the basis species\n'))
          okayelements <- FALSE
        }
      }
      if(!okayelements) stop('makeup: one or more elements not contained by the basis species.')
      nelements <- numeric(0)
      for(i in 1:length(basis.elements)) {
        ielement <- match(basis.elements[i],colnames(compound.elements))
        if(!is.na(ielement)) nelements[i] <- compound.elements[,ielement] 
        else nelements[i] <- 0
      }
      nbasis <- solve(t(thermo$basis[,1:nrow(thermo$basis)]),nelements)
      t <- data.frame(matrix(nbasis,nrow=1))
      colnames(t) <- rownames(thermo$basis)
      # basis species in the formation reaction: the negative of composition
      # send it through this function again to get rid of small numbers
      return(makeup(-t,FALSE))
    }
  }

  # if second argument is character,
  # return either a dataframe (elements on columns) or formula string
  # (with selected or default elements)
  #if(is.data.frame(compound) & is.character(component)) {
  if(is.character(component)) {
    # get a makeup if it appears a formula was supplied
    if(!is.data.frame(compound)) compound <- makeup(compound)
    # formula data frame to string
    #if(ncol(compound)>1) {
    if(is.data.frame(compound) & colnames(compound)[1]!='count') {
      f <- character()
      if(identical(component,'')) component <- colnames(compound)
      for(i in 1:ncol(compound)) {
        element <- colnames(compound)[i]
        #if(! element %in% component | element=='Z') {
        if(! element %in% component ) {
          element <- ''
          next
        }
        numcoeff <- compound[1,i]
        coeff <- as.character(numcoeff)
        if(i==ncol(compound) & element=='Z') {
          # discard a final Z
          element <- ''
          if(as.numeric(coeff) > 0) coeff <- paste('+',coeff,sep='')
          if(coeff=='1') coeff <- '+'
          if(coeff=='-1') coeff <- '-'
        } else {
          if(coeff=='1') coeff <- ''
          else coeff <- format(numcoeff,scientific=FALSE,digits=16)
        }
        f <- paste(f,element,coeff,sep='')
      }
      # 20081103: if we only have charge, explicitly write the Z
      if(identical(colnames(compound),'Z')) f <- paste('Z',f,sep="")
      return(f)
    # makeup to formula data frame (elements)
    } else {
      # character = '': default to all elements in compound
      if(identical(component,'')) component <- rownames(compound)
      g <- as.data.frame(matrix(rep(0,length(component)),nrow=1))
      colnames(g) <- component
      #e <- makeup(compound)
      for(i in 1:length(g)) {
        h <- match(colnames(g)[i],rownames(compound))
        if(!is.na(h)) g[1,i] <- compound[h,1] else g[1,i] <- 0
      }
      # strip zeros
      if(ncol(g)>1) {
        notzero <- g[1,]!=0
        g <- as.data.frame(g[,notzero])
        colnames(g) <- rownames(compound)[notzero]
      }
      rownames(g) <- ''
      return(g)
    }
  }

  # if more than one compound specified
  # (either list of dataframes or character of length > 1),
  # return the names of the unique elements
  if(length(compound>1) & (is.list(compound) | is.character(compound))) {
    elements <- character(0)
    for(i in 1:length(compound)) {
      if(is.list(compound)) c <- compound[[i]] else c <- makeup(compound[i])
      elements <- c(elements,row.names(c))
    }
    return(unique(elements))
  }

}
  
