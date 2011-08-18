# CHNOSZ/util.data.R
# add or change entries in the thermodynamic database

mod.obigt <- function(species,...,missingvalues=NA) {
  # add or modify species in thermo$obigt
  args <- list(...)
  for(i in 1:length(args)) {
    args[[i]] <- rep(args[[i]],length(species))
  }
  # a function to write dates in a specific format
  mydate <- function() {
    t <- date()
    tt <- s2c(t,sep=" ",keep.sep=FALSE)
    # for single-digit days there is an extra space
    tt <- tt[!tt==""]
    tday <- tt[3]
    tmonth <- tt[2]
    tyear <- substr(tt[5],start=3,stop=4)
    return(paste(tday,tmonth,tyear,sep='.'))
  }
  inew <- numeric()
  for(i in 1:length(species)) {
    is <- NULL
    sp <- species[i]
    if(is.factor(sp)) sp <- as.character(sp)
    if(!is.numeric(sp)) {
      if('state' %in% names(args)) ii <- info(sp,args$state,quiet=TRUE,return.approx=FALSE)
      else ii <- info(sp,quiet=TRUE,return.approx=FALSE)
    } else ii <- sp
    # 20090203 use return.approx in info calls above, so
    # the test is if length of return is zero
    #if(!is.na(ii) & !is.list(ii)) {
    if(length(ii) != 0) {
      #is <- which(species[i]==thermo$obigt$name)
      is <- ii
      if('state' %in% names(args)) mystate <- args$state[i] else mystate <- thermo$opt$state
      if(mystate %in% thermo$obigt$state[is]) is <- is[match(mystate,thermo$obigt$state[is])]
      else {if('state' %in% names(args)) is <- NULL else is <- is[1]}
    }
    if(!is.null(is)) {
      # to modify a species
      newrow <- thermo$obigt[is,]
      #return(newrow)
      for(j in 1:ncol(newrow)) {
        cnames <- s2c(colnames(newrow)[j],sep='.',keep.sep=FALSE)
        if(any(tolower(cnames) %in% tolower(names(args)))) {
          # use a provided value
          newrow[,j] <- args[[which(tolower(names(args)) %in% tolower(cnames))]][i]
        } else {
          # use default values
          if(any(cnames %in% c('name','formula'))) next
          #else if(is.na(missingvalues)) newrow[,j] <- missingvalues
          #else if(!missingvalues=='') newrow[,j] <- missingvalues
          if(!missing(missingvalues) & !any(cnames %in% c('state','ref1','ref2'))) 
            newrow[,j] <- missingvalues
          if(any(cnames %in% 'ref1')) newrow[,j] <- 'USER'
          if(any(cnames %in% 'ref2')) newrow[,j] <- NA
        }
      }
      newrow$date <- mydate()
      r1 <- as.character(newrow)
      r2 <- as.character(thermo$obigt[is,])
      r2[is.na(r2)] <- 'NA'
      r1[is.na(r1)] <- 'NA'
      if(!identical(r1,r2)) {
        cat(paste('mod.obigt: updating ',newrow$name,' ',newrow$state,' (',is,').\n',sep=''))
      } else cat(paste('mod.obigt: no change for ',newrow$name,' ',newrow$state,' (',is,').\n',sep=''))
      thermo$obigt[is,] <<- newrow
      inew <- c(inew,is)
    } else {
      # add a species
      newrow <- thermo$obigt[1,]
      for(j in 1:ncol(newrow)) {
        cnames <- s2c(colnames(newrow)[j],sep='.',keep.sep=FALSE)
        if(any(cnames %in% names(args))) {
          # use a provided value
          newrow[,j] <- args[[which(names(args) %in% cnames)]][i]
        } else {
          # use default values
          if(any(cnames %in% c('name'))) newrow[,j] <- sp
          else if(any(cnames=='date')) newrow[,j] <- mydate()
          else if(is.na(missingvalues) & colnames(newrow)[j]=='state')
            newrow[,j] <- thermo$opt$state
          else if(!is.na(missingvalues)) newrow[,j] <- missingvalues
          #else if(!any(cnames=='formula')) newrow[,j] <- missingvalues
          else if(any(cnames=='ref1')) newrow[,j] <- 'USER'
          else newrow[,j] <- NA
        }
      }
      if(is.na(newrow$formula)) warning(paste('mod.obigt: formula of ',newrow$name,
        ' ',newrow$state,' is NA.',sep=''),call.=FALSE)
      cat(paste('mod.obigt: adding ',newrow$name,' ',newrow$state,' (',nrow(thermo$obigt)+1,').\n',sep=''))
      thermo$obigt <<- rbind(thermo$obigt,newrow)
      inew <- c(inew,nrow(thermo$obigt))
    }
  }
  return(invisible(inew))
}

change <- function(name,...) {
  # a wrapper for mod.obigt and mod.buffer
  if(substr(name[1],1,1)=='_') {
    name <- substr(name,2,nchar(name))
    return(mod.buffer(name,...))
  } else {
    return(mod.obigt(species=name,...))
  }
}

add.obigt <- function(file=system.file("extdata/thermo/OBIGT-2.csv",package="CHNOSZ"),
  force=FALSE,E.units="cal") {
  # add/replace entries in thermo$obigt from values saved in a file
  # only replace if force==TRUE
  if(missing(file)) {
    # we use force=TRUE for the default data file
    if(missing(force)) force <- TRUE
  }
  to1 <- thermo$obigt
  id1 <- paste(to1$name,to1$state)
  to2 <- read.csv(file,as.is=TRUE)
  id2 <- paste(to2$name,to2$state)
  # check if the file is compatible with thermo$obigt
  tr <- try(rbind(to1,to2),silent=TRUE)
  if(identical(class(tr),'try-error')) stop(paste(file,"is not compatible with thermo$obigt data table."))
  # identify duplicated
  idup1 <- which(id1 %in% id2)
  idup2 <- which(id2 %in% id1)
  ndup <- length(idup2)
  nnew <- nrow(to2) - ndup
  iadd <- 1:nrow(to2)
  # convert from J if necessary
  if(tolower(E.units)=="j") {
    # loop over each row
    for(i in 1:nrow(to2)) {
      # GHS and EOS parameters
      # don't touch volume (in column 12)
      icol <- (8:18)[-12]
      # if it's aqueous, also include omega
      if(to2$state[i]=="aq") icol <-( 8:19)[-12]
      # convert to calories
      to2[i,icol] <- convert(to2[i,icol],"cal")
    }
  }
  if(force) {
    # drop entries from original
    if(length(idup1) > 0) to1 <- to1[-idup1,]
  } else {
    if(length(idup2) > 0) iadd <- iadd[-idup2]
    ndup <- 0
  }
  inew <- numeric()
  if(length(iadd) > 0) {
    inew <- nrow(to1) + 1:length(iadd)
    to1 <- rbind(to1,to2[iadd,])
  }
  thermo$obigt <<- to1
  cat(paste("add.obigt: from",file,"added",length(iadd),"of",nrow(to2),
    "species","(",ndup,"replacements,",nnew,"new, units =",E.units,")\n"))
  cat("add.obigt: to restore default database, use data(thermo)\n")
  return(invisible(inew))
}

browse.refs <- function(key=NULL) {
  # browse to web page associated with a given source
  # of thermodynamic data. first version: 20110615
  # 'key' can be
  # NULL: show a table of all sources in a browser
  # character: open a web page for each listed source
  # numeric: open one or two web pages for each listed species
  # list: the output of subcrt()

  # retrieve the sources table
  x <- thermo$refs
  
  if(is.null(key)) {
    # create the html links
    cite <- x$citation
    x$citation <- sprintf("<a href='%s' target='_blank'>%s</a>", x$URL, cite)
    notlinked <- x$URL=="" | is.na(x$URL)
    x$citation[notlinked] <- cite[notlinked]
    # remove the last (URL) component
    #x$URL <- NULL
    x <- x[1:4]
    # count the times each source is listed in OBIGT.csv
    ns1 <- sapply(x$key, function(x) length(which(thermo$obigt$ref1==x)) )
    ns1.2 <- sapply(x$key, function(x) length(which(thermo$obigt$ref2==x)) )
    ns1 <- ns1 + ns1.2
    ns1[ns1==0] <- ""
    # count the times each source is listed in OBIGT-2.csv
    o2 <- read.csv(system.file("extdata/thermo/OBIGT-2.csv", package = "CHNOSZ"))
    ns2 <- sapply(x$key, function(x) length(which(o2$ref1==x)) )
    ns2.2 <- sapply(x$key, function(x) length(which(o2$ref2==x)) )
    ns2 <- ns2 + ns2.2
    ns2[ns2==0] <- ""
    # count the times each source is listed in protein.csv
    npr <- sapply(x$key, function(x) length(which(thermo$protein$ref==x)) )
    npr[npr==0] <- ""
    # count the times each source is listed in stress.csv
    nst <- sapply(x$key, function(x) length(which(thermo$stress[2,]==x)) )
    nst[nst==0] <- ""
    # append the counts to the table to be shown
    x <- c(x,list(ns1=ns1,ns2=ns2,npr=npr,nst=nst))
    # title to display for web page
    title <- "Sources of Thermodynamic Data in CHNOSZ"
    ### the following is adapted from print.findFn in package 'sos'
    f0 <- tempfile()
    File <- paste(f0, ".html", sep="")
    Dir <- dirname(File)
    js <- system.file("extdata/js", "sorttable.js", package = "CHNOSZ")
    file.copy(js, Dir)
    ## Sundar's original construction:
    con <- file(File, "wt")
    on.exit(close(con))
    .cat <- function(...)
      cat(..., "\n", sep = "", file = con, append = TRUE)
    ## start
    cat("<html>", file = con)
    .cat("<head>")
    .cat("<title>", title, "</title>")
    .cat("<script src=sorttable.js type='text/javascript'></script>")
    .cat("</head>")
    ### boilerplate text
    .cat("<h1>Listing of all entries in thermo$refs</h1>")
    .cat("<h3>Click on hyperlinked references to open URL in new window</h3>")
    .cat("<h3>Click on column headers to sort</h3>")
    .cat("<h3>Columns 'n..' give number of times each reference appears in data tables:</h3>")
    .cat("ns1: 'ref1' and 'ref2' in data/OBIGT.csv<br>")
    .cat("ns2: 'ref1' and 'ref2' in extdata/thermo/OBIGT-2.csv<br>")
    .cat("npr: 'ref' in data/protein.csv<br>")
    .cat("nst: second row in data/stress.csv<br><p>")
    ### start table and headers
    .cat("<table class='sortable' border='1'>\n<thead>")
    .cat("<tr>")
    .cat(sprintf("  <th>%s</th>\n</tr>",
                 paste(names(x), collapse = "</th>\n  <th>")))
    .cat("</thead>\n<tbody>")
    ### now preparing the body of the table
    paste.list <- c(lapply(x, as.character), sep = "</td>\n  <td>")
    tbody.list <- do.call("paste", paste.list)
    tbody <- sprintf("<tr>\n  <td>%s</td>\n</tr>", tbody.list)
    tbody <- sub("<td><a", "<td class=link><a", tbody, useBytes = TRUE)
    .cat(tbody)
    ### finish it!
    .cat("</tbody></table></body></html>")
    ### end adaptation from print.findFn
    # show table in browser
    browseURL(File)
    cat("browse.refs: table of references is shown in browser\n")
  } else if(is.character(key)) {
    # open the URL(s) of the given source(s)
    for(i in seq_along(key)) {
      ix <- match(key[i],x$key)
      if(is.na(ix)) {
        cat(paste("browse.refs: reference key",key[i],"not found\n"))
        next
      } 
      URL <- x$URL[ix]
      if(URL=="" | is.na(URL)) {
        cat(paste("browse.refs: no URL available for reference key",key[i],"\n"))
        next
      }
      cat(paste("browse.refs: opening URL for ",key[i]," (",x$author[ix],", ",x$year[ix],")\n",sep=""))
      browseURL(x$URL[ix])
    }
    return(invisible(URL))
  } else if(is.numeric(key)) {
    # open the URL(s) of sources associated with the indicated species
    sinfo <- info(key,quiet=TRUE)
    mysources <- unique(c(sinfo$ref1,sinfo$ref2))
    mysources <- mysources[!is.na(mysources)]
    return(browse.refs(mysources))
  } else if(is.list(key)) {
    if("species" %in% names(key)) ispecies <- key$species$ispecies
    else if("reaction" %in% names(key)) ispecies <- key$reaction$ispecies
    else stop("list does not appear to be a result from subcrt()")
    if(is.null(ispecies)) stop("list does not appear to be a result from subcrt()")
    return(browse.refs(ispecies))
  }
}

obigt2eos <- function(obigt,state,fixGHS=FALSE) {
  # remove scaling factors from EOS parameters
  # and apply column names depending on the EOS
  if(state=="aq") {
    obigt[,13:20] <- t(t(obigt[,13:20]) * 10^c(-1,2,0,4,0,4,5,0))
    colnames(obigt)[13:20] <- c('a1','a2','a3','a4','c1','c2','omega','Z') 
  } else {
    obigt[,13:20] <- t(t(obigt[,13:20]) * 10^c(0,-3,5,0,-5,0,0,0))
    colnames(obigt)[13:20] <- c('a','b','c','d','e','f','lambda','T')
  }
  if(fixGHS) {
    # fill in one of missing G, H, S
    # for use esp. by subcrt because NA for one of G, H or S 
    # will hamper calculations at high T
    # which entries are missing just one
    imiss <- which(rowSums(is.na(obigt[,8:10]))==1)
    if(length(imiss) > 0) {
      for(i in 1:length(imiss)) {
        # calculate the missing value from the others
        ii <- imiss[i]
        GHS <- GHS(as.character(obigt$formula[ii]),DG=obigt[ii,8],DH=obigt[ii,9],S=obigt[ii,10])
        icol <- which(is.na(obigt[ii,8:10]))+7
        obigt[ii,icol] <- GHS
      }
    }
  }
  return(obigt)
}

checkEOS <- function(eos, state, prop, ret.diff=FALSE) {
  # compare calculated properties from equation-of-state
  # parameters with reference (tabulated) values
  # print message and return the calculated value
  # if tolerance is exceeded
  # or NA if the difference is within the tolerance
  # 20110808 jmd
  # get calculated value based on EOS
  if(state=="aq") {
    if(prop=="Cp") {
      # value of X consistent with IAPWS95
      X <- -2.773788E-7
      # we use the value of X consistent with SUPCRT
      X <- -3.055586E-7
      refval <- eos$Cp
      calcval <- eos$c1 + eos$c2/(298.15-228)^2 + eos$omega*298.15*X
      tol <- 1
    } else if(prop=="V") {
      # value of Q consistent with IAPWS95
      Q <- 0.00002483137
      # value of Q consistent with SUPCRT92
      Q <- 0.00002775729
      refval <- eos$V
      calcval <- 41.84*eos$a1 + 41.84*eos$a2/2601 + 
        (41.84*eos$a3 + 41.84*eos$a4/2601) / (298.15-228) - Q * eos$omega
      tol <- 1
    }
  } else {
    # all other states
    if(prop=="Cp") {
      refval <- eos$Cp
      Tr <- thermo$opt$Tr
      calcval <- eos$a + eos$b*Tr + eos$c*Tr^-2 + eos$d*Tr^-0.5 + eos$e*Tr^2 + eos$f*Tr^eos$lambda
      tol <- 1
    }
  }
  # calculate the difference
  diff <- calcval - refval
  if(ret.diff) return(diff)
  else {
    # return the calculated value
    # if the difference is higher than tol
    if(!is.na(calcval)) {
      if(!is.na(refval)) {
        if(abs(diff) > tol) {
          cat(paste('checkEOS:',prop,'of', eos$name, eos$state,
            'differs by', round(diff,2), 'from tabulated value\n'))
          return(calcval)
        }
      } else return(calcval)
    }
  }
  # return NA in most cases
  return(NA)
}

checkGHS <- function(ghs,ret.diff=FALSE) {
  # compare calculated G from H and S
  # with reference (tabulated) values
  # print message and return the calculated value
  # if tolerance is exceeded
  # or NA if the difference is within the tolerance
  # 20110808 jmd
  # get calculated value based on EOS
  refval <- ghs[,8]
  #calcval <- GHS(as.character(ghs$formula),DG=ghs[,8],DH=ghs[,9],S=ghs[,10])
  Se <- element(as.character(ghs$formula),'entropy')[,1]
  DH <- ghs[,9]
  S <- ghs[,10]
  calcval <- DH - thermo$opt$Tr * (S - Se)
  tol <- 500
  # now on to the comparison
  # calculate the difference
  diff <- calcval - refval
  if(ret.diff) return(diff)
  else if(!is.na(calcval)) {
    if(!is.na(refval)) {
      diff <- calcval - refval
      if(abs(diff) > tol) {
        cat(paste('checkGHS: G of', ghs$name, ghs$state,
          'differs by', round(diff), 'from tabulated value\n'))
        return(calcval)
      }
    } else return(calcval)
  } else {
    # calculating a value of G failed, perhaps b/c of missing elements
    return(NULL)
  }
  # return NA in most cases
  return(NA)
}

check.obigt <- function() {
  # function to check self-consistency between
  # values of Cp and V vs. EOS parameters
  # and among G, H, S values
  # 20110808 jmd replaces 'check=TRUE' argument of info()
  checkfun <- function(what) {
    # looking at thermo$obigt or OBIGT-2.csv
    if(what=="OBIGT") to <- thermo$obigt
    else if(what=="OBIGT-2") {
      file <- system.file("extdata/thermo/OBIGT-2.csv",package="CHNOSZ")
      to <- read.csv(file,as.is=1:7)
    }
    ntot <- nrow(to)
    # where to keep the results
    DCp <- DV <- DG <- rep(NA,ntot)
    # first get the aqueous species
    isaq <- to$state=="aq"
    eos.aq <- obigt2eos(to[isaq,],"aq")
    DCp.aq <- checkEOS(eos.aq,"aq","Cp",ret.diff=TRUE)
    DV.aq <- checkEOS(eos.aq,"aq","V",ret.diff=TRUE)
    cat(paste("check.obigt: GHS for",sum(isaq),"aq species in",what,"\n"))
    DG.aq <- checkGHS(eos.aq,ret.diff=TRUE)
    # store the results
    DCp[isaq] <- DCp.aq
    DV[isaq] <- DV.aq
    DG[isaq] <- DG.aq
    # then other species, if they are present
    if(sum(!isaq) > 0) {
      eos.cgl <- obigt2eos(to[!isaq,],"cgl")
      DCp.cgl <- checkEOS(eos.cgl,"cgl","Cp",ret.diff=TRUE)
      cat(paste("check.obigt: GHS for",sum(!isaq),"c,g,l species in",what,"\n"))
      DG.cgl <- checkGHS(eos.cgl,ret.diff=TRUE)
      DCp[!isaq] <- DCp.cgl
      DG[!isaq] <- DG.cgl
    }
    # put it all together
    out <- data.frame(table=what,ispecies=1:ntot,name=to$name,state=to$state,DCp=DCp,DV=DV,DG=DG)
    return(out)
  }
  # check both databases in CHNOSZ
  out1 <- checkfun("OBIGT")
  out2 <- checkfun("OBIGT-2")
  out <- rbind(out1,out2)
  # set differences within a tolerance to NA
  out$DCp[out$DCp < 1] <- NA
  out$DV[out$DV < 1] <- NA
  out$DG[out$DG < 500] <- NA
  # take out species where all reported differences are NA
  ina <- is.na(out$DCp) & is.na(out$DV) & is.na(out$DG)
  out <- out[!ina,]
  # round the values
  out$DCp <- round(out$DCp,2)
  out$DV <- round(out$DV,2)
  out$DG <- round(out$DG)
  # how to make the file at extdata/thermo/obigt_check.csv
  # write.csv(out,"obigt_check.csv",na="",row.names=FALSE)
  # return the results
  return(out)
}

