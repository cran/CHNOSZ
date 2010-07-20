# CHNOSZ/protein.R
# Copyright (C) 2006-2009 Jeffrey M. Dick
# calculate properties of proteins
# 20061109 jmd

protein <- function(protein,organism=NULL,online=thermo$opt$online) {

  protein.download <- function(protein,organism,online) {
    # download protein sequence information from
    # the saccharomyces genome database
    # or e. coli database or SWISS-PROT
    # or get it from local object (SGD,ECO,)
    if(is.na(online)) {
      # don't make an online search without
      # confirmation from the user
      ANSWER <- readline(paste('Shall I try an online search for ',protein,'_',organism,'? '))
      if(tolower(substr(ANSWER, 1, 1)) == "n") {
        thermo$opt$online <<- online <-  FALSE
      } else thermo$opt$online <<- online <- TRUE
    }
    if(!online) return(NA)
    iprotein <- numeric()
    if(organism=='SGD') {
      ## as of CHNOSZ_0.8 this has been superseded by the SGD
      ## amino acid compositions included in the package
      # try to get protein composition from yeastgenome.org
      url <- paste('http://db.yeastgenome.org/cgi-bin/locus.pl?locus=',protein,sep='')
      cat(paste('protein: trying ',url,'\n',sep=''))
      myt <- readLines(url)
      if(length(grep('did not return',myt[[3]])) > 0) {
        cat(paste('protein: unmatched locus name ',protein,'.\n',sep=''))
        return(NA)
      }
      else cat(paste('protein: scanning for URLs of proteinPage ... '))
      # we only need to look at some of the characters to get our URL
      url2 <- substr(myt[[3]],4500,nchar(myt[[3]]))
      url2 <- s2c(url2,sep=c('http','\"'))
      url2 <- url2[grep('http',url2)]
      url2 <- url2[grep('proteinPage',url2)]
      cat(' loading ...')
      myt <- readLines(url2)
      # where the sequence starts
      line.start <- grep('<pre>',myt)
      line.stop  <- grep('</pre>',myt)
      # build the sequence string
      s <- s2c(myt[[line.start]],sep='<pre>',keep.sep=FALSE)[[2]]
      for(i in (line.start+1):line.stop) {
        s <- c(s,myt[[i]])
      }
      s[length(s)] <- s2c(s[length(s)],sep='</pre>',keep.sep=FALSE)[[1]]
      s <- aminoacids(s)
      cat(paste(' length is ',sum(s[1,]),'.\n',sep=''))
      # 20090316 moved "chains" to column 5
      colnames(s) <- colnames(thermo$protein)[6:25]
      s <- data.frame(protein=protein,organism='SGD',source=NA,abbrv=NA,s)
      thermo$protein <<- rbind(thermo$protein,s)
      rownames(thermo$protein) <<- seq(1:nrow(thermo$protein))
      return(nrow(thermo$protein))
    }
    if(organism=='SWISS') {
      url <- paste('http://www.expasy.org/uniprot/',protein,sep='')
      cat(paste('protein: trying ',url,' ... ',sep=''))
      oldopt <- options(warn=-1)
      myt <- try(readLines(url),TRUE)
      options(oldopt)
      if(class(myt)=='try-error') {
        cat('failed.\n')
        return(NA)
      } else cat('got it!\n')
      # 20091102: the following is too sensitive to changes in 
      # page format so instead we just look for a link to a fasta file
      #tt <- myt[[4]]
      #if(tt=='<head>') tt <- myt[[5]]
      #accession.number <- s2c(tt,sep=c('entry',' '),keep.sep=FALSE)[[3]]
      tt <- myt[[grep("/uniprot/.*fasta",myt)[1]]]
      # also, s2c is slow so let's use strsplit
      tt <- strsplit(tt,".fasta",fixed=TRUE)[[1]][1]
      accession.number <- tail(strsplit(tt,"/uniprot/",fixed=TRUE)[[1]],1)
      cat(paste('protein: found ',accession.number,' ... ',sep=''))
      url <- paste('http://www.uniprot.org/uniprot/',accession.number,'.fasta',sep='')
      tt <- readLines(url)
      cat(paste(s2c(tt[[1]],sep=c(protein,'(EC','OS'),keep.sep=FALSE)[[2]],sep=''))
      # get rid of the header before counting letters
      tt[[1]] <- ""
      # 20080919 the following does not work, as we later get
      # Error: 'getEncChar' must be called on a CHARSXP
      #tt[[1]] <- character()
      s <- aminoacids(c2s(tt))
      cat(paste(' (length ',sum(s[1,]),').\n',sep=''))
      colnames(s) <- colnames(thermo$protein)[6:25]
      po <- s2c(protein,sep='_',keep.sep=FALSE)
      s <- data.frame(protein=po[1],organism=po[2],source=NA,abbrv=NA,chains=1,s)
      thermo$protein <<- rbind(thermo$protein,s)
      rownames(thermo$protein) <<- seq(1:nrow(thermo$protein))
      return(nrow(thermo$protein))
    }
    return(iprotein)
  }  # end protein.download

  # if first arg. is a data frame, sum that
  # composition into a new protein
  if(is.data.frame(protein[1])) {
    if(!is.null(organism)) {
      pname <- s2c(organism,sep='_',keep.sep=FALSE)
      if(identical(pname,organism)) stop('invalid name for new protein (needs underscore)')
    } else stop('missing name for new protein')
    newprotein <- protein[1,]
    newprotein[5:25] <- colSums(protein[5:25])
    newprotein$protein <- pname[1]
    newprotein$organism <- pname[2]
    inew <- protein(pname[1],pname[2],online=FALSE)
    if(length(inew)==0) {
      thermo$protein <<- rbind(thermo$protein,newprotein)
      inew <- nrow(thermo$protein)
    }
    else thermo$protein[inew,] <<- newprotein
    return(inew)
  }
  # if first arg. looks like protein name
  if(is.character(protein[1])) {
    if(length(grep('_',protein[1])>0)) {
      # return rows of thermo$protein
      iprotein <- numeric()
      for(i in 1:length(protein)) {
        po <- s2c(protein[i],sep='_',keep.sep=FALSE)
        inew <- which(thermo$protein$protein==po[[1]] & thermo$protein$organism==po[[2]])
        if(length(inew)==0) {
          # maybe we'll find it online
          #if(po[[2]]=="SGD") inewprotein <- protein.download(po[[1]],po[[2]],online=online)
          inewprotein <- protein.download(protein[i],"SWISS",online=online)
          if(!is.na(inewprotein)) inew <- inewprotein
        }
        iprotein <- c(iprotein,inew)
      }
      return(thermo$protein[iprotein,])
    }
  }
  # if second arg. looks like protein name...
  # then we think the first argument is a sequence
  if(length(grep('_',organism)>0)) {
    t <- aminoacids(protein)
    if(!FALSE %in% (as.numeric(t)==0)) stop('those do not look like amino acids.')
    po <- s2c(organism,sep='_',keep.sep=FALSE)
    pname <- po[1]; oname <- po[2]
    if(TRUE %in% (thermo$protein$protein==pname & thermo$protein$organism==oname))
      stop('a protein with that name (',organism,') already exists.')
    newrow <- data.frame(protein=pname,organism=oname,source=NA,abbrv=NA)
    newrow <- cbind(newrow,chains=1,t)
    colnames(newrow) <- colnames(thermo$protein)
    thermo$protein <<- rbind(thermo$protein,newrow)
    cat(paste('protein: added ',pname,'_',oname,' (length=',protein.length(organism),').\n',sep=''))
    return(thermo$protein[nrow(thermo$protein),])
  }
  # argument handling: make organism and protein the same length
  if(!missing(organism))
    if(length(organism) > length(protein)) protein <- rep(protein,length.out=length(organism))
  else if(length(protein) > length(organism)) organism <- rep(organism,length.out=length(protein))
    # numeric first argument: display and return the properties of
    # proteins calculated from amino acid composition
    # character second argument denotes the state
    if(is.numeric(protein[1])) {
      for(j in 1:length(protein)) {
        # get the composition
        seq <- as.numeric(thermo$protein[protein[j],6:25])
        chains <- as.numeric(thermo$protein$chains[protein[j]])
        groups <- c(colnames(thermo$protein)[6:25],'AABB','UPBB')
        length <- sum(seq)
        # seq, AABB, UPBB
        seq <- c(seq,chains,length-chains)  # n[UPBB] is n-1 if chains=1
        if(!missing(organism)) state <- organism else state <- thermo$opt$state
        # limit to state
        to <- thermo$obigt[thermo$obigt$state==state,]
        # put brackets around amino acid sidechain group names
        bracket <- function(x) paste('[',x,']',sep='')
        aa <- sapply(groups,bracket)
        to <- to[match(aa,to$name),]
        tf <- as.character(to$formula)
        to <- to[,8:20]
        # the actual adding and multiplying of thermodynamic properties
        # hmm. seems like we have to split up the multiplication/transposition
        # operations to get the result into multiple columns. 20071213
        to <- colSums(to * seq)
        to <- as.data.frame(t(to))
        # *don't* convert units by scaling factors, but rename columns
        if(state=='aq') {
          #to[,6:12] <- to[,6:12] * 10^c(-1,2,0,4,0,4,0,5)
          colnames(to)[6:13] <- c('a1','a2','a3','a4','c1','c2','omega','Z')
        } else {
          #to[,6:12] <- to[,6:11] * 10^c(0,-3,5,0,-5,0,0)
          colnames(to)[6:13] <- c('a','b','c','d','e','f','lambda','T')
        }
        # 20090331 replaced loop formula construction with following...
        # to get the formula, add up and round the group compositions
        f.in <- round(as.data.frame(t(colSums(seq*t(thermo$groups)))),2)
        # take out any elements that don't appear (sometimes S)
        f.in <- f.in[,f.in!=0]
        # turn it into a formula
        f <- makeup(f.in,"")
        name <- paste(thermo$protein$protein[protein[j]],'_',thermo$protein$organism[protein[j]],sep='')
        # make some noise
        cat(paste('protein: found '))
        cat(paste(name,' (',f,', ',sep=''))
        cat(paste(round(as.numeric(protein.length(name)),3),' residues)\n',sep=''))
        source <- thermo$protein$source[protein[j]]
        to <- data.frame(name=name,abbrv=NA,formula=f,state=state,source1=source,source2=NA,date=NA,to)
        if(j==1) jo <- to else jo <- rbind(jo,to)
      }
      return(jo)
    }
    # character arguments: display the amino acid composition, and
    # return the indices of the proteins in thermo$protein
    if(is.character(protein[1])) {
        if(is.null(organism[1])) stop('protein: please supply an organism')
        iprotein <- numeric()
        for(i in 1:length(protein)) {
          ip <- which(thermo$protein$protein==protein[i] & 
            thermo$protein$organism==organism[i])[1]
          if(is.na(ip)) {
            if(organism[i] %in% c('SGD','ECO')) ip <- protein.download(protein[i],organism[i],online=online)
            else ip <- protein.download(paste(protein[i],'_',organism[i],sep=''),
              organism='SWISS',online=online)
          } 
          if(is.null(ip)) ip <- NA
          if(is.na(ip)) if(sys.nframe() < 2) cat(paste('protein: no match for ',protein[i],'_',organism[i],'\n',sep=''))
          iprotein <- c(iprotein,ip)
        }
        # strip NA results
        iprotein <- iprotein[!is.na(iprotein)]
        return(invisible(iprotein))
    }
}

yeastgfp <- function(location,exclusive=TRUE) {
  # return a list of ORFs and protein abundances
  # for the given subcellular location
  ygfp <- thermo$yeastgfp
  ncol <- match(location,colnames(ygfp)[6:28]) + 5
  if(is.na(ncol)) ncol <- agrep(location,colnames(ygfp)[6:28])[1] + 5
  if(is.na(ncol)) stop(paste(location,'not matched in yeastgfp.'))
  do.me <- ygfp[,ncol]
  if(exclusive) {
    # find the number of localizations of each ORF
    localizations <- numeric(nrow(ygfp))
    for(i in 6:28) localizations <- localizations + as.logical(ygfp[,i])
    if(all(localizations[do.me] > 1)) warning(paste("yeastgfp: no exclusive localization found for ",location,
      " ... using non-exclusive localizations\n",sep=""))
    else do.me <- do.me & ! localizations > 1
  }
  orf <- as.character(ygfp$yORF[do.me])
  abundance <- ygfp$abundance[do.me]
  return(list(yORF=orf,abundance=abundance))
}

get.protein <- function(protein,organism,abundance=NULL,pname=NULL,average=TRUE,digits=1) {
  # a replacement for the 'proteome' function 20090311
  # return the composition of one or more proteins: 
  # from E. coli or S. cerevisiae (organism = ECO or SGD)
  # or representing various stress response experiments 
  # (organism in colnames of thermo$stress)
  ## step 0: call ourself to get data for multiple stress experiments
  if(length(protein) > 1 & all(protein %in% colnames(thermo$stress)) & !is.null(abundance)) {
    for(i in 1:length(protein)) {
      gp <- get.protein(protein[i],organism,abundance,pname,average)
      if(i==1) out <- gp else out <- rbind(out,gp)
    }
    return(out)
  }
  ## step 1: organism should be "SGD" or "ECO" (length 1)
  if(length(organism) > 1) {
    cat("get.protein: only using first value for 'organism'\n")
    organism <- organism[1]
  }
  ## step 2: get protein if we're looking at a stress experiment
  if(protein[1] %in% colnames(thermo$stress)) {
    cat(paste("get.protein: looking up proteins for",protein[1],organism))
    imatch <- which(colnames(thermo$stress) == protein[1] & thermo$stress[1,] == organism)
    if(length(imatch)==0) {
      cat("\n")
      stop(paste("didn't find proteins for",organism))
    }
    # the source of the data
    cat(paste(" from",thermo$stress[2,imatch],"\n"))
    if(is.null(pname)) pname <- protein[1]
    protein <- thermo$stress[3:nrow(thermo$stress),imatch]
    protein <- protein[!is.na(protein)]
  }
  ## step 2.1: make abundance same length as number of proteins
  abundance <- rep(abundance,length.out=length(protein))
  ## step 2.5: make the names unique
  length.in <- length(protein)
  protein <- unique(protein)
  length.unique <- length(protein)
  ## step 3: retrieve protein sequences from SGD or ECO dataframes
  im <- match(organism,names(thermo))
  if(is.na(im)) stop(paste("missing thermo$",organism," data frame",sep=""))
  mydata <- thermo[[im]]
  if(organism=="SGD") {
    # which columns to search for matches
    searchcols <- c("OLN","OLN")
    # which columns have the amino acids in the order of thermo$protein 
    iaa <- c(1,5,4,7,14,8,9,10,12,11,13,3,15,6,2,16,17,20,18,19) + 2
  } else if(organism=="ECO") {
    searchcols <- c("OLN","Name")
    iaa <- 1:20 + 4
  }
  icols <- match(searchcols,colnames(mydata))
  # find the matches
  imatch <- match(protein,mydata[,icols[1]])
  imatch2 <- match(protein,mydata[,icols[2]])
  imatch[!is.na(imatch2)] <- imatch2[!is.na(imatch2)]
  # report and remember the unsuccessful matches
  inotmatch <- which(is.na(imatch)) 
  if(length(inotmatch) > 0) {
    if(length(inotmatch)==1) verb <- "was" else verb <- "were"
    cat(paste("get.protein:",c2s(protein[inotmatch]),verb,"not matched\n"))
    imatch <- imatch[-inotmatch]
    protein <- protein[-inotmatch]
    abundance <- abundance[-inotmatch]
  }
  if(length.unique==length.in) cat(paste("get.protein: found ",length(protein)," of ",length.in," proteins",sep=""))
  else cat(paste("get.protein: found ",length(protein)," of (",length.unique," unique of ",length.in," proteins)",sep=""))
  if(length(imatch)==0) stop("no proteins found!")
  aa <- data.frame(mydata[imatch,iaa])
  ## step 4: sum up the compositions if abundance(s) are given
  if(!is.null(abundance)) {
    # strip NAs
    ina <- is.na(abundance)
    abundance <- abundance[!ina]
    aa <- aa[!ina,]
    cat(paste(" with abundances for",length(abundance),"\n"))
    # the name of the new protein
    if(is.null(pname)) stop("please give pname")
    protein <- pname
    if(average) aa <- data.frame(t(round(colSums(aa*abundance)/sum(abundance),digits)))
    else aa <- data.frame(t(colSums(aa*abundance)))
  } else cat("\n")
  ## step 5: the identifying columns
  organism <- rep(organism,length(protein))
  abbrv <- rep(NA,length(protein))
  source <- rep(date(),length(protein))
  chains <- rep(1,length(protein))
  precols <- data.frame(protein,organism,source,abbrv,chains,stringsAsFactors=FALSE)
  colnames(aa) <- aminoacids(nchar=3)
  aa <- cbind(precols,aa)
  ## step 6: return
  return(aa)
}

protein.residue <- function(proteins) {
  # make new proteins corresponding to residues
  # of the original proteins
  matches.protein <- function(name) {
    n <- s2c(name,sep='_',keep.sep=FALSE)
    p <- n[1]
    o <- n[2]
    return(thermo$protein$protein==p & thermo$protein$organism==o)
  }
  for(i in 1:length(proteins)) {
    myname <- paste(proteins[i],'RESIDUE',sep='.')
    if(!any(matches.protein(proteins[i]))) {
      cat(paste('protein.residue: no such protein',proteins[i],'\n'))
      next
    }
    p.row <- thermo$protein[matches.protein(proteins[i]),]
    n <- s2c(myname,sep='_',keep.sep=FALSE)
    p <- n[1]; o <- n[2]
    p.row$protein <- p
    p.row$organism <- o
    p.row[5:25] <- p.row[5:25]/sum(p.row[6:25])
    if(i==1) p.out <- p.row else p.out <- rbind(p.out,p.row)
  }
  return(p.out)
}

add.protein <- function(file="protein.csv") {
  # add entries to thermo$protein
  tp1 <- thermo$protein
  # can also supply a dataframe (e.g. made using get.protein)
  ftext <- ""
  if(is.data.frame(file)) tp2 <- file else {
    # if its a fasta file, read the sequences
    if(is.fasta(file)) {
      p <- read.fasta(file)
      cat(paste("add.protein: first line in fasta file is\n"))
      cat(readLines(file,n=1))
      cat("\n")
      # we recurse, using the data frame this time
      ip <- add.protein(p)
      return(ip)
    }
    # 20090428 added colClasses here
    tp2 <- read.csv(file,colClasses=c(rep("character",4),rep("numeric",21)))
    ftext <- paste("from",file)
  }
  tp1.names <- paste(tp1$protein,tp1$organism,sep="_")
  tp2.names <- paste(tp2$protein,tp2$organism,sep="_")
  iadd <- which(!tp2.names %in% tp1.names)
  if(length(iadd) > 0) tp1 <- rbind(tp1,tp2[iadd,])
  rownames(tp1) <- 1:nrow(tp1)
  thermo$protein <<- tp1
  if(length(iadd)==1) {
    cat(paste("add.protein: added ",length(iadd)," protein ",ftext,"with length ",
      protein.length(paste(tp2$protein[iadd],tp2$organism[iadd],sep="_")),"\n",sep=""))
  } else { 
    cat(paste("add.protein: added",length(iadd),"of",nrow(tp2),"proteins",ftext,"\n"))
  }
  # return the indices of the proteins
  tp1.names <- paste(tp1$protein,tp1$organism,sep="_")
  return(invisible(match(tp2.names,tp1.names)))
}

protein.info <- function(T=25) {
  # make a table of selected properties for proteins
  # what proteins are present?
  ip <- grep("_",species()$name)
  if(length(ip)==0) stop("no proteins are defined")
  protein <- character()
  length <- numeric()
  formula <- character()
  G <- numeric()
  Z <- numeric()
  G.Z <- numeric()
  # nominal carbon oxidation state
  ZC <- numeric()
  
  # add the ionizable groups
  if("H+" %in% colnames(thermo$species)) {
    is.species <- thermo$species$ispecies[ip]
    is.ion <- ionize(affinity=NULL)
    is.only.ion <- is.ion[!is.ion %in% is.species]
    a <- affinity(T=T)
    this.Z <- ionize(a$values,a$values)
    ag <- affinity(property='G.species',T=T)
    this.G.Z <- ionize(a$values,ag$values)
  } else {
    ag <- affinity(property='G.species',T=T)
    this.Z <- rep(NA,length(ip))
    this.G.Z <- rep(NA,length(ip))
  }
  for(i in 1:length(ip)) {
    p <- as.character(thermo$species$name[ip[i]])
    protein <- c(protein,p)
    length <- c(length,protein.length(p)[[1]])
    ii <- info(info(p))
    formula <- c(formula,as.character(ii$formula))
    #G <- c(G,ii$G)
    # 20090827 use G values from G.species above (at T) 
    G <- c(G,ag$values[[ip[i]]])
    Z <- c(Z,this.Z[[ip[i]]])
    G.Z <- c(G.Z,this.G.Z[[ip[i]]])
    # make sure all the elements we want are listed (even zeros)
    m <- as.data.frame(t(makeup(c(formula[i],'C0H0N0O0S0'))))
    ZC <- c(ZC,(-1*m$H+3*m$N+2*m$O+2*m$S)/m$C)
  }
  # remove ionizable groups
  if("H+" %in% colnames(thermo$species)) species(is.only.ion,delete=TRUE)
  cat('protein.info: converting things ...\n')
  length <- round(length,1)
  Z <- round(Z,2)
  ZC <- round(ZC,3)
  G <- round(G/1000,2)
  G.Z <- round(G.Z/1000,2)
  for(i in 1:length(formula)) {
    f <- makeup(formula[i])
    formula[i] <- makeup(makeup(round(f[order(rownames(f)),,drop=FALSE]),''),'')
  }
  return(data.frame(protein=protein,length=length,formula=formula,G=G,Z=Z,G.Z=G.Z,ZC=ZC))
}

residue.info <- function(T=25) {
  # 20090902 calculate the per-residue composition
  # of proteins that are in the species list
  # which species are proteins
  ip <- grep("_",thermo$species$name)
  if(length(ip)==0) stop("no proteins are defined")
  # grab the species definitions
  ts <- thermo$species[ip,]
  # calculate ionization states if H+ is a basis species
  if("H+" %in% colnames(thermo$species)) {
    pH <- -thermo$basis$logact[rownames(thermo$basis)=="H+"]
    # load ionizable residues
    iz <- ionize()
    # calculate ionization states
    a <- affinity(T=T)
    z <- as.numeric(ionize(a$values,a$values)[ip])
    # save the results
    ts[,colnames(ts)=="H+"] <- z
    # unload ionizable groups
    species(iz,delete=TRUE)
  }
  # get lengths of proteins
  pl <- protein.length(ts$name)
  # calculate numbers of components per residue
  nb <- nrow(thermo$basis)
  for(i in 1:nb) ts[,i] <- ts[,i]/pl
  # return the result
  return(ts[,c(1:nb,ncol(ts))])
}

