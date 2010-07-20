# CHNOSZ/util.fasta.R
# Copyright (C) 2009-2010 Jeffrey M. Dick
# read and manipulate FASTA sequence files

is.fasta <- function(file) {
  # check if the file is in FASTA format
  # read two lines in case the first one is blank
  l <- readLines(file,n=2)
  if(length(grep("^>",l)) == 0) return(FALSE) else return(TRUE)
}

grep.file <- function(file,x="",y=NULL,ignore.case=TRUE,startswith=">",lines=NULL,grep="grep") {
  # return the line numbers of the file that contain
  # the search term x and optionally don't contain y
  # use the system's grep if available and y is NULL
  # TODO: better test for system (include eg MacOSX etc)
  if(is.null(y) & Sys.info()[[1]]=="Linux" & is.null(lines)) {
    if(is.null(startswith)) startswith <- "" else startswith <- paste("^",startswith,".*",sep="")
    if(ignore.case) ic <- "-i" else ic <- ""
    # 20091021 changed grep to egrep
    sysexp <- paste(grep,' -n ',ic,' "',startswith,x,'" "',file,'" | cut -f 1 -d ":"',sep="")
    ix <- as.integer(system(sysexp,intern=TRUE))
  } else {
    if(is.null(lines)) lines <- readLines(file)
    ix <- grep(x,lines,ignore.case=ignore.case)
    if(!is.null(y)) {
      iy <- grep(y,lines,ignore.case=ignore.case)
      ix <- ix[!ix %in% iy] 
    }
    if(!is.null(startswith)) {
      ihead <- which(substr(lines,1,1)==startswith)
      ix <- ix[ix %in% ihead]
    }
  }
  return(ix)
}

read.fasta <- function(file,i=NULL,ret="count",lines=NULL,ihead=NULL) {
  # read sequences from a fasta file
  # all of them or only those indicated by i
  # if aa=TRUE compile a data frame of the amino acid
  # compositions suitable for add.protein
  # some of the following code was adapted from 
  # read.fasta in package seqinR
  # TODO: better test for type of system
  # value of 'ret' determines format of return value:
  # aa: amino acid composition (same columns as thermo$protein)
  # seq: amino acid sequence
  # fas: fasta entry
  is.nix <- Sys.info()[[1]]=="Linux"
  if(is.nix & is.null(lines)) {
    nlines <- as.integer(system(paste('grep -c "**" "',file,'"',sep=""),intern=TRUE))
    ihead <- as.integer(system(paste('grep -n "^>" "',file,'" | cut -f 1 -d ":"',sep=""),intern=TRUE))
    #linefun <- function(i1,i2) as.character(system(paste('sed -n ',i1,',',i2,'p ',file,sep=""),intern=TRUE))
    lines <- system(paste('cat "',file,'"',sep=""),intern=TRUE)
    linefun <- function(i1,i2) lines[i1:i2]
  } else {
    if(is.null(lines)) lines <- readLines(file)
    nlines <- length(lines)
    if(is.null(ihead)) ihead <- which(substr(lines,1,1)==">")
    linefun <- function(i1,i2) lines[i1:i2]
  }
  if(is.null(i)) {
    i <- ihead
    start <- i + 1
    end <- i - 1
    end <- c(end[-1], nlines)
  } else {
    start <- i + 1
    iend <- match(i,ihead)
    # we have to be careful about the last record
    iend[iend==ihead[length(ihead)]] <- NA
    end <- ihead[iend+1] - 1
    end[is.na(end)] <- nlines
  } 
  # just return the lines from the file
  if(ret=="fas") {
    iline <- numeric()
    for(i in 1:length(start)) iline <- c(iline,(start[i]-1):end[i])
    return(lines[iline])
  }
  seqfun <- function(i) paste(linefun(start[i],end[i]),collapse="")
  sequences <- mylapply(1:length(i), seqfun)
  nomfun <- function(j) {
    firstword <- strsplit(linefun(i[j],i[j]), " ")[[1]][1]
    substr(firstword, 2, nchar(firstword))
  }
  nomseq <- mylapply(1:length(i),nomfun) 
  if(ret=="count") {
    aa <- aminoacids(sequences)
    colnames(aa) <- aminoacids(nchar=3)
    protein <- as.character(nomseq)
    organism <- s2c(file,sep=".",keep.sep=FALSE)[1]
    source <- abbrv <- NA
    chains <- 1
    # 20090507 made stringsAsFactors FALSE
    return(cbind(data.frame(protein=protein,organism=organism,
      source=source,abbrv=abbrv,chains=chains,stringsAsFactors=FALSE),aa))
  } else return(sequences)
}

splitline <- function(line,length) {
  # to split a line into multiple lines with a specified length
  out <- character()
  count <- 0
  n <- nchar(line)
  while(count < n) {
    split <- substr(line,count+1,count+length)
    out <- c(out,split)
    count <- count + length
  }
  return(out)
}

trimfas <- function(file,start,stop) {
  # to extract certain positions from an (aligned) fasta file
  lines <- readLines(file)
  fas <- read.fasta(file="",lines=lines,ret="seq")
  # the length of lines to use
  ll <- nchar(lines[2])
  ihead <- grep("^>",lines)
  head <- lines[ihead]
  out <- character()
  for(i in 1:length(head)) {
    extract <- substr(fas[i],start,stop)
    out <- c(out,head[i],splitline(extract,ll),"")
  }
  #write.table(out,paste(file,"trim",sep="."),row.names=FALSE,col.names=FALSE,quote=FALSE)
  return(out)
}

