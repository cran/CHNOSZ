# CHNOSZ/util.fasta.R
# read and manipulate FASTA sequence files

is.fasta <- function(file) {
  # check if the file is in FASTA format
  # read two lines in case the first one is blank
  l <- readLines(file,n=2)
  if(length(grep("^>",l)) == 0) return(FALSE) else return(TRUE)
}

grep.file <- function(file,pattern="",y=NULL,ignore.case=TRUE,startswith=">",lines=NULL,grep="grep") {
  # return the line numbers of the file that contain
  # the search term x and optionally don't contain y
  sysgrep <- function(i) {
    # 20091021 changed grep to egrep
    sysexp <- paste(mycat,' "',file,'" | ',grep,' -n ',ic,' "',startswith,pattern[i],'"  | cut -f 1 -d ":"',sep="")
    ix <- as.integer(system(sysexp,intern=TRUE))
    return(ix)
  }
  Rgrep <- function(i) {
    ix <- grep(pattern[i],lines,ignore.case=ignore.case)
    if(!is.null(y)) {
      iy <- grep(y,lines,ignore.case=ignore.case)
      ix <- ix[!ix %in% iy] 
    }
    if(!is.null(startswith)) {
      ihead <- which(substr(lines,1,1)==startswith)
      ix <- ix[ix %in% ihead]
    }
    return(ix)
  }
  # use the system's grep if available and y is NULL
  # TODO: include other *nix
  if(is.null(y) & Sys.info()[[1]]=="Linux" & is.null(lines)) {
    # figure out whether to use 'cat', 'zcat' or 'xzcat'
    suffix <- substr(file,nchar(file)-2,nchar(file))
    if(suffix==".gz") mycat <- "zcat"
    else if(suffix==".xz") mycat <- "xzcat"
    else mycat <- "cat"
    # use the system grep
    if(is.null(startswith)) startswith <- "" else startswith <- paste("^",startswith,".*",sep="")
    if(ignore.case) ic <- "-i" else ic <- ""
    out <- palply(1:length(pattern),sysgrep)
  } else {
    # use R grep
    if(is.null(lines)) lines <- readLines(file)
    out <- palply(1:length(pattern),Rgrep)
  }
  # make numeric (NA for ones that aren't matched)
  out <- as.numeric(sapply(out,as.numeric))
  return(as.numeric(out))
}

read.fasta <- function(file,i=NULL,ret="count",lines=NULL,ihead=NULL,pnff=FALSE) {
  # read sequences from a fasta file
  if(file != "") msgout("read.fasta: reading ",basename(file),"\n")
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
  # pnff: take the protein name from filename (TRUE) or entry name (FALSE)
  is.nix <- Sys.info()[[1]]=="Linux"
  if(is.nix & is.null(lines)) {
    # figure out whether to use 'cat', 'zcat' or 'xzcat'
    suffix <- substr(file,nchar(file)-2,nchar(file))
    if(suffix==".gz") mycat <- "zcat"
    else if(suffix==".xz") mycat <- "xzcat"
    else mycat <- "cat"
    nlines <- as.integer(system(paste(mycat,' "',file,'" | wc -l',sep=""),intern=TRUE))
    ihead <- as.integer(system(paste(mycat,' "',file,'" | grep -n "^>" | cut -f 1 -d ":"',sep=""),intern=TRUE))
    #linefun <- function(i1,i2) as.character(system(paste('sed -n ',i1,',',i2,'p ',file,sep=""),intern=TRUE))
    lines <- system(paste(mycat,' "',file,'"',sep=""),intern=TRUE)
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
  sequences <- palply(1:length(i), seqfun)
  # process the header line for each entry
  # (strip the ">" and go to the first space or underscore)
  nomfun <- function(befund) {
    nomnomfun <- function(j,pnff) {
      # get the text of the line
      f1 <- linefun(i[j],i[j])
      # stop if the first character is not ">"
      # or the first two charaters are "> "
      if(substr(f1,1,1)!=">" | length(grep("^> ",f1)>0))
        stop(paste("file",basename(file),"line",j,"doesn't begin with FASTA header '>'."))
      # discard the leading '>'
      f2 <- substr(f1, 2, nchar(f1))
      # keep everything before the first space
      f3 <- strsplit(f2," ")[[1]][1]
      # then before or after the first underscore
      if(befund) f4 <- strsplit(f3,"_")[[1]][1]
      else f4 <- strsplit(f3,"_")[[1]][2]
      return(f4)
    }
    noms <- as.character(palply(1:length(i),nomnomfun))
    return(noms)
  }
  # process the file name
  # (basename minus extension)
  bnf <- strsplit(basename(file),split=".",fixed=TRUE)[[1]][1]
  if(pnff) {
    # protein name is from file name
    # organism name is from entry
    protein <- bnf
    organism <- nomfun(befund=FALSE)
  } else {
    # vice versa
    protein <- nomfun(befund=TRUE)
    organism <- bnf
  }
  if(ret=="count") {
    aa <- count.aa(sequences)
    colnames(aa) <- aminoacids(3)
    ref <- abbrv <- NA
    chains <- 1
    # 20090507 made stringsAsFactors FALSE
    return(cbind(data.frame(protein=protein,organism=organism,
      ref=ref,abbrv=abbrv,chains=chains,stringsAsFactors=FALSE),aa))
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

