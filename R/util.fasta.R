# CHNOSZ/util.fasta.R
# read and manipulate FASTA sequence files

read.fasta <- function(file, iseq=NULL, ret="count", lines=NULL, ihead=NULL,
  start=NULL, stop=NULL, type="protein", id=NULL) {
  # read sequences from a fasta file
  # some of the following code was adapted from 
  # read.fasta in package seqinR
  # value of 'iseq' is what sequences to read (default is all)
  # value of 'ret' determines format of return value:
  #   count: amino acid composition (same columns as thermo$protein, can be used by add.protein)
  #        or nucleic acid base composition (A-C-G-T)
  #   seq: amino acid sequence
  #   fas: fasta entry
  # value of 'id' is used for 'protein' in output table,
  #   otherwise ID is parsed from FASTA header (can take a while)
  is.nix <- Sys.info()[[1]]=="Linux"
  if(is.nix & is.null(lines)) {
    message("read.fasta: reading ",basename(file))
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
    if(is.null(lines)) {
      lines <- readLines(file)
      message("read.fasta: reading ",basename(file))
    }
    nlines <- length(lines)
    if(is.null(ihead)) ihead <- which(substr(lines,1,1)==">")
    linefun <- function(i1,i2) lines[i1:i2]
  }
  # identify the lines that begin and end each sequence
  begin <- ihead + 1
  end <- ihead - 1
  end <- c(end[-1], nlines)
  # use all or selected sequences
  if(is.null(iseq)) iseq <- seq_along(begin)
  # just return the lines from the file
  if(ret=="fas") {
    iline <- numeric()
    for(i in iseq) iline <- c(iline,(begin[i]-1):end[i])
    return(lines[iline])
  }
  # get each sequence from the begin to end lines
  seqfun <- function(i) paste(linefun(begin[i],end[i]),collapse="")
  sequences <- lapply(iseq, seqfun)
  # organism name is from file name
  # (basename minus extension)
  bnf <- strsplit(basename(file),split=".",fixed=TRUE)[[1]][1]
  organism <- bnf
  # protein/gene name is from header line for entry
  # (strip the ">" and go to the first space)
  missid <- missing(id)
  if(is.null(id)) id <- as.character(palply("", iseq, function(j) {
    # get the text of the line
    f1 <- linefun(ihead[j],ihead[j])
    # stop if the first character is not ">"
    # or the first two charaters are "> "
    if(substr(f1,1,1)!=">" | length(grep("^> ",f1)>0))
      stop(paste("file",basename(file),"line",j,"doesn't begin with FASTA header '>'."))
    # discard the leading '>'
    f2 <- substr(f1, 2, nchar(f1))
    # keep everything before the first space
    return(strsplit(f2," ")[[1]][1])
  } ))
  if(ret=="count") {
    counts <- count.aa(sequences, start, stop, type)
    ref <- abbrv <- NA
    chains <- 1
    if(type=="protein") {
      colnames(counts) <- aminoacids(3)
      # 20090507 made stringsAsFactors FALSE
      out <- cbind(data.frame(protein=id, organism=organism,
        ref=ref, abbrv=abbrv, chains=chains, stringsAsFactors=FALSE), counts)
      # 20170117 extra processing for files from UniProt
      isUniProt <- grepl("\\|......\\|.*_", out$protein[1])
      if(isUniProt & missid) {
        p1 <- sapply(strsplit(out$protein, "\\|"), "[", 1)
        p2 <- sapply(strsplit(out$protein, "\\|"), "[", 2)
        p3 <- sapply(strsplit(out$protein, "\\|"), "[", 3)
        out$abbrv <- sapply(strsplit(p3, "_"), "[", 1)
        out$organism <- sapply(strsplit(p3, "_"), "[", 2)
        out$protein <- paste0(p1, "|", p2)
      }
      out
    } else if(type %in% c("DNA", "RNA")) {
      cbind(data.frame(gene=id, organism=organism,
        ref=ref, abbrv=abbrv, chains=chains, stringsAsFactors=FALSE), counts)
    }
  } else return(sequences)
}

uniprot.aa <- function(protein, start=NULL, stop=NULL) {
  # download protein sequence information from UniProt
  iprotein <- numeric()
  # construct the initial URL
  proteinURL <- paste("http://www.uniprot.org/uniprot/", protein, sep="")
  message("uniprot.aa: trying ", proteinURL, " ...", appendLF=FALSE)
  # try loading the URL, hiding any warnings
  oldopt <- options(warn=-1)
  URLstuff <- try(readLines(proteinURL),TRUE)
  options(oldopt)
  if(class(URLstuff)=="try-error") {
    message(" ::: FAILED :::")
    return(NA)
  }
  # 20091102: look for a link to a fasta file
  link <- grep("/uniprot/.*fasta", URLstuff)
  if(length(link) > 0) linkline <- URLstuff[[link[1]]]
  else {
    message(" ::: FAILED :::")
    return(NA)
  }
  # extract accession number from the link
  linkhead <- strsplit(linkline, ".fasta", fixed=TRUE)[[1]][1]
  accession.number <- tail(strsplit(linkhead, "/uniprot/", fixed=TRUE)[[1]], 1)
  message(" accession ", accession.number, " ...")
  # now download the fasta file
  fastaURL <- paste("http://www.uniprot.org/uniprot/", accession.number, ".fasta", sep="")
  URLstuff <- readLines(fastaURL)
  # get the header information / show  the user
  header <- URLstuff[[1]]
  header3 <- strsplit(header, "|", fixed=TRUE)[[1]][3]
  headerP_O <- strsplit(header3, " ")[[1]][1]
  header.id <- strsplit(header, headerP_O)[[1]][1]
  header.id <- substr(header.id, 2, nchar(header.id)-1)
  header.organism <- strsplit(headerP_O, "_")[[1]][2]
  message(paste0(header), appendLF=FALSE)
  # 20130206 use read.fasta with lines, start, stop arguments
  aa <- read.fasta(file="", lines=URLstuff, start=start, stop=stop)
  message(" (length ", sum(aa[1, 6:25]), ")", sep="")
  aa$protein <- header.id
  aa$organism <- header.organism
  return(aa)
}

count.aa <- function(seq, start=NULL, stop=NULL, type="protein") {
  # count amino acids or DNA bases in one or more sequences given as elements of the list seq
  # put them in alphabetical order (amino acids: same order as in thermo$protein)
  if(type=="protein") letts <- aminoacids(1)
  else if(type=="DNA") letts <- c("A", "C", "G", "T")
  else stop(paste("unknown sequence type", type))
  # to count the letters in each sequence
  countfun <- function(seq, start, stop) {
    # get a substring if one or both of start or stop are given
    # if only one of start or stop is given, get a default value for the other
    if(!is.null(start)) {
      if(is.null(stop)) stop <- nchar(seq)
      seq <- substr(seq, start, stop)
    } else if(!is.null(stop)) {
      seq <- substr(seq, 1, stop)
    }
    # the actual counting
    nnn <- table(strsplit(toupper(seq), "")[[1]])
    # get them in alphabetical order
    ilett <- match(names(nnn), letts)
    # in case any letters aren't in our alphabet
    ina <- is.na(ilett)
    if(any(ina)) message(paste("count.aa: unrecognized letter(s) in", type, "sequence:",
      paste(names(nnn)[ina], collapse=" ")))
    count <- numeric(length(letts))
    count[ilett[!ina]] <- nnn[!ina]
    return(count)
  }
  # counts for each sequence
  a <- palply("", seq, countfun, start, stop)
  a <- t(as.data.frame(a, optional=TRUE))
  # clean up row/column names
  colnames(a) <- letts
  rownames(a) <- 1:nrow(a)
  return(a)
}

