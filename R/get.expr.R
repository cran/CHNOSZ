# CHNOSZ/get.expr.R
# get abundance data from a protein expression experiment
# and add the proteins to the working instance of CHNOSZ
# 20101006 jmd

get.expr <- function(file,idcol,abundcol,seqfile,filter=NULL,is.log=FALSE,loga.total=0) {
  # extracted from findit.R 20100926 jmd
  # file: the name of the file with sequence ids and abundance data
  # idcol: the column of the data file that has the sequence ids
  # abundcol: the column of the data file that has the abundances
  # seqfile: the name of the file with the sequences
  # filter: optional column names/search terms to filter the results
  # is.log: are the abundances in logarithmic units?
  # loga.tot: logarithm of total activity

  # the name of the data file
  edata <- read.csv(file)
  # which columns for IDs and abundances
  if(is.numeric(idcol)) iid <- idcol else iid <- match(idcol,colnames(edata))
  if(is.numeric(abundcol)) ia <- abundcol else ia <- match(abundcol,colnames(edata))
  # first clean up the data file: duplicated sequences ids
  idup <- duplicated(edata[,iid])
  edata <- edata[!idup,]
  # remove NA sequence ids
  ina <- is.na(edata[,iid])
  edata <- edata[!ina,]
  # apply a filter if requested
  if(!is.null(filter)) {
    jfilter <- match(names(filter),colnames(edata))
    ifilter <- grep(filter[[1]],edata[,jfilter])
    edata <- edata[ifilter,]
  }
  # tell the user what we're looking for
  cat(paste("get.expr: searching for",nrow(edata),"entries... "))
  # read the fasta file
  if(seqfile %in% c("SGD","ECO","HUM")) {
    p <- get.protein(edata[,iid],seqfile)
    # in case we didn't get all the proteins
    j <- match(p$protein,edata[,iid])
    edata <- edata[j,]
  } else {
    flines <- readLines(seqfile)
    ihead <- grep("^>",flines)
    # we need fixed=TRUE so that accession numbers with dots in them
    # don't get matched to the wrong entries
    iihead <- as.numeric(mylapply(edata[,iid],function(x) grep(x,flines[ihead],fixed=TRUE)[1]))
    # if any were not matched, take them out
    if(any(is.na(iihead))) {
      ina <- which(is.na(iihead))
      iihead <- iihead[-ina]
      cat("got",length(iihead),"\n")
      cat(paste("get.expr: the following IDs were not found:",c2s(edata[,iid][ina]),"\n"))
      edata <- edata[-ina,]
    } else cat("got them all!\n")
    j <- ihead[iihead]
    p <- read.fasta(seqfile,j,lines=flines,ihead=ihead) 
  }
  # now store the experimental abundances

  loga.target <- edata[,ia]
  # take their logarithms if they're not already taken
  if(missing(is.log)) {
    # make a guess: if the column name has "log" don't take the logarithm
    if(length(grep("log",colnames(edata)[ia]))==0) loga.target <- log10(loga.target)
  } else if(!is.log) loga.target <- log10(loga.target)
  # scale the abundances so that total activity of residues is unity
  if(!is.null(loga.total)) {
    pl <- rowSums(p[,6:25])
    loga.target <- unitize(loga.target,pl,loga.total)
  }
  # now get the proteins ids used in CHNOSZ,
  # in the same order as the experimental list
  ip <- add.protein(p)
  # return our results
  out <- list(id=edata[,iid],iprotein=ip,loga.ref=loga.target)
  return(out)
}
