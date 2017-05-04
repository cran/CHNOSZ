# CHNOSZ/read.expr.R
# protein abundance and localization data from experiments
# yeastgfp: protein localization and abundance from yeastgfp.csv
# read.expr: protein abundance from other files (including in extdata/abundance)

yeastgfp <- function(location=NULL, exclusive=TRUE) {
  # return a list of ORFs and protein abundances for a subcellular location
  # using data from the YeastGFP project 
  # (yeastgfp.csv data file added to CHNOSZ_0.8, 20090422)
  ypath <- "extdata/abundance/yeastgfp.csv.xz"
  yfile <- system.file(ypath, package="CHNOSZ")
  # yeastgfp preprocessing
  ygfp <- read.csv(yfile)
  # convert factors to numeric w/o NA coercion warnings
  ygfp$abundance <- as.numeric.nowarn(as.character(ygfp$abundance))
  # if location is NULL, just report on the content of the file
  # and return the names of the locations
  if(is.null(location)) {
    message("yeastgfp: ", ypath, " has ", nrow(ygfp), " localizations and ",
      length(ygfp$abundance[!is.na(ygfp$abundance)]), " abundances")
    return(invisible(colnames(ygfp)[6:28]))
  }
  # iterate over multiple locations
  out <- list()
  for(i in 1:length(location)) {
    # what location do we want?
    ncol <- match(location[i], colnames(ygfp)[6:28]) + 5
    if(is.na(ncol)) ncol <- agrep(location[i], colnames(ygfp)[6:28])[1] + 5
    if(is.na(ncol)) stop(paste(location[i], "is not one of the subcellular locations in", ypath))
    thisygfp <- ygfp[, ncol]
    if(exclusive) {
      # find the number of localizations of each ORF
      localizations <- numeric(nrow(ygfp))
      for(j in 6:28) localizations <- localizations + as.logical(ygfp[,j])
      if(all(localizations[thisygfp] > 1)) message("yeastgfp: no exclusive localization found for ",location[i],
        " ... using non-exclusive localizations",sep="")
      else thisygfp <- thisygfp & ! localizations > 1
    }
    protein <- as.character(ygfp$yORF[thisygfp])
    abundance <- ygfp$abundance[thisygfp]
    if(length(location)==1) out <- list(protein=protein, abundance=abundance)
    else {
      out$protein <- c(out$protein, list(protein))
      out$abundance <- c(out$abundance, list(abundance))
    }
  }
  return(out)
}

read.expr <- function(file, idcol, abundcol, filter=NULL) {
  ## read protein expression data from files
  ## extracted from findit.R 20100926 jmd
  ## file: the name of the file with sequence ids and abundance data
  ## idcol: the column of the data file that has the sequence ids
  ## abundcol: the column of the data file that has the abundances
  ## filter: optional column names/search terms to filter the results
  # the name of the data file
  edata <- read.csv(file, stringsAsFactors=FALSE, check.names=FALSE)
  # which columns for IDs and abundances
  if(is.numeric(idcol)) iid <- idcol else iid <- match(idcol, colnames(edata))
  if(is.na(iid)) stop("unidentified protein ID column in", file)
  if(is.numeric(abundcol)) ia <- abundcol else ia <- match(abundcol, colnames(edata))
  if(is.na(ia)) stop("unidentified protein abundance column in", file)
  # first clean up the data file: duplicated sequences ids
  idup <- duplicated(edata[, iid])
  edata <- edata[!idup, ]
  # remove NA sequence ids
  ina <- is.na(edata[, iid])
  edata <- edata[!ina, ]
  # apply a filter if requested
  if(!is.null(filter)) {
    ifilter <- 1:nrow(edata)
    for(i in 1:length(filter)) ifilter <- intersect(ifilter, grep(filter[[i]], edata[, names(filter)[[i]]]))
    edata <- edata[unique(ifilter), ]
  }
  # that should be it
  protein <- edata[, iid]
  abundance <- edata[, ia]
  return(list(protein=protein, abundance=abundance))
#  # take their logarithms if they're not already taken
#  if(missing(is.log)) {
#    # make a guess: if the column name has "log" don't take the logarithm
#    if(length(grep("log",colnames(edata)[ia]))==0) loga.target <- log10(loga.target)
#  } else if(!is.log) loga.target <- log10(loga.target)
#  # scale the abundances so that total activity of residues is unity
#  if(!is.null(loga.total)) {
#    pl <- rowSums(p[,6:25])
#    loga.target <- unitize(loga.target,pl,loga.total)
#  }
}
