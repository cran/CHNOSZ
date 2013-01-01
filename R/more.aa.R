# CHNOSZ/more.aa.R
# get amino acid compositions of proteins from 
# model organisms really exciting!
# (Eco.csv or Sce.csv)

more.aa <- function(protein=NULL, organism) {
  # return the composition of one or more proteins from
  # a "model organism", E. coli (Eco) or S. cerevisiae (Sce)
  # extracted from get.protein 20120519
  datapath <- paste("extdata/protein/", organism, ".csv.xz", sep="")
  datafile <- system.file(datapath, package="CHNOSZ")
  if(datafile=="") stop(paste("missing", datapath))
  mydata <- read.csv(datafile)
  # if protein is not supplied, just give some information about the datafile
  if(is.null(protein)) {
    msgout("more.aa: ", datapath, " has data for ", nrow(mydata), " proteins\n")
    return(invisible())
  }
  if(organism=="Sce") {
    # which columns to search for matches
    searchcols <- c("OLN", "OLN")
    # which columns have the amino acids in the order of thermo$protein 
    iaa <- c(1,5,4,7,14,8,9,10,12,11,13,3,15,6,2,16,17,20,18,19) + 2
  } else if(organism=="Eco") {
    # which columns to search for matches
    searchcols <- c("protein", "abbrv")
    # which columns have the amino acids in the order of thermo$protein 
    iaa <- 1:20 + 5
  }
  # iterate over a list
  waslist <- TRUE
  out <- list()
  if(!is.list(protein)) {
    waslist <- FALSE
    protein <- list(protein)
  }
  for(i in 1:length(protein)) {
    # find the matches
    icols <- match(searchcols, colnames(mydata))
    imatch <- match(protein[[i]], mydata[, icols[1]])
    imatch2 <- match(protein[[i]], mydata[, icols[2]])
    # use not-NA matches for "abbrv" in Eco.csv
    imatch[!is.na(imatch2)] <- imatch2[!is.na(imatch2)]
    # report and remember the unsuccessful matches
    if(all(is.na(imatch))) stop("no proteins found!")
    inotmatch <- which(is.na(imatch)) 
    if(length(inotmatch) > 0) {
      if(length(inotmatch)==1) verb <- " was" else verb <- " were"
      msgout("more.aa: ", paste(protein[[i]][inotmatch], collapse=" "), verb, " not matched\n")
    }
    aa <- data.frame(mydata[imatch, iaa])
    # add the identifying columns
    organism <- rep(organism[[1]], length(protein[[i]]))
    abbrv <- rep(NA, length(protein[[i]]))
    ref <- rep(NA, length(protein[[i]]))
    chains <- rep(1, length(protein[[i]]))
    chains[inotmatch] <- NA
    precols <- data.frame(protein[[i]], organism, ref, abbrv, chains, stringsAsFactors=FALSE)
    colnames(precols)[1] <- "protein"
    colnames(aa) <- aminoacids(3)
    aa <- cbind(precols, aa)
    out <- c(out, list(aa))
  }
  # done!
  if(!waslist) return(out[[1]])
  else return(out)
}


