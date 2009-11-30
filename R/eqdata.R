# eqdata.R
# extract data from eq3/6 output files
# 20091028 jmd

eqdata <- function(file,species,prop="log act") {

  # string constants that show up in the file
  zistring <- " stepping to zi="
  aqstring <- "   species                moles        grams"
  Tstring <- "                     temperature    ="

  # to read the file
  eqlines <- function(file) readLines(file)

  # to find the lines identifying each step of xi
  zilines <- function(lines) grep(zistring,lines)

  # to find the lines indicating the temperature
  Tlines <- function(lines) grep(Tstring,lines)

  # to find the lines heading the moles etc of aqueous species
  aqlines <- function(lines) grep(aqstring,lines)

  # to read the moles etc of aqueous species
  getaqprop <- function(lines,iaq,prop,species) {
    # the names of the properties
    props <- c("species","moles","grams","conc","log conc","log g","log act")
    iprop <- match(prop,props)
    if(is.na(iprop)) stop(paste('property "',prop,'" is not available',sep=''))
    # the property of a species is NA unless we find it
    # so we make a list with an entry for each species
    # the entries are NAs repeated to length of iaq
    prop.out <- lapply(1:length(species),function(x) rep(NA,length(iaq)))
    for(i in 1:length(iaq)) {
      # the line number of the header
      n <- iaq[i]
      # start reading two lines below the header
      dn <- 2
      myline <- lines[n+dn]
      # read until we reach the end of the block
      while(!identical(myline,"")) {
        myvalues <- strsplit(myline," ")[[1]]
        myvalues <- myvalues[!myvalues==""]
        # did we hit one of the desired species?
        ispecies <- match(myvalues[1],species)
        if(!is.na(ispecies)) prop.out[[ispecies]][i] <- as.numeric(myvalues[iprop])
        dn <- dn + 1
        myline <- lines[n+dn]
      }
    }
    prop.out <- as.data.frame(prop.out)
    colnames(prop.out) <- species
    return(prop.out)
  }

  # to get the values of zi for each of the aqueous species blocks
  getzi <- function(lines,izi,iaq) {
    # the zi linenumbers before each of the aqueous species blocks
    myzi <- sapply(iaq,function(x,y) {max(izi[x > y])},izi)
    # to get values of zi:
    # first split at the zistring (beginning of line)
    # then at the comma (after value)
    # then convert to numeric
    sapply(myzi,function(x) {as.numeric(strsplit(strsplit(lines[x],zistring)[[1]][2],",")[[1]][1])} )
  }

  # to get the values of T for each of the aqueous species blocks
  getT <- function(lines,iT,iaq) {
    # the T linenumbers before each of the aqueous species blocks
    myT <- sapply(iaq,function(x,y) {max(iT[x > y])},iT)
    # to get values of T:
    # first split at the Tstring (beginning of line)
    # then at " degrees" (after value)
    # then convert to numeric
    sapply(myT,function(x) {as.numeric(strsplit(strsplit(lines[x],Tstring)[[1]][2]," degrees")[[1]][1])} )
  }

  # put it all together
  # read each line of the file
  lines <- eqlines(file)
  # get line numbers for zi
  izi <- zilines(lines)
  # get line numbers for T
  iT <- Tlines(lines)
  # get line numbers for headers of aqueous species blocks
  iaq <- aqlines(lines)
  # get values of zi
  zi <- getzi(lines,izi,iaq)
  # get values of T
  T <- getT(lines,iT,iaq)
  # get properties of aqueous species
  props <- getaqprop(lines,iaq,prop,species)
  # make a single table
  out <- cbind(zi,T,props)
  # done!
  write.csv(out,paste(file,prop,"csv",sep="."),row.names=FALSE)
  return(invisible(out))

}
