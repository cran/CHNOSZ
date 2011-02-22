# CHNOSZ/blast.R
# functions to analyze BLAST output files
# 20100320 jmd

## process a blast tabular output file, counting 
## representation of each phylum
count.taxa <- function(file,gi.taxid,taxid.phylum,
  similarity=30,evalue=1e-5,max.hits=10,min.query=0,min.phylum=0,min.taxon=0) {
  # read the blast tabular file
  cat(paste("count.taxa: reading",file,"\n"))
  blast <- read.csv(file,header=FALSE,sep="\t",stringsAsFactors=FALSE)
  cat(paste("  read",nrow(blast),"lines with",length(unique(blast$V1)),"query sequences\n"))
  # take out rows that don't meet our desired similarity
  is <- which(blast$V3 >= similarity)
  blast <- blast[is,]
  cat(paste("  similarity filtering leaves",length(is),"lines and",length(unique(blast$V1)),"query sequences\n"))
  # take out rows that don't meet our desired e-value
  ie <- which(blast$V11 <= evalue)
  blast <- blast[ie,]
  cat(paste("  evalue filtering leaves",length(ie),"lines and",length(unique(blast$V1)),"query sequences\n"))
  # now take only max hits for each query sequence
  query.shift <- query <- blast$V1
  lq <- length(query)
  # for short (i.e., 1 query sequence) files, make sure that
  # the hits get counted
  query.shift[max((lq-max.hits+1),1):lq] <- -1
  #query.shift <- query.shift[c((max.hits+1):lq,1:max.hits)]
  query.shift <- query.shift[c((lq-max.hits+1):lq,1:(lq-max.hits))]
  ib <- which(query!=query.shift)
  blast <- blast[ib,]
  cat(paste("  max hits filtering leaves",length(ib),"lines and",length(unique(blast$V1)),"query sequences\n"))
  # what are gi numbers of the hits
  gi <- blast$V2
  ugi <- unique(gi)
  # what taxid do they hit
  cat("  getting taxids ... ")
  # we use def2gi to extract just the gi numbers
  def2gi <- function(def) {
    # extract gi numbers from FASTA deflines 20110131
    stuff <- strsplit(def,"\\|")
    gi <- sapply(1:length(stuff),function(x) {
      # the gi number should be in the 2nd position (after "gi")
      if(length(stuff[[x]])==1) return(stuff[[x]][1])
      else return(stuff[[x]][2])
    })
    return(gi)
  }
  imatch <- match(def2gi(ugi),def2gi(gi.taxid[[1]]))
  utaxid <- gi.taxid[[2]][imatch]
  # what phyla are these
  cat("getting taxon names ... ")
  iphy <- match(utaxid,taxid.phylum$taxid)
  uphyla <- taxid.phylum$phylum[iphy]
  uspecies <- taxid.phylum$species[iphy]
  cat(paste(length(unique(uphyla)),"phyla,",length(unique(utaxid)),"taxa\n"))
  # now expand phyla into our blast table
  igi <- match(gi,ugi)
  taxid <- utaxid[igi]
  phylum <- uphyla[igi]
  species <- uspecies[igi]
  tax.out <- data.frame(taxid=taxid,phylum=as.character(phylum),species=species)
  blast.out <- blast[,c(1,2,3,11)]
  colnames(blast.out) <- c("query","subject","similarity","evalue")
  blast <- cbind(blast.out,tax.out)
  # drop taxa that do not appear a certain number of times
  blast$taxid <- as.character(blast$taxid)
  nt <- tapply(blast$taxid,blast$taxid,length)
  it <- which(nt/sum(nt) >= min.taxon)
  itt <- which(blast$taxid %in% names(nt)[it])
  blast <- blast[itt,]
  cat(paste("  min taxon abundance filtering leaves",length(unique(blast$query)),
    "query sequences,",length(unique(blast$phylum)),"phyla,",length(unique(blast$taxid)),"taxa\n"))
  # only take phylum assignments that make up at least a certain 
  # fraction ('amin') of hits to the query sequence
  uquery <- unique(blast$query)
  iquery <- match(uquery,blast$query)
  # function to select the (highest) represented phylum for each query
  iqfun <- function(i) {
    if((i-1)%%1000==0) cat(paste(i,""))
    myiq <- which(blast$query==uquery[i])
    # we don't have to do the calculation if there's only one hit
    if(length(myiq)==1) {
      iq <- iquery[i]
    } else {
      # take those hits, count each phyla, take the highest abundance,
      # check if it's above the minimum proportion of hits
      p <- as.character(blast$phylum[myiq])
      np <- tapply(p,p,length)
      pp <- np/sum(np)
      ip <- which.max(pp)
      if(pp[ip] < min.query) return(numeric())
      # use the first hit to that phylum
      myphy <- blast$phylum[myiq]
      iphy <- which(names(pp[ip])==myphy)
      iq <- myiq[iphy[1]]
    }
    return(iq)
  }
  # using lapply is tempting, but haven't got it working
  #iq <- as.numeric(lapply(1:length(uquery),iqfun))
  #iq <- iq[!is.na(iq)]
  if(min.query!=0) {
    cat("  filtering by phylum representation per query... ")
    iq <- numeric()
    for(i in 1:length(uquery)) iq <- c(iq,iqfun(i))
    cat("done!\n")
    blast <- blast[iq,]
    cat(paste("  min query representation filtering leaves",length(iq),"query sequences,",length(unique(blast$phylum)),
      "phyla,",length(unique(blast$taxid)),"taxa\n"))
  } else {
    # we'll just take the first hit for each query sequence
    blast <- blast[iquery,]
  }
  # now on to drop those phyla that are below a certain relative abundance
  blast$phylum <- as.character(blast$phylum)
  np <- tapply(blast$phylum,blast$phylum,length)
  ip <- which(np/sum(np) >= min.phylum)
  ipp <- which(blast$phylum %in% names(np)[ip])
  blast <- blast[ipp,]
  cat(paste("  min phylum abundance filtering leaves",length(ipp),"query sequences,",length(unique(blast$phylum)),
    "phyla,",length(unique(blast$taxid)),"taxa\n"))
  return(blast)
} 

