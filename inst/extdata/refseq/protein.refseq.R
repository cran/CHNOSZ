# protein.refseq.R
# calculate the average amino acid
# composition of proteins for each taxid
# 20100704 jmd

# we need CHNOSZ for read.fasta
require(CHNOSZ)

# where the microbial*.protein.faa.gz files are kept
faadir <- "protein"

# read gi.taxid if it already doesn't exist
if(!exists("gi.taxid")) {
  cat("reading gi.taxid.txt ... could take a while!\n")
  gi.taxid <- read.table("gi.taxid.txt")
  # we should make sure that all 3 columns are numeric (not factors)
}

protein.refseq <- function(n=NULL) {
  # the list of sequence files (microbial44.protein.fa.gz etc)
  files <- dir(faadir)
  # loop over each one, getting the contents
  if(is.null(n)) n <- length(files)
  out <- NULL
  for(i in 1:n) {
    cat(paste("(",i,"of",n,") "))
    # check if it's a fasta file
    if(length(grep(".protein.faa.gz$",files[i]))==0) {
      cat(paste("skipping \'",files[i],"\' ...\n"))
      next
    } else { 
      cat(paste("reading \'",files[i],"\' ... "))
    }
    # read the fasta file
    fa <- readLines(paste(faadir,files[i],sep="/"))
    # get the gi's for this file
    ihead.orig <- ihead <- grep("^>",fa)
    cat(paste("file has",length(ihead),"sequences\n"))
    gi <- def2gi(fa[ihead])
    # now find the taxids that belong to these gi's
    cat("finding matching taxids in catalog ...")
    gi.order <- order(as.character(gi))
    gi <- gi[gi.order]
    ihead <- ihead[gi.order]
    write.table(gi,"gi.txt",row.names=FALSE,col.names=FALSE,quote=FALSE)
    system("join -1 1 -2 1 gi.txt gi.taxid.txt > gi.taxid.match")
    gi.taxid.match <- read.table("gi.taxid.match",header=FALSE)
    tax <- gi.taxid.match$V2
    utax <- unique(tax)
    cat(paste(" found",length(utax),"unique taxa\n"))
    # now loop over each taxon
    for(j in 1:length(utax)) {
      jtax <- which(tax==utax[j])
      nseq.read <- length(jtax)
      cat(paste("taxon",utax[j],"has",nseq.read,"sequences"))
      # count the number of sequences in the catalog
      jjtax <- which(gi.taxid$V2==utax[j])
      nseq.cat <- length(jjtax)
      if(nseq.cat==nseq.read) cat(" (same as catalog)")
      else cat(paste(" WARNING: catlog has",nseq.cat))
      # get the amino acid compositions of the sequences
      aa <- read.fasta(file="",i=ihead[jtax],lines=fa,ihead=ihead.orig)
      naa.read <- sum(aa[,6:25])
      # count the number of amino acids in the catalog
      naa.cat <- sum(gi.taxid$V3[jjtax])
      if(naa.cat==naa.read) cat(paste("... all of",naa.cat,"amino acids "))
      else cat(paste("...",naa.read,"of",naa.cat,"amino acids "))
      # prepare the dataframe
      # average amino acid composition
      aa[1,6:25] <- colSums(aa[,6:25])/nrow(aa)
      aa <- aa[1,]
      aa$protein <- "refseq"
      aa$organism <- utax[j]
      aa$chains <- 1
      # the number of amino acids read from files
      aa$abbrv <- naa.read
      # what is the name of this organism
      itax1 <- match(utax[j],gi.taxid.match$V2)
      gi1 <- gi.taxid.match$V1[itax1]
      igi1 <- match(gi1,gi)
      iihead <- ihead[igi1]
      orgname <- s2c(fa[iihead],sep=" [",keep.sep=TRUE)[2]
      orgname <- substr(orgname,2,nchar(orgname))
      cat(paste(orgname,"\n"))
      # numbers of sequences, amino acids, organism name
      aa <- cbind(aa,data.frame(nseq=nseq.read,naa=naa.read,nseq.cat=nseq.cat,naa.cat=naa.cat,orgname=orgname))
      # keep track of source file name 
      # (append number of sequences according to the catalog)
      aa$ref <- paste(sub(".protein.faa.gz","",files[i]),"(",nseq.read,",",naa.read,")",sep="")
      if(j==1) aa.out <- aa else aa.out <- rbind(aa.out,aa)
    }
    if(is.null(out)) out <- aa.out else out <- rbind(out,aa.out)
  }
  # now we have to combine taxids that showed up more than once
  # e.g. microbial10;microbial9(2231743)
  duptax <- unique(out$organism[duplicated(out$organism)])
  if(length(duptax) > 0) {
    for(i in 1:length(duptax)) {
      id <- which(out$organism==duptax[i])
      # add them together, weighted by numbers of sequences
      aa <- colSums(out$nseq[id]*out[id,6:25])/sum(out$nseq[id])
      # keep the result in the first row
      out[id[1],6:25] <- aa
      # total number of sequences found in the files
      out$nseq[id[1]] <- sum(out$nseq[id])
      # total number of amino acids found in the sequence files
      out$naa[id[1]] <- sum(out$naa[id])
      out$abbrv[id[1]] <- out$naa[id[1]]
      # write down the file names 
      out$ref[id[1]] <- c2s(out$ref[id],sep=";")
      # drop the duplicated ones
      out <- out[-id[2:length(id)],]
    }
  }
  # summarize the organism name and numbers of
  # sequences/amino acids read vs. the catalog
  for(i in 1:nrow(out)) {
    dnseq <- out$nseq[i] - out$nseq.cat[i]
    dnaa <- out$naa[i] - out$naa.cat[i]
    if(dnseq != 0 | dnaa != 0) out$ref[i] <- paste("[",dnseq,",",dnaa,"]",out$ref[i],sep="")
    out$ref[i] <- paste(out$ref[i],out$orgname[i],sep="")
  }
  # we can drop the extra columns now
  out <- out[,1:25]
  # and do a little rounding
  out[,6:25] <- round(out[,6:25],3)
  write.csv(out,"protein_refseq.csv",row.names=FALSE)
  return(out)
}
