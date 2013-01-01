# CHNOSZ/iprotein.R
# calculate properties of proteins 20061109 jmd
# reorganize protein functions 20120513

# iprotein - find rownumber in thermo$protein
# ip2aa - select amino acid counts (data frame) from thermo$protein
# aa2eos - perform group additivity calculations
# seq2aa - calculate amino acid counts from a sequence
# dl.aa - get amino acid counts from SWISS-PROT
# aasum - combine amino acid counts (sum, average, or weighted sum by abundance)
# read.aa - read amino acid counts from a file
# add.protein - add amino acid counts to thermo$protein (returns iprotein)

iprotein <- function(protein, organism=NULL) {
  # find the rownumber(s) of thermo$protein that matches
  # 'protein' numeric (the rownumber itself)
  # 'protein' character, e.g. LYSC_CHICK
  # 'protein' and 'organism', e.g. 'LYSC', 'CHICK'
  if(is.numeric(protein)) {
    iproteins <- 1:nrow(thermo$protein)
    protein[!protein %in% iproteins] <- NA
    iprotein <- protein
  } else {
    # from here we'll search by protein/organism pairs
    tp.po <- paste(thermo$protein$protein, thermo$protein$organism, sep="_")
    if(is.null(organism)) my.po <- protein
    else my.po <- paste(protein, organism, sep="_")
    iprotein <- match(my.po, tp.po)
  }
  # tell the user about NA's
  if(any(is.na(iprotein))) {
    nNA <- sum(is.na(iprotein))
    if(nNA==1) ptext <- "" else ptext <- "s"
    msgout("iprotein: ", sum(is.na(iprotein)), " protein", ptext, " not matched\n")
  }
  return(iprotein)
}

ip2aa <- function(protein, organism=NULL, residue=FALSE) {
  # return amino acid counts (rows from thermo$protein)
  # or 'protein' if it is a data frame
  if(is.data.frame(protein)) return(protein)
  iprotein <- iprotein(protein, organism)
  # drop NA matches
  iprotein <- iprotein[!is.na(iprotein)]
  out <- thermo$protein[iprotein, ]
  # compute per-residue counts
  if(residue) out[, 5:25] <- out[, 5:25]/rowSums(out[, 6:25])
  return(out)
}

aa2eos <- function(aa, state=thermo$opt$state) {
  # display and return the properties of
  # proteins calculated from amino acid composition
  # the names of the protein backbone groups depend on the state
  # [UPBB] for aq or [PBB] for cr
  if(state=="aq") bbgroup <- "UPBB" else bbgroup <- "PBB"
  # names of the AABB, sidechain and protein backbone groups
  groups <- c("AABB", colnames(aa)[6:25], bbgroup)
  # put brackets around the group names
  groups <- paste("[", groups, "]", sep="")
  # the rownumbers of the groups in thermo$obigt
  groups_state <- paste(groups, state)
  obigt_state <- paste(thermo$obigt$name, thermo$obigt$state)
  igroup <- match(groups_state, obigt_state)
  # the properties are in columns 8-20 of thermo$obigt
  groupprops <- thermo$obigt[igroup, 8:20]
  # the elements in each of the groups
  groupelements <- i2A(igroup)
  # a function to work on a single row of aa
  eosfun <- function(aa) {
    # numbers of groups: chains [=AABB], sidechains, protein backbone
    nchains <- as.numeric(aa[, 5])
    length <- sum(as.numeric(aa[, 6:25]))
    npbb <- length - nchains
    ngroups <- c(nchains, as.numeric(aa[, 6:25]), npbb)
    # the actual adding and multiplying of thermodynamic properties
    # hmm. seems like we have to split up the multiplication/transposition
    # operations to get the result into multiple columns. 20071213
    eos <- t(data.frame(colSums(groupprops * ngroups)))
    # to get the formula, add up and round the group compositions 20090331
    f.in <- round(colSums(groupelements * ngroups), 3)
    # take out any elements that don't appear (sometimes S)
    f.in <- f.in[f.in!=0]
    # turn it into a formula
    f <- as.chemical.formula(f.in)
    # now the species name
    name <- paste(aa$protein, aa$organism, sep="_")
    # make some noise for the user
    msgout("aa2eos: found ")
    msgout(name, " (", f, ", ")
    msgout(round(length, 3), " residues)\n")
    ref <- aa$ref
    header <- data.frame(name=name, abbrv=NA, formula=f, state=state, ref1=ref, ref2=NA, date=NA, stringsAsFactors=FALSE)
    eosout <- cbind(header, eos)
    return(eosout)
  }
  # loop over each row of aa
  out <- lapply(1:nrow(aa), function(i) eosfun(aa[i, ]))
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  return(out)
}

seq2aa <- function(protein, sequence) {
  # make a data frame from counting the amino acids in the sequence
  caa <- count.aa(sequence)
  colnames(caa) <- aminoacids(3)
  # a protein with no amino acids is sort of boring
  if(all(caa==0)) stop("no characters match an amino acid")
  ip <- suppressMessages(iprotein(protein))
  # now make the data frame
  po <- strsplit(protein, "_")[[1]]
  aa <- data.frame(protein=po[1], organism=po[2], ref=NA, abbrv=NA)
  aa <- cbind(aa, chains=1, caa)
  return(aa)
}

dl.aa <- function(protein) {
  # download protein sequence information from SWISS-PROT
  iprotein <- numeric()
  # construct the initial URL
  proteinURL <- paste("http://www.uniprot.org/uniprot/", protein, sep="")
  msgout("dl.aa: trying ", proteinURL, " ...")
  # try loading the URL, hiding any warnings
  oldopt <- options(warn=-1)
  URLstuff <- try(readLines(proteinURL),TRUE)
  options(oldopt)
  if(class(URLstuff)=="try-error") {
    msgout(" failed\n")
    return(NA)
  }
  # 20091102: look for a link to a fasta file
  linkline <- URLstuff[[grep("/uniprot/.*fasta", URLstuff)[1]]]
  # extract accession number from the link
  linkhead <- strsplit(linkline, ".fasta", fixed=TRUE)[[1]][1]
  accession.number <- tail(strsplit(linkhead, "/uniprot/", fixed=TRUE)[[1]], 1)
  msgout(" accession ", accession.number, " ...\n")
  # now download the fasta file
  fastaURL <- paste("http://www.uniprot.org/uniprot/", accession.number, ".fasta", sep="")
  URLstuff <- readLines(fastaURL)
  # show the name of the protein to the user
  header <- URLstuff[[1]]
  header2 <- strsplit(header, paste(protein, ""))[[1]][2]
  header3 <- strsplit(header2, " OS=")[[1]]
  protein.name <- header3[1]
  header4 <- strsplit(header3[2], " GN=")[[1]][1]
  header5 <- strsplit(header4[1], " PE=")[[1]]
  organism.name <- header5[1]
  msgout("dl.aa: ", protein.name, " from ", organism.name)
  # get rid of the header before counting amino acid letters
  URLstuff[[1]] <- ""
  aa <- count.aa(c2s(URLstuff, sep=""))
  msgout(" (length ", sum(aa[1,]), ")\n", sep="")
  colnames(aa) <- colnames(thermo$protein)[6:25]
  po <- strsplit(protein, "_")[[1]]
  out <- data.frame(protein=po[1], organism=po[2], ref=NA, abbrv=NA, chains=1, aa)
  return(out)
}

aasum <- function(aa, abundance=1, average=FALSE, protein=NULL, organism=NULL) {
  # returns the sum of the amino acid counts in aa,
  # multiplied by the abundances of the proteins
  abundance <- rep(abundance, length.out=nrow(aa))
  # drop any NA rows or abundances
  ina.aa <- is.na(aa$chains)
  ina.ab <- is.na(abundance)
  ina <- ina.aa | ina.ab
  if(any(ina)) {
    aa <- aa[!ina, ]
    abundance <- abundance[!ina]
    msgout("aasum: dropped ", sum(ina), " proteins with NA composition and/or abundance\n")
  }
  # we don't know how to deal with different numbers of polypeptide chains
  if(!all(aa$chains==aa$chains[1])) stop("different numbers of polypeptide chains")
  # multiply
  aa[, 6:25] <- aa[, 6:25] * abundance
  # sum
  out <- aa[1, ]
  out[, 5:25] <- colSums(aa[, 5:25])
  # average if told to do so
  if(average) {
    # polypeptide chains by number of proteins, residues by frequence
    out[, 5] <- out[, 5]/nrow(aa)
    out[, 6:25] <- out[, 6:25]/sum(abundance)
  }
  # add protein and organism names if given
  if(!is.null(protein)) out$protein <- protein
  if(!is.null(organism)) out$organism <- organism
  return(out)
}

read.aa <- function(file="protein.csv") {
  # if its a fasta file, read the sequences
  if(is.fasta(file)) {
    aa <- read.fasta(file)
    msgout("read.aa: first line in FASTA file is\n")
    msgout(readLines(file, n=1), "\n")
  } else {
    # 20090428 added colClasses here
    aa <- read.csv(file,colClasses=c(rep("character",4),rep("numeric",21)))
  }
  if(!identical(colnames(aa), colnames(thermo$protein))) 
    stop("format of", file, "is incompatible with thermo$protein")
  return(aa)
}

add.protein <- function(aa, print.existing=FALSE) {
  # add a properly constructed data frame of 
  # amino acid counts to thermo$protein
  if(!identical(colnames(aa), colnames(thermo$protein))) 
    stop("the value of 'aa' is not a data frame with the same columns as thermo$protein")
  # find any protein IDs that are duplicated
  po <- paste(aa$protein, aa$organism, sep="_")
  ip <- suppressMessages(iprotein(po))
  ipdup <- !is.na(ip)
  # now we're ready to go
  thermo$protein <<- rbind(thermo$protein, aa[!ipdup, ])
  rownames(thermo$protein) <<- NULL
  # return the new rownumbers
  ip <- iprotein(po)
  # make some noise
  msgout("add.protein: added ", nrow(aa)-sum(ipdup), " of ", nrow(aa), " proteins\n")
  if(!all(is.na(ipdup)) & print.existing) {
    potext <- paste(aa$protein[ipdup], aa$organism[ipdup], sep="_", collapse=" ")
    msgout("add.protein: skipped existing ", potext, "\n")
  }
  return(ip)
}
