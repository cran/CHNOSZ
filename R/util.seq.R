# CHNOSZ/util.seq.R
# Copyright (C) 2006-2009 Jeffrey M. Dick
# functions to work with sequence data

aminoacids <- function(seq=NULL,nchar=1) {
  # if no sequence is given just return the
  # abbreviations of the amino acids
  # these are in the same order as in thermo$protein
  aa <- c('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y')
  aa.3 <- c('Ala','Cys','Asp','Glu','Phe','Gly','His','Ile','Lys','Leu',
            'Met','Asn','Pro','Gln','Arg','Ser','Thr','Val','Trp','Tyr')
  aa.NA <- c("alanine","cysteine","aspartic acid","glutamic acid","phenylalanine","glycine",
    "histidine","isoleucine","lysine","leucine","methionine","asparagine","proline","glutamine",
    "arginine","serine","threonine","valine","tryptophan","tyrosine")
  aa.Z <- c("alanine","cysteine","aspartate","glutamate","phenylalanine","glycine",
    "histidine","isoleucine","lysinium","leucine","methionine","asparagine","proline","glutamine",
    "argininium","serine","threonine","valine","tryptophan","tyrosine")
  if(is.null(seq)) {
    if(is.na(nchar)) return(aa.NA)
    if(nchar=="Z") return(aa.Z)
    if(as.numeric(nchar)==1) return(aa)
    if(as.numeric(nchar)==3) return(aa.3)
  }
  # sequences are given as elements of the list seq
  # to count the number of each amino acids in a sequence
  count.aa <- function(aa,seq) sum(seq==aa)
  count.i <- function(i,seq) as.numeric(lapply(aa,count.aa,strsplit(toupper(seq[i]),"")[[1]]))
  # count amino acids in each sequence
  a <- t(as.data.frame(mylapply(1:length(seq),count.i,seq),optional=TRUE))
  # clean up row/column names
  colnames(a) <- aa
  rownames(a) <- 1:nrow(a)
  return(a)
}

nucleicacids <- function(seq=NULL,type="DNA",comp=NULL,comp2=NULL) {
  # count bases or compute the formula, e.g.
  # n <- nucleicacids(list("AGCT","TTTT"))  # a dataframe of counts
  # f <- nucleicacids(n)  # a series of formulas
  # 20090926 jmd
  if(is.null(seq)) stop("please provide a sequence")
  if(type=="DNA") {
    na <- c("A","C","G","T")
    na.NA <- c("adenine","cytosine","guanine","thymine")
  } else if(type=="RNA") {
    na <- c("U","G","C","A")
    na.NA <- c("uracil","guanine","cytosine","adenine")
  } else stop(paste("invalid type:",type))
  if(is.data.frame(seq)) {
    # compute the chemical formula of bases
    if(!all(na %in% colnames(seq))) {
      nabases <- c2s(na[which(!na %in% colnames(seq))],sep=" ")
      stop(paste("requested type is",type,"but",nabases,"is/are not in the colnames of supplied dataframe"))
    }
    f.base <- thermo$obigt$formula[info(na.NA)]
    f0 <- makeup("C0H0N0O0")
    for(i in 1:4) {
      m.base <- makeup(f0,makeup(f.base[i]))
      if(i==1) m <- m.base else m <- cbind(m,m.base)
    }
    m <- as.matrix(m)
    f.out <- character()
    for(i in 1:nrow(seq)) {
      s <- t(seq[i,])
      c <- m %*% s
      c <- as.data.frame(c,row.names=row.names(c))
      colnames(c) <- "count"
      f <- makeup(makeup(c,""),"")
      f.out <- c(f.out,f)
    }
    return(f.out)
  } else {
    # count the numbers of nucleic acid bases in a sequence
    # sequences are given as elements of the list seq
    # to count the number of each amino acids in a sequence
    count.na <- function(na,seq) sum(seq==na)
    count.i <- function(i,seq) as.numeric(lapply(na,count.na,strsplit(toupper(seq[i]),"")[[1]]))
    # count bases in each sequence
    n <- t(as.data.frame(mylapply(1:length(seq),count.i,seq),optional=TRUE))
    n <- as.data.frame(n)
    # clean up row/column names
    colnames(n) <- na
    rownames(n) <- 1:nrow(n)
    # return the complement if requested e.g.
    # nucleicacids(x,type,"DNA")  # DNA complement
    # nucleicacids(x,type,"RNA")  # RNA complement
    # nucleicacids(x,type,"DNA","RNA")  # DNA, then RNA complement
    if(!is.null(comp)) {
      if(comp=="DNA") colnames(n) <- c("T","G","C","A")
      else if(comp=="RNA") colnames(n) <- c("U","G","C","A")
      else stop(paste("invalid complement request:",comp))
    }
    if(!is.null(comp2)) {
      if(comp2=="DNA") colnames(n) <- c("A","C","G","T")
      else if(comp2=="RNA") colnames(n) <- c("A","C","G","U")
      else stop(paste("invalid complement request:",comp))
    }
    return(n)
  }
}


