# CHNOSZ/util.seq.R
# functions to work with sequence data

aminoacids <- function(nchar=1, which=NULL) {
  # return the abbreviations or names of the amino acids
  # the following are all in the same order as thermo$protein
  # the single-letter codes
  aa1 <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", 
           "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
  # the 3-letter codes
  aa3 <- c("Ala", "Cys", "Asp", "Glu", "Phe", "Gly", "His", "Ile", "Lys", "Leu",
            "Met", "Asn", "Pro", "Gln", "Arg", "Ser", "Thr", "Val", "Trp", "Tyr")
  # names of the neutral amino acids
  aaneutral <- c("alanine", "cysteine", "aspartic acid", "glutamic acid", "phenylalanine", 
    "glycine", "histidine", "isoleucine", "lysine", "leucine", "methionine", "asparagine", 
    "proline", "glutamine", "arginine", "serine", "threonine", "valine", "tryptophan", "tyrosine")
  # names of the amino acids with ionized side chains
  aacharged <- c("alanine", "cysteinate", "aspartate", "glutamate", "phenylalanine", 
    "glycine", "histidinium", "isoleucine", "lysinium", "leucine", "methionine", "asparagine", 
    "proline", "glutamine", "argininium", "serine", "threonine", "valine", "tryptophan", "tyrosinate")
  # defaults are in the same order as in thermo$protein
  if(is.null(which)) which <- aa1
  # figure out which amino acids are wanted (can use 1- or 3-letter codes, or neutral names)
  if(all(nchar(which)==1)) iaa <- match(which, aa1)
  else if(all(nchar(which)==3)) iaa <- match(which, aa3)
  else iaa <- match(which, aaneutral)
  # return the desired abbreviations or names
  if(nchar==1) return(aa1[iaa])
  else if(nchar==3) return(aa3[iaa])
  else if(nchar=="") return(aaneutral[iaa])
  else if(nchar=="Z") return(aacharged[iaa])
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
    # the formulas of each of the bases
    f.base <- get("thermo")$obigt$formula[info(na.NA[match(colnames(seq),na)])]
    # loop over the base counts
    f.out <- character()
    for(i in 1:nrow(seq)) {
      # use makeup() with multipliers and sum=TRUE  20120119 jmd
      f <- as.chemical.formula(makeup(f.base, multiplier=as.numeric(seq[i,]), sum=TRUE))
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
    n <- t(as.data.frame(palply(1:length(seq),count.i,seq),optional=TRUE))
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


