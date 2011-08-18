# CHNOSZ/util.stat.R
# miscellaneous statistical functions

qqr <- function(a) {
  # to calculate the correlation coefficient on a q-q plot
  qqn <- qqnorm(a,plot.it=FALSE)
  return(cor(qqn[[1]],qqn[[2]]))
}

rmsd <- function(a,b) {
  # to calculate the root mean square deviation
  d <- b - a
  sd <- d^2
  msd <- mean(sd)
  rmsd <- sqrt(msd)
  return(rmsd)
}

cvrmsd <- function(a,b) {
  # to calculate the coefficient of variation of the RMSD
  rmsd <- rmsd(a,b)
  # the observed values are in a
  cvrmsd <- rmsd/mean(a)
  return(cvrmsd)
}

spearman <- function(a,b) {
  # calculate Spearman's rho (rank correlation coefficient)
  # based on help(dSpearman) in package SuppDists
  if(length(a)!=length(b)) stop("a and b must be same length")
  if(any(is.na(a)) | any(is.na(b))) return(NA)
  ra <- rank(a)
  rb <- rank(b)
  dr <- rb - ra
  d <- sum(dr^2)
  r <- length(a)
  x <- 1-6*d/(r*(r^2-1))
  return(x)
}

lograt <- function(a,b) {
  # calculate the difference between a and b
  # (i.e. log10 of the ratio of activities b/a)
  # a - list of single values
  # b - list of values, any dimension
  mydim <- dim(b[[1]])
  out <- b
  for(i in 1:length(b)) {
    out[[i]] <- as.vector(b[[i]]) - a[[i]]
    dim(out[[i]]) <- mydim
  }
  return(out)
}
