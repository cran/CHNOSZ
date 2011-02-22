\name{util.blast}
\alias{util.blast}
\alias{count.taxa}
\title{Functions to Work with BLAST Output Files}

\description{
  Summarize the representation of taxa in a BLAST tabular output file, after filtering for similarity, E-value, and number of hits.
}

\usage{
  count.taxa(file, gi.taxid, taxid.phylum, similarity=30, evalue=1e-5,
    max.hits=10, min.query=0, min.phylum=0, min.taxon=0)
}

\arguments{
  \item{file}{character, name of BLAST tabular output file.}
  \item{gi.taxid}{list, first component is sequence identifiers (gi numbers), second is taxon ids.}
  \item{taxid.phylum}{dataframe, with columns taxid (taxon id), phylum (name of phylum), species (name of species).}
  \item{similarity}{numeric, hits above this similarity score are kept.}
  \item{evalue}{character, hits below this evalue are kept.}
  \item{max.hits}{numeric, up to this many hits are processed for each query sequence.}
  \item{min.query}{numeric, query sequence is counted if a single phylum makes up this fraction of hits.}
  \item{min.phylum}{numeric, phylum total abundance is kept if it is above this value.}
  \item{min.taxon}{numeric, this taxon is kept if it makes up at least this fraction of total.}
}

\details{
  The BLAST (Altschul et al., 1997) tabular output file can be generated using the -m 8 switch to the \dQuote{blastall} command. The BLAST results are tied to taxon ids using the \code{gi.taxid} list. For this to work, the sequence identifiers in BLAST database (and that show up in the output file) must be the same as those in the \code{gi.taxid} list. 

}

\value{
  A dataframe with one row for each query sequence that remains after the filtering operations. The columns are query, subject (i.e., sequence id or gi number), similarity, evalue, taxid, phylum and species.
}

\examples{
   # process some of the BLAST output for proteins
   # from Bison Pool metagenome (JGI, 2007)
   # read the file that connects taxids with the sequence identifier
   tfile <- system.file("extdata/gi.taxid.txt",package="CHNOSZ")
   gi.taxid <- scan(tfile,what=as.list(character(2)),flush=TRUE)
   # read the file that connects names with the taxids
   nfile <- system.file("data/taxid.phylum.csv",package="CHNOSZ")
   taxid.phylum <- read.csv(nfile)
   # the BLAST files
   bfile <- paste("extdata/bison",c("N","R","P"),"_vs_refseq39.blastp",sep="")
   for(i in 1:3) {
     file <- system.file(bfile[i],package="CHNOSZ")
     # process the blast file -- get taxon names
     b <- count.taxa(file,gi.taxid,taxid.phylum,min.taxon=2^-6)
     # count each of the phyla
     bd <- as.matrix(sapply(unique(b$phylum),function(x) (length(which(x==b$phylum)))))
     # make a matrix -- each column for a different file
     if(i==1) bardata <- bd else {
       bardata <- merge(bardata,bd,all=TRUE,by="row.names")
       rownames(bardata) <- bardata$Row.names
       bardata <- bardata[,-1]
     }
   }
   # normalize the counts
   bardata[is.na(bardata)] <- 0
   bardata <- t(t(bardata)/colSums(bardata))
   colnames(bardata) <- c("N","R","P")
   # make a bar chart
   bp <- barplot(as.matrix(bardata),col=rainbow(nrow(bardata)),
     xlab="location",ylab="fractional abundance")
   # add labels to the bars
   names <- substr(row.names(bardata),1,3)
   for(i in 1:3) {
     bd <- bardata[,i]
     ib <- bd!=0
     y <- (cumsum(bd) - bd/2)[ib]
     text(bp[i],y,names[ib])
   }
   title(main=paste("Most Abundant Phyla in a Portion",
     "of the Bison Pool Metagenome",sep="\n"))
}

\references{

  Altschul, S. F., Madden, T. L., Schaffer, A. A., Zhang, J. H., Zhang, Z., Miller, W. and Lipman, D. J., 1997. Gapped BLAST and PSI-BLAST: a new generation of protein database search programs. \emph{Nucleic Acids Res.}, 25, 3389-3402. \url{http://dx.doi.org/doi:10.1093/nar/25.17.3389} 


}

\keyword{misc}


