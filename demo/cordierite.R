### 1-D property plot
## for hydrous cordierite = cordierite + H2O
## after Helgeson et al., 1978
## (Summary and critique of the thermodynamic properties of 
## rock-forming minerals. Am. J. Sci., 278-A, 1-229)
basis(c("cordierite,hydrous","Mg+2","SiO2","H2O","O2","H+"))
species("cordierite")
# water.SUPCRT92 can only get us up to 5000 bar
# (lines for 7000 and 10000 bar are in the original diagram)
P <- c(1,2,3,5)*1000
col <- rainbow(length(P))
for(i in 1:length(P)) {
  a <- affinity(property="logK",T=c(20,800),P=P[i])
  diagram(a,add=(i!=1),ylim=c(-4,2),legend.x=NULL,
    col=col[i],main="")
}
legend("topright",lty=1,col=col,legend=paste(P,"bar"))
title(main=paste("hydrous cordierite = cordierite + H2O",
  "After Helgeson et al., 1978",sep="\n"),cex.main=0.9)
