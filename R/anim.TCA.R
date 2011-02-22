# CHNOSZ/anim.TCA.R
# functions to make various animations

anim.TCA <- function(redox=list(O2=c(-95,-60)),high.T=FALSE,
  nframes.in=15,nframes.out=20,nframes.anim=140,pHlim=c(0,10),
  width = 768, height = 576) {
  # tricarboxylic acids 20100701, 20110204 jmd
  # we depend on an empty png directory
  if(!"png" %in% dir()) stop("directory 'png' not present")
  else if(length(dir("png")) > 0) stop("directory 'png' not empty")
  # initialize the system
  # add supplementary data (from default location of data/OBIGT-2.csv)
  # which includes properties from 
  add.obigt()
  # expand default logfO2 range if we're at high temperature
  if(high.T & missing(redox)) redox <- list(O2=c(-100,-40))
  # the name of 'redox' should be either O2 or H2
  basis(c("CO2","H2O",names(redox),"H+"))
  # load the species in order of increasing degree of ionization
  species(c("pyruvic acid","oxaloacetic acid","malic acid", # Z = 0
    "fumaric acid","a-ketoglutaric acid","citric acid",
    "pyruvate","H-oxaloacetate","H-malate","H-fumarate",    # Z = -1
    "H-a-ketoglutarate","H2-citrate",
    "oxaloacetate-2","malate-2","fumarate-2",               # Z = -2
    "a-ketoglutarate-2","H-citrate-2",
    "citrate-3"),rep("aq",18))                              # Z = -3
  # start the plot device - multiple png figures
  png(filename="png/Rplot%04d.png",width=width,height=height)
  par(mar=c(3,3.5,1.5,1),mgp=c(1.9,1,0))
  # leadin/out frames
  # the variable(s) that change with each frame
  pH <- seq(pHlim[1],pHlim[2],length.out=nframes.anim)
  pHres <- length(pH)
  pHlim <- range(pH)
  baseres <- 128
  # calculate the affinity
  redox[[1]] <- c(redox[[1]],baseres)
  affargs <- c(redox,list(H2O=c(-10,10,baseres),pH=c(pHlim,pHres)))
  a.cold <- do.call(affinity,affargs)
  if(high.T) {
    affargs <- c(affargs,list(T=100))
    a.hot <- do.call(affinity,affargs)
  }
  # put the movie together
  iframes <- c(rep(1,nframes.in),1:length(pH),rep(pHres,nframes.out))
  for(i in 1:length(iframes)) {
    print(paste("frame",i,"of",length(iframes),"at pH",pH[iframes[i]]))
    # don't fill fields with heat colors if we're doing a high-T overlay
    if(high.T) color <- "lightgrey" else color <- "heat"
    # take a single O2-H2O slice along pH
    diagram(slice.affinity(a.cold,3,iframes[i]),color=color,cex=1.5)
    ltext <- substitute(paste(italic(T)==x~degree*C),list(x=25))
    lcol <- "black"
    # overlay the high-T diagram
    if(high.T) {
      ltext <- c(ltext,substitute(paste(italic(T)==x~degree*C),list(x=100)))
      lcol <- c(lcol,"red")
      diagram(slice.affinity(a.hot,3,iframes[i]),add=TRUE,col="red",cex=1.5)
    }
    legend("topleft",legend=as.expression(ltext),lty=1,col=lcol)
    title(main=paste("TCA cycle reactants; pH =",
      format(round(pH[iframes[i]],2))))
  }
  # close the plot device - convert to movie
  dev.off()
  cat("anim.TCA: attempting conversion to animated GIF...\n")
  outfile <- "TCA.gif"
  # convert tool from imagemagick.org; -loop 0 for infinite loop
  #syscmd <- paste("convert -loop 2 -delay 12 png/*.png png/",outfile,sep="")
  syscmd <- paste("convert -loop 0 -delay 10 png/*.png png/",outfile,sep="")
  cat(paste(syscmd,"\n"))
  sres <- system(syscmd)
  if(sres==0) cat(paste("anim.TCA: output is in png/",outfile,"\n",sep=""))
  else cat("anim.TCA: error converting to animated GIF\n")
  #print("removing png files")
  #system("rm png/*.png")
}

