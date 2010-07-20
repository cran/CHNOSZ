# CHNOSZ/examples.R
# Copyright (C) 2007 Jeffrey M. Dick

examples <- function(do.png=FALSE) {
  # run all the examples in CHNOSZ documentation
  .ptime <- proc.time()
  # "hkf" refers to "eos"
  topics <- c("CHNOSZ","util.args","util.array","util.blast","util.character","util.data",
    "util.fasta","util.formula","util.misc","util.seq","util.units",
    "taxonomy","info","hkf","water","subcrt","examples","makeup","basis","species","affinity",
    "diagram","buffer","protein","ionize","get.protein","revisit","transfer")
  do.plot <- FALSE
  if(is.character(do.png))
    png(paste(do.png,"%d.png",sep=""),width=700,height=700,pointsize=18)
  else if(do.png) do.plot <- TRUE
  for(i in 1:length(topics)) {
    if(do.plot) png(paste(topics[i],"%d.png",sep=""),width=700,height=700,pointsize=18)
    myargs <- list(topic=topics[i],ask=FALSE)
    do.call(example,myargs)
    if(do.plot) dev.off()
  }
  if(is.character(do.png)) dev.off()
  cat("Time elapsed: ", proc.time() - .ptime, "\n")
}


longex <- function(which) {
  # extra examples for fun
  #if(is.null(which)) which <- c("sources","copper","cordierite","phosphate","nucleobase","pie")
  if(which=="sources") {
    ## cross-checking sources
    # the reference sources
    ref.source <- thermo$source$source
    # only take those that aren't journal abbreviations
    ref.source <- ref.source[-grep('_',ref.source)]
    # sources of elemental data
    element.source <- thermo$element$source
    # sources in the primary thermodynamic database
    os1 <- thermo$obigt$source1
    os2 <- thermo$obigt$source2
    # sources also in the supplemental database
    danger()
    os3 <- thermo$obigt$source1
    os4 <- thermo$obigt$source2
    data(thermo)
    # all of the thermodynamic data sources - some of them might be NA
    obigt.source <- unique(c(os1,os2,os3,os4))
    obigt.source <- obigt.source[!is.na(obigt.source)]
    # sources of protein compositions
    protein.source <- thermo$protein$source
    # sources of stress response proteins
    stress.source <- as.character(thermo$stress[2,])
    # if the sources are all accounted for 
    # these all produce character(0)
    print("missing these sources for elemental properties:")
    print(unique(element.source[!(element.source %in% ref.source)]))
    print("missing these sources (1) for thermodynamic properties:")
    print(unique(obigt.source[!(obigt.source %in% ref.source)]))
    print("missing these sources for protein compositions:")
    print(unique(protein.source[!(protein.source %in% ref.source)]))
    print("missing these sources for stress response experiments:")
    print(unique(stress.source[!(stress.source %in% ref.source)]))
    # determine if all the reference sources are cited
    my.source <- c(element.source,obigt.source,protein.source,stress.source)
    # this should produce character(0)
    print("these sources are present but not cited:")
    print(ref.source[!(ref.source %in% my.source)])
  } else if(which=="copper") {
    ## Eh-pH diagrams for copper-water-glycine
    ## After Fig. 2 of Aksu and Doyle, 2001
    ## (Aksu, S. and Doyle, F. M., 2001. Electrochemistry of copper in aqueous glycine 
    ## solutions. J. Electrochem. Soc., 148, B51-B57. doi:10.1149/1.1344532)
    ##  We need to add some species and change some Gibbs energies.
    # update rows of the database
    i <- info(c("Cu(Gly)+","glycinate","e-","H+"))
    n <- nrow(thermo$obigt <<- rbind(thermo$obigt,thermo$obigt[rep(i[1],2),]))
    thermo$obigt$name[n-1] <<- "Cu(Gly)2-"
    thermo$obigt$name[n] <<- "HCu(Gly)+2"
    thermo$obigt$formula[n-1] <<- makeup(makeup(c(i[1],i[2],i[3]),""),"")
    thermo$obigt$formula[n] <<- makeup(makeup(c(i[1],i[4]),""),"")
    # In Fig 2b, total log activities of Cu (Cu_T) 
    # and glycine (L_T) are -4 and -1
    basis(c("Cu+2","H2O","H+","e-","glycine","CO2"),c(999,0,999,999,-1,999))
    # solid species
    species(c("copper","cuprite","tenorite"))
    # aqueous species
    species(c("glycinium","glycine","glycinate","Cu+","Cu+2","CuO2-2","HCuO2-",
      "Cu(Gly)+","Cu(Gly)2","Cu(Gly)2-","HCu(Gly)+2"),-4)
    ispecies <- species()$ispecies
    # update the Gibbs energies using A&D's Table 1 and Table II
    logK <- c(convert(convert(c(0,-146,-129.7,-384.061,-370.647,-314.833,
      49.98,65.49,-183.6,-258.5,-298.2)*1000,"cal"),"logK"),15.64,10.1,2.92) 
    # do it in order so later species take account of prev. species' values
    for(i in 1:length(logK)) {
      G <- convert(logK[i],"G")
      if(i==12) G <- G + thermo$obigt$G[ispecies[8]] + 
        2*thermo$obigt$G[ispecies[6]]
      if(i==13) G <- G + thermo$obigt$G[ispecies[7]] + 
        2*thermo$obigt$G[ispecies[6]]
      if(i==14) G <- G + thermo$obigt$G[ispecies[11]]
      thermo$obigt$G[ispecies[i]] <- G
    }  # done with changing Gibbs free energies!
    # we have to get some leftovers out of there or diagram() gets confused
    species(c("glycinium","glycine","glycinate"),delete=TRUE)
    # make a plot to see if it's working
    ispecies <- ispecies[-(1:6)]
    afun <- function(cu,gly) {
      # from above: our fifth basis species is glycine(-ate,-ium)
      basis(rownames(basis())[5],gly)
      t <- match(ispecies,species()$ispecies)
      species(t,cu)
      affinity(pH=c(0,16),Eh=c(-0.6,1.0))
    }
    diagram(afun(-4,-1))
    title(main=paste("Aqueous Copper + Glycine, 25 deg C, 1 bar",
      "After Aksu and Doyle, 2001 Fig. 2b",sep="\n"))
    # What's missing? Try glycinate not glycine in reactions at ph > ~9.8
    basis(c("Cu+2","H2O","H+","e-","glycinate","CO2"),
      c(999,0,999,999,-2,999))
    species(c("copper","cuprite","tenorite","Cu+","Cu+2","CuO2-2","HCuO2-",
      "Cu(Gly)+","Cu(Gly)2","Cu(Gly)2-","HCu(Gly)+2"))
    loga_Cu <- -4
    loga_Gly <- -1
    diagram(afun(loga_Cu,loga_Gly),color=NULL,col="blue",
      names=species()$name,col.names="blue",add=TRUE)
    water.lines()
    # the glycine ionization constants could be calculated using
    # subcrt, here they are taken from A&D Table II
    abline(v=c(2.35,9.778),lty=3)
    # now put glycinium (low pH) in the basis
    basis(c("Cu+2","H2O","H+","e-","glycinium","CO2"),c(999,0,999,999,-2,999))
    species(c("copper","cuprite","tenorite","Cu+","Cu+2","CuO2-2","HCuO2-",
      "Cu(Gly)+","Cu(Gly)2","Cu(Gly)2-","HCu(Gly)+2"))
    diagram(afun(loga_Cu,loga_Gly),color=NULL,col="green",
      names=NULL,col.names="green",add=TRUE)
    # let's forget our changes to 'thermo' so the example 
    # below that uses glycine will work as expected
    data(thermo)
  } else if(which=="cordierite") {
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
        col=col[i],title="")
    }
    legend("topright",lty=1,col=col,legend=paste(P,"bar"))
    title(main=paste("hydrous cordierite = cordierite + H2O",
      "After Helgeson et al., 1978",sep="\n"),cex.main=0.9)
  } else if(which=="phosphate") {
    ## speciation of phosphate as a function of ionic strength
    basis("CHNOPS+")
    T <- c(25,100)
    species(c("HPO4-2","H2PO4-"))
    diagram(affinity(IS=c(0,0.25),T=T[1]),ylim=c(-3.2,-2.8))
    title(main=paste("Aqueous phosphate speciation, pH=7",
      "25 and 100 deg C - black and red lines",sep="\n"))
    diagram(affinity(IS=c(0,0.25),T=T[2]),ylim=c(-3.2,-2.8),add=TRUE,col="red")  
    ## phosphate predominance f(IS,pH)
    diagram(affinity(pH=c(6.8,7.2),IS=c(0,0.25),T=T[1]))
    title(main=paste("Aqueous phosphate predominance, 25 deg C",
      "and 100 deg C (dotted overlay)",sep="\n"))
    diagram(affinity(pH=c(6.8,7.2),IS=c(0,0.25),T=T[2]),add=TRUE,dotted=2,
      names=NULL)
  } else if(which=="nucleobase") {
    ## Nucleobase - Amino Acid Interaction Eh-H2O
    # for this example we try a unique basis definition
    basis(c("CO2","H2O","glutamine","e-","H+"),c(-3,0,-3,0,-7))
    species(c("uracil","cytosine","adenine","guanine",
      "phenylalanine","proline","lysine","glycine"),"aq")
    # this loaded four nucleobases and four related amino acids
    # (coded for by the homocodon triplets)
    # check out the predominance diagrams
    a.1 <- affinity(H2O=c(-5,0),Eh=c(-0.5,0))
    diagram(a.1,color=NULL)
    # overlay a different temperature
    a.2 <- affinity(H2O=c(-5,0),Eh=c(-0.5,0),T=100)
    diagram(a.2,col="red",add=TRUE,names=NULL)
    # start make a title for the plot
    tb <- thermo$basis   # includes activities of basis species
    # exclude those that are on the axes
    tb <- tb[!((rownames(tb) %in% c("e-","H2O"))),]
    title(main=paste("Nucleobases and amino acids,",
      "T=25 and 100 C at Psat\n",
      describe(tb,T=NULL,P=NULL)),cex.main=0.9)
  } else if(which=="pie") {
    #    This is an attempt to predict some characteristics of the relative
    #    abundances of organisms that are depicted in the pie charts shown by
    #    Spear et al., 2005 (Spear, J. R., Walker, J. J., McCollom, T. M. and
    #    Pace, N. R. Hydrogen and bioenergetics in the Yellowstone geothermal
    #    ecosystem. Proc. Natl. Acad. Sci. U. S. A., 102, 2555-2560, 2005.).
    #    For each type of organism present, we use a single model protein. We
    #    take the reported relative abundances of the four most abundant
    #    organisms to generate an assemblage of proteins that buffers the
    #    activies of CO2, H2O and NH3. Then we calculate relative abundances of
    #    all proteins at different values of oxygen fugacity and H2S activity
    #    and display them on pie charts that can be compared with organismal
    #    abundances. The results show that it is possible to predict conditions
    #    that foster high (many slices in the pie diagrams) or low (few slices)
    #    diversity, even if the make up of the community is not perfectly
    #    replicated.
    ### Setup
    # Bacterial species and abundances from Spear et al., 2005 (Fig. 4):
    O1 <- c('Aquificales','Thermotogales','Bacillus','Thermodesulfobacteria')
    O2 <- c('Thermus/Deinococcus','Hydrogenobacter', 'Hydrogenobaculum','Bacteroidetes')
    O3 <- c('a-proteobacterium','d-proteobacterium','g-proteobacterium','OD-1')
    ORGANISMS <- c(O1,O2,O3) 
    # Which species are most abundant in low and high sulfur environments; 
    # the first three of each are aquificales, thermotogales, thermodesulfobacteria:
    lowS.which <- c(1,2,4,3,5)
    lowS.values <- c(75,12,4,7,2)
    highS.which <- c(1,2,4,11,6,9,7,8,10,13)
    highS.values <- c(63,2,10,2,5,2,9,2,3,1)
    # What chemical species we use to represent each organism:
    P1 <- c('ACCA_AQUAE','A8F7V4_THELT','RBL1_THIDA')
    P2 <- c('Q93EV7_9BACT','ACCA_DEIRA','Q05KD0_HYDTH','A7WGI1_9AQUI')
    P3 <- c('Q64VW6_BACFR','ACCA_CAUCR','A1VC70_DESVV','ACCA_PSEAE')
    PROTEINS <- c(P1,P2,P3)
    ### Function definitions
    # To make pie charts for calculated abundances:
    plot.pie.calc <- function(which="high",T=25,main='') {
      # first clean up the buffer definition in case we have
      # been run already
      thermo$buffer <<- thermo$buffer[thermo$buffer$name!='PROTEINS',]
      # we take four predominant proteins from SWM+05
      myprot <- PROTEINS[get(paste(which,"S.which",sep=""))][1:4]
      mypercent <- get(paste(which,"S.values",sep=""))[1:4]
      # use these four proteins to create a buffer
      mybufprot <- paste(myprot,'RESIDUE',sep='.')
      mod.buffer('PROTEINS',mybufprot,'aq',log10(mypercent/100))
      # our species are the residues of all proteins
      species(delete=TRUE)
      add.protein(protein.residue(PROTEINS))
      species(paste(PROTEINS,"RESIDUE",sep="."))
      # assign the buffer to three basis species
      basis(c('CO2','H2O','NH3'),'PROTEINS')
      # calculate the buffered activities
      a <- affinity(return.buffer=TRUE,balance=1,T=T)
      # make the titles
      sub <- c2s(paste(names(a)[1:3],round(as.numeric(a)[1:3])))
      main <- paste('\n',main,'\n',sub)
      # set the total species activities to those in the buffer
      species(1:nrow(species()),-99)
      species(mybufprot,log10(mypercent/100))
      # get the activities back to numeric
      basis(names(a)[1:3],as.numeric(a)[1:3])
      thermo$basis$logact <<- as.numeric(thermo$basis$logact)
      # colors
      col <- rep('white',99)
      col[match(myprot,PROTEINS)] <- heat.colors(4)
      # calculate the distribution of species
      mylogaH2O <- thermo$basis$logact[rownames(thermo$basis)=='H2O']
      a <- affinity(H2O=c(mylogaH2O,mylogaH2O-1,2),T=T)
      logacts <- diagram(a,do.plot=FALSE,residue=TRUE,balance=1)$logact
      # assemble the names and logarithms of activities
      # of species that are above a minimum value
      names <- character()
      values <- numeric()
      cols <- character()
      logactmin <- -2
      for(i in 1:length(logacts)) {
        myvalue <- logacts[[i]][1]
        if(myvalue > logactmin) {
          names <- c(names,ORGANISMS[i])
          values <- c(values,myvalue)
          cols <- c(cols,col[i])
        }
      }
      # remove the logarithms
      values <- 10^values
      # sort them by abundance
      isort <- sort(values,index.return=TRUE,decreasing=TRUE)$ix
      names <- names[isort]
      values <- values[isort]
      cols <- cols[isort]
      # make a pie chart
      pie(values,names,clockwise=FALSE,main=main,col=cols,radius=0.7)
    }
    # To plot pie charts for observed abundances of organisms:
    plot.pie.obs <- function(which="low") {
      # the values from SWM+05
      names <- ORGANISMS[get(paste(which,"S.which",sep=""))]
      values <- get(paste(which,"S.values",sep=""))
      main <- paste("observed at",which,"H2S")
      # colors for the four dominant species
      mycol <- heat.colors(4)
      # colors for the rest
      mycol <- c(mycol,rep('white',length(names)-length(mycol)))
      # sort the values
      isort <- sort(values,index.return=TRUE,decreasing=TRUE)$ix
      values <- values[isort]
      names <- names[isort]
      mycol <- mycol[isort]
      pie(values,names,clockwise=FALSE,main=main,col=mycol,radius=0.7)
    }
    # To plot both types of pie diagrams (showing calculated protein 
    # activities and observed abundances of organisms) at a given temperature:
    plot.pie <- function(T=80) {
      # first deal with the layout
      layout(matrix(1:4,byrow=TRUE,nrow=2))
      par(mar=c(1,1,1,1))
      # basis definition
      basis(c('CO2','H2O','NH3','H2S','hydrogen','H+'))
      basis('pH',7)
      val <- function(text)
        round(as.numeric(thermo$basis$logact[rownames(thermo$basis)==text]))
      # now to plotting
      # low sulfur and relatively oxidizing
      basis(c('H2','H2S'),c(-9,-9))
      plot.pie.calc("low",T=T,main=paste("calculated for H2",val("H2"),"H2S",val("H2S")))
      label.plot('a')
      plot.pie.obs("low")
      label.plot('b')
      # high sulfur and relatively reducing
      basis(c('H2','H2S'),c(-5,-3))
      plot.pie.calc("high",T=T,main=paste("calculated for H2",val("H2"),"H2S",val("H2S")))
      label.plot('c')
      plot.pie.obs("high")
      label.plot('d')
    }
    ### Now run it!
    # at 80 degrees C:
    plot.pie(80)
  }
}
