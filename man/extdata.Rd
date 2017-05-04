\encoding{UTF-8}
\name{extdata}
\alias{extdata}
\docType{data}
\title{Extra Data}

\description{
  The files in the subdirectories of \code{extdata} support the examples in the package documentation and vignettes. 
}

\details{

  Files in \code{abundance} contain protein abundance and microbial occurrence data:
  \itemize{
    \item \code{TBD+05.csv} lists genes with transcriptomic expression changes in carbon limitation stress response experiments in yeast (Tai et al., 2005). See \code{\link{read.expr}} for an example that uses this file.
    \item \code{ISR+08.csv} has columns excerpted from Additional File 2 of Ishihama et al. (2008) for protein abundances in \emph{E. coli} cytosol. The columns in this file are ID (Swiss-Prot ID), accession (Swiss-Prot accession), emPAI (exponentially modified protein abundance index), copynumber (emPAI-derived copy number/cell), GRAVY (Kyte-Doolittel), FunCat (FunCat class description), PSORT (PSORT localisation), ribosomal (yes/no). See \code{\link{read.expr}} for an example that uses this file.
    \item \code{yeastgfp.csv.xz} Has 28 columns; the names of the first five are \code{yORF}, \code{gene name}, \code{GFP tagged?}, \code{GFP visualized?}, and \code{abundance}. The remaining columns correspond to the 23 subcellular localizations considered in the YeastGFP project (Huh et al., 2003 and Ghaemmaghami et al., 2003) and hold values of either \code{T} or \code{F} for each protein. \samp{yeastgfp.csv} was downloaded on 2007-02-01 from http://yeastgfp.ucsf.edu using the Advanced Search, setting options to download the entire dataset and to include localization table and abundance, sorted by orf number. See \code{\link{yeastgfp}} and \code{demo("yeastgfp")} for examples that use this file.
    \item \code{microbes.csv} has data for microbial occurrence (i.e. relative enrichement) in colorectal cancer and normal tissue. The file is from the Supporting Information of Dick (2016). This file is used by \code{demo("bugstab")}.
  }

  Files in \code{bison} contain BLAST results and taxonomic information for an environmental metagenome from the Bison Pool hot spring in Yellowstone National Park:
  \itemize{
    \item \code{bisonN_vs_refseq57.blast.xz}, \code{bisonS...}, \code{bisonR...}, \code{bisonQ...}, \code{bisonP...} are partial tabular BLAST results for proteins in the Bison Pool Environmental Genome. Protein sequences predicted in the metagenome were downloaded from the Joint Genome Institute's IMG/M system on 2009-05-13. The target database for the searches was constructed from microbial protein sequences in National Center for Biotechnology Information (NCBI) RefSeq database version 57, representing 7415 microbial genomes. The \sQuote{blastall} command was used with the default setting for E value cuttoff (10.0) and options to make a tabular output file consisting of the top 20 hits for each query sequence. The function \code{\link{read.blast}} was used to extract only those hits with E values less than or equal to 1e-5 and with sequence similarity (percent identity) at least 30 percent, and to keep only the first hit for each query sequence. The function \code{\link{write.blast}} was used to save partial BLAST files (only selected columns). The files provided with CHNOSZ contain the first 5,000 hits for each sampling site at Bison Pool, representing between about 7 to 15 percent of the first BLAST hits after similarity and E value filtering.
    \item \code{gi.taxid.txt.xz} is a table that lists the sequence identifiers (gi numbers) that appear in the example BLAST files (see above), together with the corresponding taxon ids used in the NCBI databases. This file is \emph{not} a subset of the complete \sQuote{gi_taxid_prot.dmp.gz} available at \url{ftp://ftp.ncbi.nih.gov/pub/taxonomy/} but instead is a subset of \sQuote{gi.taxid.txt} generated from the RefSeq release catalog using \sQuote{gencat.sh} in the \code{refseq} directory. See \code{\link{id.blast}} for an example that uses this file and the BLAST files described above.
  }


  Files in \code{cpetc} contain heat capacity and other thermodynamic data and parameters:
  \itemize{
    \item \code{PM90.csv} Heat capacities of four unfolded aqueous proteins taken from Privalov and Makhatadze, 1990. Temperature in \eqn{^{\circ}}{°}C is in the first column, and heat capacities of the proteins in J mol\eqn{^{-1}}{^-1} K\eqn{^{-1}}{^-1} in the remaining columns. See \code{\link{ionize.aa}} and the vignette \code{anintro.Rmd} for examples that uses this file.
    \item \code{RH95.csv} Heat capacity data for iron taken from Robie and Hemingway, 1995. Temperature in Kelvin is in the first column, heat capacity in J K\eqn{^{-1}}{^-1} mol\eqn{^{-1}}{^-1} in the second. See \code{\link{subcrt}} for an example that uses this file.
    \item \code{RT71.csv} pH titration measurements for unfolded lysozyme (\samp{LYSC_CHICK}) taken from Roxby and Tanford, 1971. pH is in the first column, net charge in the second. See \code{\link{ionize.aa}} for an example that uses this file.
    \item \code{SOJSH.csv} Experimental equilibrium constants for the reaction NaCl(aq) = Na+ + Cl- as a function of temperature and pressure taken from Fig. 1 of Shock et al., 1992. Data were extracted from the figure using g3data (\url{http://www.frantz.fi/software/g3data.php}). See \code{demo("NaCl")} for an example that uses this file.
    \item \code{Cp.CH4.HW97.csv}, \code{V.CH4.HWM96.csv} Apparent molar heat capacities and volumes of CH4 in dilute aqueous solutions reported by Hnědkovský and Wood, 1997 and Hnědkovský et al., 1996. See \code{\link{EOSregress}} and the vignette \code{eos-regress.Rmd} for examples that use these files.
    \item \code{SC10_Rainbow.csv} Values of temperature (\eqn{^{\circ}}{°}C), pH and logarithms of activity of \CO2, \H2, \NH4plus, \H2S and \CH4 for mixing of seawater and hydrothermal fluid at Rainbow field (Mid-Atlantic Ridge), taken from Shock and Canovas, 2010. See the vignette \code{anintro.Rmd} for an example that uses this file.
    \item \code{SS98_Fig5a.csv}, \code{SS98_Fig5b.csv} Values of logarithm of fugacity of \eqn{\mathrm{O_2}}{O2} and pH as a function of temperature for mixing of seawater and hydrothermal fluid, digitized from Figs. 5a and b of Shock and Schulte, 1998. See the vignette \code{anintro.Rmd} for an example that uses this file.
    \item \code{rubisco.csv} UniProt IDs for Rubisco, ranges of optimal growth temperature of organisms, domain and name of organisms, and URL of reference for growth temperature, from Dick, 2014. See the vignette \code{anintro.Rmd} for an example that uses this file.
  }

  Files in \code{fasta} contain protein sequences:
  \itemize{
    \item \code{EF-Tu.aln} consists of aligned sequences (394 amino acids) of elongation factor Tu (EF-Tu). The sequences correspond to those taken from UniProtKB for ECOLI (\emph{Escherichia coli}), THETH (\emph{Thermus thermophilus}) and THEMA (\emph{Thermotoga maritima}), and reconstructed ancestral sequences taken from Gaucher et al., 2003 (maximum likelihood bacterial stem and mesophilic bacterial stem, and alternative bacterial stem). See \code{\link{read.fasta}} for an example that uses this file.
    \item \code{rubisco.fasta} Sequences of Rubisco obtained from UniProt (see Dick, 2014). See the vignette \code{anintro.Rmd} for an example that uses this file.
  }

  Files in \code{protein} contain amino acid compositions for proteins.
  See \code{\link{more.aa}} and \code{\link{read.expr}} for examples that use these files.
  \itemize{
    \item \code{Sce.csv.xz}
      Data frame of amino acid composition of 6716 proteins from the \emph{Saccharomyces} Genome Database (SGD).
      Values in the first three columns are the \code{ORF} names of proteins, \code{SGDID}, and \code{GENE} names. The remaining twenty columns (\code{ALA}..\code{VAL}) contain the numbers of the respective amino acids in each protein.
      The sources of data for \samp{Sce.csv} are the files \samp{protein_properties.tab} and \samp{SGD_features.tab} (for the gene names), downloaded from \url{http://www.yeastgenome.org} on 2013-08-24.
    \item \code{Eco.csv.xz}
      Amino acid compositions of 4407 proteins in \emph{Escherichia coli} strain K12.
      Format is the one used \code{\link{thermo}$protein}, with columns \samp{protein} holding the gene name, \samp{organism} set to \samp{ECOLI}, and \code{abbrv} holding the UniProt ID.
      The source of data is the file \samp{ECOLI.fas} downloaded from the HAMAP (High-quality Automated and Manual Annotation of microbial Proteomes system) FTP site (Gattiker et al., 2003) on 2010-09-25 (old URL: ftp://ftp.expasy.org/databases/hamap/complete_proteomes/fasta/bacteria).
    \item \code{DS11.csv}, \code{DS13.csv}
      These two files contain amino acid compositions of metagenomically encoded proteins, averaged together according to functional annotation (DS11) or taxonomic affiliation (DS13).
      The data are from Dick and Shock, 2011 and 2013.
      They are used in the vignette \emph{Hot-spring proteins in CHNOSZ}.
    \item \code{microbial.aa.csv}
      Overall protein compositions of microbial species reported to be positively or negatively enriched in colorectal cancer.
      This file is taken from Dick, 2016.
      It is used by \code{demo("bugstab")}.
  }

  Files in \code{refseq} contain code and results of processing NCBI Reference Sequences (RefSeq) for microbial proteins, using RefSeq release 61 of 2013-09-09:
  \itemize{
    \item \code{README.txt} Instructions for producing the data files.
    \item \code{gencat.sh} Bash script to extract microbial protein records from the RefSeq catalog.
    \item \code{gi.taxid.txt} Output from above. The complete file is too large to distribute with CHNOSZ, but a portion is included in \code{extdata/bison} to support processing example BLAST files for the Bison Pool metagenome (based on RefSeq 57, 2013-01-08).
    \item \code{mkfaa.sh} Combine the contents of .faa.gz files into a single FASTA file (to use e.g. for making a BLAST database).
    \item \code{protein.refseq.R} Calculate average amino acid composition of all proteins for each organism identified by a taxonomic ID.
    \item \code{trim_refseq.R} Keep only selected organism names (reduces number of taxa from 6758 to 779, helps to control package size).
    \item \code{protein_refseq.csv.xz} Output from above. See example in \code{\link{pinfo}}.
    \item \code{taxid.names.R} Generate a table of scientific names for the provided taxids. Requires the complete \code{names.dmp} and \code{nodes.dmp} from NCBI taxonomy files.
    \item \code{taxid_names.csv.xz} Output from above.
      NOTE: For backward compatibility with the example BLAST files for the Bison Pool metagenome, the packaged file merges records for taxids found in either RefSeq 57 or 61.
      Certain taxids in release 57 were not located in the current RefSeq catalog, probably related to the transition to the \dQuote{WP} multispecies accessions (\url{ftp://ftp.ncbi.nlm.nih.gov/refseq/release/announcements/WP-proteins-06.10.2013.pdf}).
      See example for \code{\link{id.blast}}.
  }

  Files in \code{supcrt} contain scripts for reading and comparing SUPCRT files (including slop98.dat and newer slop files from GEOPIG (http://geopig.asu.edu)) with the database in CHNOSZ:
  \itemize{
    \item \code{read.supcrt.R} defines the function \code{read.supcrt} that can be used to read SUPCRT files.
    \item \code{compare.R} uses \code{read.supcrt} to compare data in the SUPCRT file with that in \code{thermo$obigt}.
    \item \code{newnames.csv} maps names generated by \code{read.supcrt}, based on names present in the source SUPCRT files, to names used in \code{thermo$obigt}.
  }

  Files in \code{taxonomy} contain taxonomic data files:
  \itemize{
    \item \code{names.dmp} and \code{nodes.dmp} are excerpts of the taxonomy files available on the NCBI ftp site (\url{ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz}, accessed 2010-02-15). These files contain only the entries for \emph{Escherichia coli} K-12, \emph{Saccharomyces cerevisiae}, \emph{Homo sapiens}, \emph{Pyrococcus furisosus} and \emph{Methanocaldococcus jannaschii} (taxids 83333, 4932, 9606, 186497, 243232) and the higher-ranking nodes (genus, family, etc.) in the respective lineages. See \code{\link{taxonomy}} for examples that use these files.
  }

  Files in \code{thermo} contain additional thermodynamic data and group additivity definitions:
  \itemize{
    \item \code{Ste01.csv} and \code{BZA10.csv} contain supplementary thermodynamic data taken from Stefansson (2001) and Bazarkina et al. (2010). They can be added to the database in the current session using \code{\link{add.obigt}}. See \code{\link{add.obigt}} and the vignette \code{eos-regress.Rmd} for examples that use these files.
    \item \code{obigt_check.csv} contains the results of running \code{\link{check.obigt}} to check the internal consistency of entries in the primary and supplementary databases.
    \item \code{RH98_Table15.csv} Group stoichiometries for high molecular weight crystalline and liquid organic compounds taken from Table 15 of Richard and Helgeson, 1998. The first three columns have the \code{compound} name, \code{formula} and physical \code{state} (\samp{cr} or \samp{liq}). The remaining columns have the numbers of each group in the compound; the names of the groups (columns) correspond to species in \code{\link{thermo}$obigt}. The compound named \samp{5a(H),14a(H)-cholestane} in the paper has been changed to \samp{5a(H),14b(H)-cholestane} here to match the group stoichiometry given in the table. See \code{\link{RH2obigt}} for a function that uses this file.
    \item \code{DLEN67.csv} Standard Gibbs energies of formation, in kcal/mol, from Dayhoff et al., 1967, for nitrogen (N2) plus 17 compounds shown in Fig. 2 of Dayhoff et al., 1964, at 300, 500, 700 and 1000 K. See \code{demo("wjd")} and the vignette \code{wjd.Rmd} for examples that use this file.
  }


}

\references{
Bazarkina, E. F., Zotov, A. V. and Akinfiev, N. N. (2010) Pressure-dependent stability of cadmium chloride complexes: Potentiometric measurements at 1–1000 bar and 25°C. \emph{Geol. Ore Deposits} \bold{52}, 167--178. \url{https://doi.org/10.1134/S1075701510020054}

Dayhoff, M. O. and Lippincott, E. R. and Eck, R. V. (1964) Thermodynamic Equilibria In Prebiological Atmospheres. \emph{Science} \bold{146}, 1461--1464. \url{https://doi.org/10.1126/science.146.3650.1461}

Dayhoff, M. O. and Lippincott, E. R., Eck, R. V. and Nagarajan (1967) Thermodynamic Equilibrium In Prebiological Atmospheres of C, H, O, N, P, S, and Cl. Report SP-3040, National Aeronautics and Space Administration. \url{http://ntrs.nasa.gov/search.jsp?R=19670017966}

Dick, J. M. (2014) Average oxidation state of carbon in proteins. \emph{J. R. Soc. Interface} \bold{11}, 20131095. \url{https://doi.org/10.1098/rsif.2013.1095}

Dick, J. M. (2016) Proteomic indicators of oxidation and hydration state in colorectal cancer. \emph{PeerJ} \bold{4}:e2238. \url{https://doi.org/10.7717/peerj.2238}

Dick, J. M. and Shock, E. L. (2011) Calculation of the relative chemical stabilities of proteins as a function of temperature and redox chemistry in a hot spring. \emph{PLoS ONE} \bold{6}, e22782. \url{https://doi.org/10.1371/journal.pone.0022782}

Dick, J. M. and Shock, E. L. (2013) A metastable equilibrium model for the relative abundance of microbial phyla in a hot spring. \emph{PLoS ONE} \bold{8}, e72395. \url{https://doi.org/10.1371/journal.pone.0072395}

Gattiker, A., Michoud, K., Rivoire, C., Auchincloss, A. H., Coudert, E., Lima, T., Kersey, P., Pagni, M., Sigrist, C. J. A., Lachaize, C., Veuthey, A.-L., Gasteiger, E. and Bairoch, A. (2003) Automatic annotation of microbial proteomes in Swiss-Prot. \emph{Comput. Biol. Chem.} \bold{27}, 49--58. \url{https://doi.org/10.1016/S1476-9271(02)00094-4}

Gaucher, E. A., Thomson, J. M., Burgan, M. F. and Benner, S. A (2003) Inferring the palaeoenvironment of ancient bacteria on the basis of resurrected proteins. \emph{Nature} \bold{425}(6955), 285--288. \code{https://doi.org/10.1038/nature01977}

Ghaemmaghami, S., Huh, W., Bower, K., Howson, R. W., Belle, A., Dephoure, N., O'Shea, E. K. and Weissman, J. S. (2003) Global analysis of protein expression in yeast. \emph{Nature} \bold{425}(6959), 737--741. \code{https://doi.org/10.1038/nature02046}

Huh, W. K., Falvo, J. V., Gerke, L. C., Carroll, A. S., Howson, R. W., Weissman, J. S. and O'Shea, E. K. (2003) Global analysis of protein localization in budding yeast. \emph{Nature} \bold{425}(6959), 686--691. \code{https://doi.org/10.1038/nature02026}

HAMAP system. HAMAP FTP directory, \url{ftp://ftp.expasy.org/databases/hamap/}

Hnědkovský, L., Wood, R. H. and Majer, V. (1996) Volumes of aqueous solutions of CH4, CO2, H2S, and NH3 at temperatures from 298.15 K to 705 K and pressures to 35 MPa. \emph{J. Chem. Thermodyn.} \bold{28}, 125--142. \url{https://doi.org/10.1006/jcht.1996.0011}

Hnědkovský, L. and Wood, R. H. (1997) Apparent molar heat capacities of aqueous solutions of CH4, CO2, H2S, and NH3 at temperatures from 304 K to 704 K at a pressure of 28 MPa. \emph{J. Chem. Thermodyn.} \bold{29}, 731--747. \url{https://doi.org/10.1006/jcht.1997.0192}

Ishihama, Y., Schmidt, T., Rappsilber, J., Mann, M., Hartl, F. U., Kerner, M. J. and Frishman, D. (2008) Protein abundance profiling of the \emph{Escherichia coli} cytosol. \emph{BMC Genomics} \bold{9}:102. \url{https://doi.org/10.1186/1471-2164-9-102}

Joint Genome Institute (2007) Bison Pool Environmental Genome. Protein sequence files downloaded from IMG/M (\url{http://img.jgi.doe.gov/cgi-bin/m/main.cgi?section=FindGenomes&page=findGenomes})

Privalov, P. L. and Makhatadze, G. I. (1990) Heat capacity of proteins. II. Partial molar heat capacity of the unfolded polypeptide chain of proteins: Protein unfolding effects. \emph{J. Mol. Biol.} \bold{213}, 385--391. \url{https://doi.org/10.1016/S0022-2836(05)80198-6}

Richard, L. and Helgeson, H. C. (1998) Calculation of the thermodynamic properties at elevated temperatures and pressures of saturated and aromatic high molecular weight solid and liquid hydrocarbons in kerogen, bitumen, petroleum, and other organic matter of biogeochemical interest. \emph{Geochim. Cosmochim. Acta} \bold{62}, 3591--3636. \url{https://doi.org/10.1016/S0016-7037(97)00345-1}

Robie, R. A. and Hemingway, B. S. (1995) \emph{Thermodynamic Properties of Minerals and Related Substances at 298.15 K and 1 Bar (\eqn{10^5} Pascals) Pressure and at Higher Temperatures}. U. S. Geol. Surv., Bull. 2131, 461 p. \url{http://www.worldcat.org/oclc/32590140}

Roxby, R. and Tanford, C. (1971) Hydrogen ion titration curve of lysozyme in 6 M guanidine hydrochloride. \emph{Biochemistry} \bold{10}, 3348--3352. \url{https://doi.org/10.1021/bi00794a005}

SGD project. \emph{Saccharomyces} Genome Database, \url{http://www.yeastgenome.org}

Shock, E. and Canovas, P. (2010) The potential for abiotic organic synthesis and biosynthesis at seafloor hydrothermal systems. \emph{Geofluids} \bold{10}, 161--192. \url{https://doi.org/10.1111/j.1468-8123.2010.00277.x}

Shock, E. L., Oelkers, E. H., Johnson, J. W., Sverjensky, D. A. and Helgeson, H. C. (1992) Calculation of the thermodynamic properties of aqueous species at high pressures and temperatures: Effective electrostatic radii, dissociation constants and standard partial molal properties to 1000 \eqn{^{\circ}}{°}C and 5 kbar. \emph{J. Chem. Soc. Faraday Trans.} \bold{88}, 803--826. \url{https://doi.org/10.1039/FT9928800803}

Shock, E. L. and Schulte, M. D. (1998) Organic synthesis during fluid mixing in hydrothermal systems. \emph{J. Geophys. Res.} \bold{103}, 28513--28527. \url{https://doi.org/10.1029/98JE02142}

Stefánsson, A. (2001) Dissolution of primary minerals of basalt in natural waters. I. Calculation of mineral solubilities from 0°C to 350°C \emph{Chem. Geol.} \bold{172}, 225--250. \url{https://doi.org/10.1016/S0009-2541(00)00263-1}

Tai, S. L., Boer, V. M., Daran-Lapujade, P., Walsh, M. C., de Winde, J. H., Daran, J.-M. and Pronk, J. T. (2005) Two-dimensional transcriptome analysis in chemostat cultures: Combinatorial effects of oxygen availability and macronutrient limitation in \emph{Saccharomyces cerevisiae}. \emph{J. Biol. Chem.} \bold{280}, 437--447. \url{https://doi.org/10.1074/jbc.M410573200}

YeastGFP project. Yeast GFP Fusion Localization Database, http://yeastgfp.ucsf.edu; Current location: \url{http://yeastgfp.yeastgenome.org}
}

\concept{Thermodynamic data}
