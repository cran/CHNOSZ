\encoding{UTF-8}
\name{extdata}
\alias{extdata}
\docType{data}
\title{Extra Data}

\description{
  The files in the subdirectories of \code{extdata} support the examples in the package documentation and vignettes. 
}

\details{

  Files in \code{abundance} contain protein abundance data:
  \itemize{
    \item \code{stress} is a data frame listing proteins identified in selected proteomic stress response experiments. The names of proteins begin at row 3, and columns are all the same length (padded as necessary at the bottom by \code{NA}s). Names correspond to ordered locus names (for \samp{Sce}) or gene names (for \samp{Eco}). The column names identify the experiments, the first row contains the name of the organism (\samp{Sce} or \samp{Eco}) and the third row has the reference key for the source of the data (listed in \code{\link{thermo}$refs}).
    \item \code{AA03.csv} has reference abundances for 71 proteins taken from Fig. 3 of Anderson and Anderson, 2002 (as corrected in Anderson and Anderson, 2003). The columns with data taken from these sources are type (hemoglobin, plasma, tissue, or interleukin), description (name used in the original figure), log10(pg/ml) (\emph{upper limit} of abundance interval shown in Anderson and Anderson, 2003, log10 of concentration in pg/ml). The additional columns are data derived from a search of the SWISS-PROT/UniProtKB database based on the descriptions of the proteins: name (nominal UniProtKB name for this protein), name2 (other UniProtKB names(s) that could apply to the protein), and note (notes based on searching for a protein of this description). The amino acid compositions of all proteins whose names are not NA are included in \code{thermo$protein}. The \code{abbrv} column for the proteins contains the description given by Anderson and Anderson, 2003, followed by (in parentheses) the UniProtKB accession number. Annotated initiator methionines (e.g. for ferritin, myoglobin, ENOG), signal peptides or propeptides were removed from the proteins (except where they are not annotated in UniProtKB: IGHG1, IGHA1, IGHD, MBP). In cases were multiple isoforms are present in UniProtKB (e.g. Albumin) only the first isoform was taken. In the case of C4 Complement (CO4A) and C5 Complement (CO5), the amino acid composition of only the alpha chains are listed. In the case of the protein described as iC3b, the amino acid sequence is taken to be that of Complement C3c alpha' chain fragment 1 from CO3, and is given the name CO3.C3c. The non-membrane (soluble) chains of TNF-binding protein (TNR1A) and TNF-alpha (TNFA) were used. Rantes, MIP-1 beta and MIP-1 alpha were taken from C-C motif chemokines (CCL5, CCL4, CCL3 respectively). C-peptide was taken from the corresponding annotation for insulin and here is named INS.C. See \code{\link{protein}} and \code{\link{read.expr}} for examples that use this file.
    \item \code{ISR+08.csv} has columns excerpted from Additional File 2 of Ishihama et al. (2008) for protein abundances in \emph{E. coli} cytosol. The columns in this file are ID (Swiss-Prot ID), accession (Swiss-Prot accession), emPAI (exponentially modified protein abundance index), copynumber (emPAI-derived copy number/cell), GRAVY (Kyte-Doolittel), FunCat (FunCat class description), PSORT (PSORT localisation), ribosomal (yes/no). See \code{\link{read.expr}} for examples that use this file.
%%      \item \code{GLL+98.csv} has columns "oln" for ordered locus name and "ratio" for change in expression of yeast proteins in response to H2O2 treatment, from Godon et al., 1998. One protein, YMR108W, was listed as both induced and repressed in the original data set and is not included in this table.
    \item \code{yeastgfp.csv.xz} Has 28 columns; the names of the first five are \code{yORF}, \code{gene name}, \code{GFP tagged?}, \code{GFP visualized?}, and \code{abundance}. The remaining columns correspond to the 23 subcellular localizations considered in the YeastGFP project (Huh et al., 2003 and Ghaemmaghami et al., 2003) and hold values of either \code{T} or \code{F} for each protein. \samp{yeastgfp.csv} was downloaded on 2007-02-01 from http://yeastgfp.ucsf.edu using the Advanced Search, setting options to download the entire dataset and to include localization table and abundance, sorted by orf number. See \code{\link{yeastgfp}} for examples that use this file.
  }

  Files in \code{bison} contain BLAST results and taxonomic information for a metagenome:
  \itemize{
    \item \code{bisonN_vs_refseq55.blast.xz}, \code{bisonS_vs_refseq55.blast.xz}, \code{bisonR_vs_refseq55.blast.xz}, \code{bisonQ_vs_refseq55.blast.xz}, \code{bisonP_vs_refseq55.blast.xz} are partial tabular BLAST results for proteins in the Bison Pool Environmental Genome. Protein sequences predicted in the metagenome were downloaded from the Joint Genome Institute's IMG/M system on 2009-05-13. The target database for the searches was constructed from microbial protein sequences in National Center for Biotechnology Information (NCBI) RefSeq database version 55, representing 4567 microbial genomes. The \sQuote{blastall} command was used with the default setting for E value cuttoff (10.0) and options to make a tabular output file consisting of the top 20 hits for each query sequence. The function \code{\link{read.blast}} was used to extract only those hits with E values less than or equal to 1e-5 and with sequence similarity (percent identity) at least 30 percent, and to keep only the first hit for each query sequence. The function \code{\link{write.blast}} was used to save partial BLAST files (only selected columns). The files provided with CHNOSZ contain the first 5,000 hits for each sampling site at Bison Pool, representing between about 7 to 15 percent of the first BLAST hits after similarity and E value filtering.
    \item \code{gi.taxid.txt.xz} is a table that lists the sequence identifiers (gi numbers) that appear in the example BLAST files (see above), together with the corresponding taxon ids used in the NCBI databases. This file is \emph{not} a subset of the complete \sQuote{gi_taxid_prot.dmp.gz} available at \url{ftp://ftp.ncbi.nih.gov/pub/taxonomy/} but instead is a subset of \sQuote{gi.taxid.txt} generated from the RefSeq release catalog using \sQuote{gencat.sh} in the \code{refseq} directory. See \code{\link{id.blast}} for an example that uses this file and the BLAST files described above.
  }


  Files in \code{cpetc} contain heat capacity data and other thermodynamic properties:
  \itemize{
    \item \code{PM90.csv} Heat capacities of four unfolded aqueous proteins taken from Privalov and Makhatadze, 1990. Names of proteins are in the first column, temperature in \eqn{^{\circ}}{°}C in the second, and heat capacities in J mol\eqn{^{-1}}{^-1} K\eqn{^{-1}}{^-1} in the third. See \code{\link{ionize.aa}} for an example that uses this file.
    \item \code{RH95.csv} Heat capacity data for iron taken from Robie and Hemingway, 1995. Temperature in Kelvin is in the first column, heat capacity in J K\eqn{^{-1}}{^-1} mol\eqn{^{-1}}{^-1} in the second. See \code{\link{subcrt}} for an example that uses this file.
    \item \code{RT71.csv} pH titration measurements for unfolded lysozyme (\samp{LYSC_CHICK}) taken from Roxby and Tanford, 1971. pH is in the first column, net charge in the second. See \code{\link{ionize.aa}} for an example that uses this file.
    \item \code{SOJSH.csv} Experimental equilibrium constants for the reaction NaCl(aq) = Na+ + Cl- as a function of temperature and pressure taken from Fig. 1 of Shock et al., 1992. Data were extracted from the figure using g3data (\url{http://www.frantz.fi/software/g3data.php}). See \code{\link{water}} for an example that uses this file.
    \item \code{Cp.CH4.HW97.csv}, \code{V.CH4.HWM96.csv} Apparent molar heat capacities and volumes of CH4 in dilute aqueous solutions reported by Hnedkovsky and Wood, 1997 and Hnedkovsky et al., 1996. See \code{\link{EOSregress}} for examples that use these files.
    \item \code{BKM60_Fig7.dat} Eh-pH values for normal, wet and waterlogged soils from Fig. 7 of Baas Becking et al., 1960. See the \sQuote{anintro} vignette for an example that uses this file.
    \item \code{SC10_Rainbow.csv} Values of temperature (\eqn{^{\circ}}{°}C), pH and logarithms of activity of \eqn{\mathrm{CO_2}}{CO2}, \eqn{\mathrm{H_2}}{H2}, \eqn{\mathrm{NH_4^+}}{NH4+}, \eqn{\mathrm{H_2S}}{H2S} and \eqn{\mathrm{CH_4}}{CH4} for mixing of seawater and hydrothermal fluid at Rainbow field (Mid-Atlantic Ridge), taken from Shock and Canovas, 2010.
  }

  Files in \code{fasta} contain protein sequences:
  \itemize{
    \item \code{HTCC1062.faa.xz} is a FASTA file of 1354 protein sequences in the organism \emph{Pelagibacter ubique} HTCC1062 downloaded from the NCBI RefSeq collection on 2009-04-12. The search term was Protein: txid335992[Organism:noexp] AND "refseq"[Filter]. See \code{\link{util.fasta}} and \code{\link{revisit}} for examples that use this file.
    \item \code{EF-Tu.aln} consists of aligned sequences (394 amino acids) of elongation factor Tu (EF-Tu). The sequences correspond to those taken from UniProtKB for ECOLI (\emph{Escherichia coli}), THETH (\emph{Thermus thermophilus}) and THEMA (\emph{Thermotoga maritima}), and reconstructed ancestral sequences taken from Gaucher et al., 2003 (maximum likelihood bacterial stem and mesophilic bacterial stem, and alternative bacterial stem). See the \sQuote{formation} vignette for an example that uses this file.
  }

  Files in \code{protein} contain protein composition data for model organisms:
  \itemize{
    \item \code{Sce.csv.xz} Data frame of amino acid composition of proteins from the \emph{Saccharomyces} Genome Database. Contains twenty-two columns. Values in the first column are the rownumbers, the second column (\code{OLN}) has the ordered locus names of proteins, and the remaining twenty columns (\code{Ala}..\code{Val}) contain the numbers of the respective amino acids in each protein; the columns are arranged in alphabetical order based on the three-letter abbreviations for the amino acids. The source of data for \samp{Sce.csv} is the file \samp{protein_properties.tab} found on the FTP site of the SGD project on 2008-08-04. Blank entries were replaced with "NA" and column headings were added. See \code{\link{read.expr}} for examples that use this file.
    \item \code{Eco.csv.xz} Contains 24 columns. Values in the first column correspond to rownumbers, the second column {\code{AC}} holds the accession numbers of the proteins, the third column (\code{Name}) has the names of the corresponding genes, and the fourth column {\code{OLN}} lists the ordered locus names of the proteins. The remaining twenty columns (\code{A}..\code{Y}) give the numbers of the respective amino acids in each protein and are ordered alphabetically by the one-letter abbreviations of the amino acids. The sources of data for \samp{Eco.csv} are the files \samp{ECOLI.dat} \url{ftp://ftp.expasy.org/databases/hamap/complete_proteomes/entries/bacteria} and \samp{ECOLI.fas} \url{ftp://ftp.expasy.org/databases/hamap/complete_proteomes/fasta/bacteria} downloaded from the HAMAP (High-quality Automated and Manual Annotation of microbial Proteomes system) FTP site (Gattiker et al., 2003) on 2007-12-20. The proteins can be included in calculations using \code{\link{more.aa}} as well as \code{\link{read.expr}}.
%%    \item \code{HUM.csv.xz} Downloaded the file \code{uniprot_sprot_human.dat.gz}, dated 2010-08-10, from \url{ftp://ftp.expasy.org/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/}, converted from UniProt to FASTA format using the \code{seqret} tool from EMBOSS (\url{http://emboss.sourceforge.net/}). Counted amino acid frequencies using \code{\link{read.fasta}}. Columns are as described in \code{thermo$protein}, except column \code{protein} and \code{abbrv} contain Swiss-Prot name and accession number, respectively (both taken from the header lines in the FASTA file).
  }


  Files in \code{refseq} contain code and results of processing NCBI Reference Sequences (RefSeq) for microbial proteins, updated for RefSeq release 55 of 2012-09-17:
  \itemize{
    \item \code{README.txt} Instructions for producing the data files.
    \item \code{gencat.sh} Bash script to extract microbial protein records from the RefSeq catalog.
    \item \code{gi.taxid.txt} Output from above. The complete file is too large to distribute with CHNOSZ, but a portion is included in \code{extdata/bison} to support processing example BLAST files for the Bison Pool metagenome.
    \item \code{mkfaa.sh} Combine the contents of .faa.gz files into a single FASTA file (to use e.g. for making a BLAST database).
    \item \code{protein.refseq.R} Calculate average amino acid composition of all proteins for each organism identified by a taxonomic ID.
    \item \code{protein_refseq.csv.xz} Output from above. See example for \code{\link{ZC}}.
    \item \code{taxid.names.R} Generate a table of scientific names for the provided taxids. Requires the complete \code{names.dmp} and \code{nodes.dmp} from NCBI taxonomy files.
    \item \code{taxid_names.csv.xz} Output from above. See example for \code{\link{id.blast}}.
  }

  Files in \code{taxonomy} contain example taxonomic data files:
  \itemize{
    \item \code{names.dmp} and \code{nodes.dmp} are excerpts of the taxonomy files available on the NCBI ftp site (\url{ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz}, accessed 2010-02-15). These example files contain only the entries for \emph{Escherichia coli} K-12, \emph{Saccharomyces cerevisiae}, \emph{Homo sapiens}, \emph{Pyrococcus furisosus} and \emph{Methanocaldococcus jannaschii} (taxids 83333, 4932, 9606, 186497, 243232) and the higher-ranking nodes (genus, family, etc.) in the respective lineages. See \code{\link{taxonomy}} for examples that use this file.
  }

  Files in \code{thermo} contain additional thermodynamic data and group additivity definitions:
  \itemize{
    \item \code{OBIGT-2.csv} contains supplementary thermodynamic data in the same format as the primary database in \code{data/OBIGT.csv}. Data for some entries in the primary database are taken from different literature sources in this file. The default action of \code{\link{add.obigt}} is to add the contents of this file to CHNOSZ's working database in \code{thermo$obigt}. See \code{\link{diagram}} and the code of \code{\link{anim.TCA}} for examples that use this file.
    \item \code{obigt_check.csv} contains the results of running \code{\link{check.obigt}} to check the internal consistency of entries in the primary and supplementary databases.
    \item \code{groups_big.csv} Group contribution matrix: five structural groups on the columns ([-CH3],[-CH2-],[-CH2OH],[-CO-],[-COOH]) and 24 compounds on the rows (alkanes, alcohols, ketones, acids, multiply substituted compounds).
    \item \code{groups_small.csv} Group contribution matrix: twelve bond-specific groups on the columns, and 25 compounds on the rows (as above, plus isocitrate). Group identity and naming conventions adapted from Benson and Buss (1958) and Domalski and Hearing (1993). See the \sQuote{xadditivity} vignette for examples that use this file and \code{groups_big.csv}.
    \item \code{RH98_Table15.csv} Group stoichiometries for high molecular weight crystalline and liquid organic compounds taken from Table 15 of Richard and Helgeson, 1998. The first three columns have the \code{compound} name, \code{formula} and physical \code{state} (\samp{cr} or \samp{liq}). The remaining columns have the numbers of each group in the compound; the names of the groups (columns) correspond to species in \code{\link{thermo}$obigt}. The compound named \samp{5a(H),14a(H)-cholestane} in the paper has been changed to \samp{5a(H),14b(H)-cholestane} here to match the group stoichiometry given in the table. See \code{\link{RH2obigt}} for a function that uses this file.
    \item \code{DLEN67.csv} Standard Gibbs energies of formation, in kcal/mol, from Dayhoff et al., 1967, for nitrogen (N2) plus 17 compounds shown in Fig. 2 of Dayhoff et al., 1964, at 300, 500, 700 and 1000 K.
  }


}

\references{

  Anderson, N. L. and Anderson, N. G. (2002) The human plasma proteome: History, character and diagnostic prospects. \emph{Molecular and Cellular Proteomics} \bold{1}, 845--867. \url{http://dx.doi.org/10.1074/mcp.R200007-MCP200}

  Anderson, N. L. and Anderson, N. G. (2003) The human plasma proteome: History, character and diagnostic prospects (Vol. 1 (2002) 845-867). \emph{Molecular and Cellular Proteomics} \bold{2}, 50. \url{http://dx.doi.org/10.1074/mcp.A300001-MCP200}

  Baas Becking, L. G. M., Kaplan, I. R. and Moore, D. (1960) Limits of the natural environment in terms of pH and oxidation-reduction potentials. \emph{Journal of Geology} \bold{68}(3), 243--284. \url{http://www.jstor.org/stable/30059218} 

  Benson, S. W. and Buss, J. H. (1958) Additivity rules for the estimation of molecular properties. Thermodynamic properties. \emph{J. Chem. Phys.} \bold{29}, 546--572. \url{http://dx.doi.org/10.1063/1.1744539}

  Dayhoff, M. O. and Lippincott, E. R. and Eck, R. V. (1964) Thermodynamic Equilibria In Prebiological Atmospheres. \emph{Science} \bold{146}, 1461--1464. \url{http://dx.doi.org/10.1126/science.146.3650.1461}

  Dayhoff, M. O. and Lippincott, E. R., Eck, R. V. and Nagarajan (1967) Thermodynamic Equilibrium In Prebiological Atmospheres of C, H, O, N, P, S, and Cl. Report SP-3040, National Aeronautics and Space Administration. \url{http://ntrs.nasa.gov/search.jsp?R=19670017966}

  Domalski, E. S. and Hearing, E. D. (1993) Estimation of the thermodynamic properties of C-H-N-O-S-Halogen compounds at 298.15 K \emph{J. Phys. Chem. Ref. Data} \bold{22}, 805--1159. \url{http://dx.doi.org/10.1063/1.555927}

  Gattiker, A., Michoud, K., Rivoire, C., Auchincloss, A. H., Coudert, E., Lima, T., Kersey, P., Pagni, M., Sigrist, C. J. A., Lachaize, C., Veuthey, A.-L., Gasteiger, E. and Bairoch, A. (2003) Automatic annotation of microbial proteomes in Swiss-Prot. \emph{Comput. Biol. Chem.} \bold{27}, 49--58. \url{http://dx.doi.org/10.1016/S1476-9271(02)00094-4}

  Gaucher, E. A., Thomson, J. M., Burgan, M. F. and Benner, S. A (2003) Inferring the palaeoenvironment of ancient bacteria on the basis of resurrected proteins. \emph{Nature} \bold{425}(6955), 285--288. \url{http://dx.doi.org/10.1038/nature01977}

  Ghaemmaghami, S., Huh, W., Bower, K., Howson, R. W., Belle, A., Dephoure, N., O'Shea, E. K. and Weissman, J. S. (2003) Global analysis of protein expression in yeast. \emph{Nature} \bold{425}(6959), 737--741. \url{http://dx.doi.org/10.1038/nature02046}

  Huh, W. K., Falvo, J. V., Gerke, L. C., Carroll, A. S., Howson, R. W., Weissman, J. S. and O'Shea, E. K. (2003) Global analysis of protein localization in budding yeast. \emph{Nature} \bold{425}(6959), 686--691. \url{http://dx.doi.org/10.1038/nature02026}

  HAMAP system. HAMAP FTP directory, \url{ftp://ftp.expasy.org/databases/hamap/}

  Hnedkovsky, L., Wood, R. H. and Majer, V. (1996) Volumes of aqueous solutions of CH4, CO2, H2S, and NH3 at temperatures from 298.15 K to 705 K and pressures to 35 MPa. \emph{J. Chem. Thermodyn.} \bold{28}, 125--142. \url{http://dx.doi.org/10.1006/jcht.1996.0011}

  Hnedkovsky, L. and Wood, R. H. (1997) Apparent molar heat capacities of aqueous solutions of CH4, CO2, H2S, and NH3 at temperatures from 304 K to 704 K at a pressure of 28 MPa. \emph{J. Chem. Thermodyn.} \bold{29}, 731--747. \url{http://dx.doi.org/10.1006/jcht.1997.0192}

  Ishihama, Y., Schmidt, T., Rappsilber, J., Mann, M., Hartl, F. U., Kerner, M. J. and Frishman, D. (2008) Protein abundance profiling of the \emph{Escherichia coli} cytosol. \emph{BMC Genomics} \bold{9}:102. \url{http://dx.doi.org/10.1186/1471-2164-9-102}

  Joint Genome Institute (2007) Bison Pool Environmental Genome. Protein sequence files downloaded from IMG/M (\url{http://img.jgi.doe.gov/cgi-bin/m/main.cgi?section=FindGenomes&page=findGenomes})

  Privalov, P. L. and Makhatadze, G. I. (1990) Heat capacity of proteins. II. Partial molar heat capacity of the unfolded polypeptide chain of proteins: Protein unfolding effects. \emph{J. Mol. Biol.} \bold{213}, 385--391. \url{http://dx.doi.org/10.1016/S0022-2836(05)80198-6}

  Richard, L. and Helgeson, H. C. (1998) Calculation of the thermodynamic properties at elevated temperatures and pressures of saturated and aromatic high molecular weight solid and liquid hydrocarbons in kerogen, bitumen, petroleum, and other organic matter of biogeochemical interest. \emph{Geochim. Cosmochim. Acta} \bold{62}, 3591--3636. \url{http://dx.doi.org/10.1016/S0016-7037(97)00345-1}

  Robie, R. A. and Hemingway, B. S. (1995) \emph{Thermodynamic Properties of Minerals and Related Substances at 298.15 K and 1 Bar (\eqn{10^5} Pascals) Pressure and at Higher Temperatures}. U. S. Geol. Surv., Bull. 2131, 461 p. \url{http://www.worldcat.org/oclc/32590140}

  Roxby, R. and Tanford, C. (1971) Hydrogen ion titration curve of lysozyme in 6 M guanidine hydrochloride. \emph{Biochemistry} \bold{10}, 3348--3352. \url{http://dx.doi.org/10.1021/bi00794a005}

  SGD project. \emph{Saccharomyces} Genome Database, \url{http://www.yeastgenome.org}

  Shock, E. L., Oelkers, E. H., Johnson, J. W., Sverjensky, D. A. and Helgeson, H. C. (1992) Calculation of the thermodynamic properties of aqueous species at high pressures and temperatures: Effective electrostatic radii, dissociation constants and standard partial molal properties to 1000 \eqn{^{\circ}}{°}C and 5 kbar. \emph{J. Chem. Soc. Faraday Trans.} \bold{88}, 803--826. \url{http://dx.doi.org/10.1039/FT9928800803}

  Shock, E. and Canovas, P. (2010) The potential for abiotic organic synthesis and biosynthesis at seafloor hydrothermal systems. \emph{Geofluids} \bold{10}, 161--192. \url{http://dx.doi.org/10.1111/j.1468-8123.2010.00277.x}

  YeastGFP project. Yeast GFP Fusion Localization Database, http://yeastgfp.ucsf.edu; Current location: \url{http://yeastgfp.yeastgenome.org}

}

\keyword{datasets}
