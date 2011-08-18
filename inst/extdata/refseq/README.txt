# the following data files support calculations using the 
# RefSeq database (release 47, 2011-05-07)
protein_refseq.csv.gz: overall (average) amino acid composition of all proteins for each
  microbial genome in the RefSeq collection (n=3266); generated using the steps below
taxon_names.csv: taxid, phylum name and species name for 3266 microbial taxa;
  generated from the taxids using sciname() and parent()

# these functions/scripts have the following purpose (output files listed in parentheses):
mkfaa.sh - combine gzipped sequence files into one big FASTA file (refseq47.faa)
gencat.sh - extract gi number, taxid, sequence length from RefSeq release catalog (gi.taxid.txt)
protein.refseq.R - get average amino acid composition for each taxid from gzipped sequence files (protein_refseq.csv)

# bash scripts assume a GNU/Linux-like operating system
# timings were made for processing RefSeq 47 on a recent (2009) intel laptop

# get the list of files and entries in the database
1. download 'release47.files.installed' and 'RefSeq-release47.catalog' from NCBI
   (ftp://ftp.ncbi.nih.gov/refseq/release/release-catalog)

# download stuff
3. list URLS for the microbial protein sequence files:
     grep microbial.*.protein.faa* release47.files.installed | \
       sed -e "s/^/ftp\:\/\/ftp.ncbi.nih.gov\/refseq\/release\/microbial\//g" > urllist
4. download the files using 'wget -i urllist' [1.9 GB, ~2 hours]
5. move the .gz files to a directory named 'protein'
6. run ls protein/*.gz > filelist
7. use 'mkfaa.sh' to combine the sequences into a single file 
   (change the OUTFILE first if needed) [~2 minutes]

# protein stuff
8. use 'gencat.sh' to generate gi.taxid.txt (first change RELEASE if needed)
   note that the intermediate file gi.taxid.unsrt may have to be edited manually 
     -- see instructions in gencat.sh
   when done, the output of 'cat gi.taxid.txt | wc -l'  
   should be equal to 'grep "^>" refseq47.faa | wc -l'
   (for microbial proteins in RefSeq 47, the number is 9920861)
9. generate protein_refseq.csv in R:  [~5.5 hours]
   > source("protein.refseq.R")
   > protein.refseq()
   note that this depends on gi.taxid.txt and the .faa.gz files in the 'protein' directory

# taxonomy stuff
10. edit 'microbial.taxa.R' so that 'taxdir' points to the directory where the files
    'names.dmp' and 'nodes.dmp' are present. these files can be downloaded from
    ftp://ftp.ncbi.nih.gov/pub/taxonomy/ .
11. source 'taxon.names.R' to generate the file 'taxon_names.csv' [~1.5 hours]

# BLAST stuff (optional)
12. to create a BLAST database use a command like
    formatdb -t refseq47 -i refseq47.faa -p T -o T  [~15 minutes]
    (tested with BLAST version 2.2.24)
