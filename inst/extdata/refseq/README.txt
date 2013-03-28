# the following data files support calculations using the 
# RefSeq database (release 57, 2013-01-08)
protein_refseq.csv: overall (average) amino acid composition of all proteins for each
  microbial genome in the RefSeq collection (n=7415)
taxid_names.csv: taxid, phylum name and species name for 7415 microbial taxa

# these functions/scripts have the following purpose (output files listed in parentheses):
mkfaa.sh - combine gzipped sequence files into one big FASTA file (refseq57.faa)
gencat.sh - extract gi number, taxid, sequence length from RefSeq release catalog (gi.taxid.txt)
protein.refseq.R - get average amino acid composition for each taxid from gzipped sequence files (protein_refseq.csv)
taxid.names.R - get taxonomic names for each taxid represented (taxid_names.csv)

# bash scripts assume a GNU/Linux-like operating system
# timings were made for processing RefSeq 57 on a recent (2009) intel laptop

# download stuff
1. download 'release57.files.installed' and 'RefSeq-release57.catalog.gz' from NCBI
   (ftp://ftp.ncbi.nih.gov/refseq/release/release-catalog)
2. list URLS for the microbial protein sequence files:
     grep microbial.*.protein.faa* release57.files.installed | \
       sed -e "s/^/ftp\:\/\/ftp.ncbi.nih.gov\/refseq\/release\/microbial\//g" > urllist
3. download the files using 'wget -i urllist' [3227 files, 4.6 GB]
4. move the .gz files to a directory named 'protein'

# protein stuff
5. gzip -d RefSeq-release57.catalog.gz [3.1 GB]
6. use 'gencat.sh' to generate gi.taxid.txt for microbial proteins in the catalog [7 minutes]
   for RefSeq57, 'cat gi.taxid.txt | wc -l' is 24488527
7. generate protein_refseq.csv in R:  [~19 hours]
   > source("protein.refseq.R")
   > protein.refseq()
   note that this depends on gi.taxid.txt and the .faa.gz files in the 'protein' directory
8. trim entries from protein_refseq.csv (smaller size, better for package distribution)
   > source("trim_refseq.R")

# taxonomy stuff
9. edit 'taxid.names.R' so that 'taxdir' points to the directory where the files
    'names.dmp' and 'nodes.dmp' are present. these files can be downloaded from
    ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz (accessed on 2013-01-15)
10. source 'taxid.names.R' to generate the file 'taxid_names.csv' [~4.5 hours]

# BLAST stuff (optional)
11. run ls protein/*.gz > filelist
12. use 'mkfaa.sh' to combine the sequences into a single file 'refseq57.faa' [9.3 GB, 11 minutes]
    for RefSeq57, 'grep "^>" refseq57.faa | wc -l' is 24477649
    the difference from the catalog (step 6 above), is 10878 sequences:
    taxid 1211777 (Rhizobium mesoamericanum STM3625) (6356 sequences) 
    taxid 313627 (Bacillus sp. NRRL B-14911) (4522 sequences)
13. make a BLAST database, e.g. formatdb -t refseq57 -i refseq57.faa -p T

