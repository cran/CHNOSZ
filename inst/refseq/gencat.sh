#/bin/sh
# extract microbial, genomic records from the RefSeq catalog
RELEASE=45
ORG=microbial
MOL=protein
BASENAME=RefSeq-release$RELEASE.catalog 

# extract the microbial records
grep \|$ORG $BASENAME  > $BASENAME.$ORG

# extract the protein records
# alternatively, could use egrep:
#egrep "[[:space:]]AP_ | [[:space:]]NP_ | [[:space:]]XP_ | \
#  [[:space:]]YP_ | [[:space:]]ZP_" $BASENAME.$ORG  > $BASENAME.$ORG.$MOL
grep "[[:space:]]AP_" $BASENAME.$ORG  > $BASENAME.$ORG.$MOL
grep "[[:space:]]NP_" $BASENAME.$ORG >> $BASENAME.$ORG.$MOL
grep "[[:space:]]XP_" $BASENAME.$ORG >> $BASENAME.$ORG.$MOL
grep "[[:space:]]YP_" $BASENAME.$ORG >> $BASENAME.$ORG.$MOL
grep "[[:space:]]ZP_" $BASENAME.$ORG >> $BASENAME.$ORG.$MOL


# to save only the gi, taxid and sequence length columns
cat $BASENAME.$ORG.$MOL | awk '{FS="\t"} {print $4,$1,$7}' > gi.taxid.unsrt

# for some reason the first line needs to be corrected manually
# str. 316407 W3110 --> 89106885 316407 21

# sort the file on gi so that it can be used with e.g. the unix 'join' command
cat gi.taxid.unsrt | sort > gi.taxid.txt
