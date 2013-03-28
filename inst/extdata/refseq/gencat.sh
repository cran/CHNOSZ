#/bin/sh
# extract microbial, genomic records from the RefSeq catalog
RELEASE=57
ORG=microbial
MOL=protein
BASENAME=RefSeq-release$RELEASE.catalog 

# extract the microbial records
grep \|$ORG $BASENAME  > $BASENAME.$ORG

# extract the protein records
# alternatively, could use egrep:
#egrep "[[:space:]]AP_ | [[:space:]]NP_ | [[:space:]]XP_ | \
#  [[:space:]]YP_ | [[:space:]]ZP_" $BASENAME.$ORG  > $BASENAME.$ORG.$MOL
grep "[[:space:]]AP_" $BASENAME.$ORG  > $BASENAME.$ORG.$MOL  # 0 records
grep "[[:space:]]NP_" $BASENAME.$ORG >> $BASENAME.$ORG.$MOL  # 450218 records
grep "[[:space:]]XP_" $BASENAME.$ORG >> $BASENAME.$ORG.$MOL  # 0 records
grep "[[:space:]]YP_" $BASENAME.$ORG >> $BASENAME.$ORG.$MOL  # 6786156 records
grep "[[:space:]]ZP_" $BASENAME.$ORG >> $BASENAME.$ORG.$MOL  # 17252153 records

# to save only the gi, taxid and sequence length columns, in that order
# the field separator (tab) is defined in command line, not in awk program,
#   otherwise the first line gets processed incorrectly
cat $BASENAME.$ORG.$MOL | awk -F\\t '{print $4,$1,$7}' > gi.taxid.unsrt

# sort the file on gi so that it can be used with e.g. the unix 'join' command
cat gi.taxid.unsrt | sort > gi.taxid.txt
