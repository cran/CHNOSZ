# trim the protein_refseq.csv, removing some entries with highly-represented names
# (to keep the file size down for CHNOSZ package)
# 20130320 jmd

# the original file (7415 rows)
pr <- read.csv("protein_refseq.csv")
# the common names (comments show number in RefSeq 57
names <- c(
  "Escherichia coli",  # 662
  "Streptococcus",     # 432
  "Salmonella",        # 299
  "Staphylococcus",    # 290
  "Enterococcus",      # 218
  "Vibrio",            # 190
  "Lactobacillus",     # 179
  "Helicobacter",      # 164
  "Pseudomonas",       # 155
  "Mycobacterium",     # 138
  "Campylobacter",     # 131
  "Neisseria",         # 121
  "Clostridium",       # 118
  "Yersinia",          # 111
  "Bacillus cereus",   # 105
  "Acinetobacter",     # 103
  "Propionibacterium", # 93
  "Burkholderia",      # 91
  "Candidatus",        # 83
  "Bacteroides",       # 80
  "Mycoplasma",        # 73
  "Streptomyces",      # 72
  "Corynebacterium",   # 70
  "Listeria",          # 60
  "Leptospira",        # 57
  "Klebsiella",        # 57
  "Bifidobacterium",   # 55
  "Brucella",          # 53
  "Shigella",          # 50
  "Haemophilus",       # 47
  "Rickettsia",        # 46
  "Prevotella",        # 44
  "Chlamydia",         # 42
  "Francisella",       # 37
  "Bacillus thuringiensis", # 36
  "Borrelia",          # 34
  "Fusobacterium",     # 33
  "Xanthomonas",       # 31
  "Rhizobium",         # 27
  "Bartonella",        # 26
  "Pseudoalteromonas", # 24
  "Bacillus anthracis", # 23
  "Actinomyces",       # 23
  "Treponema",         # 21
  "Actinobacillus",    # 20
  "Aggregatibacter",   # 20
  "Gardnerella"        # 19
)

# loop over the names, identify rows, leave the first row
for(name in names) {
  iname <- grep(name, pr$ref)
  iname <- tail(iname, -1)
  pr <- pr[-iname, ]
}

# we're left with 2600 rows 
write.csv(pr, "protein_refseq.csv", row.names=FALSE)
