## cross-checking sources
# the reference sources
ref.source <- thermo$refs$key
# sources in the primary thermodynamic database
# we omit the [S92] in "HDNB78 [S92]" etc.
os1 <- gsub("\ .*", "", thermo$obigt$ref1)
os2 <- gsub("\ .*", "", thermo$obigt$ref2)
# all of the thermodynamic data sources - some of them might be NA
obigt.source <- unique(c(os1,os2))
obigt.source <- obigt.source[!is.na(obigt.source)]
# these all produce character(0) if the sources are all accounted for
print("missing these sources for thermodynamic properties:")
print(unique(obigt.source[!(obigt.source %in% ref.source)]))
# determine if all the reference sources are cited
# this should produce character(0)
print("these sources are present but not cited:")
print(ref.source[!(ref.source %in% obigt.source)])
