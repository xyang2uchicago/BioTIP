##
## cod dataset
##
## Loading cod_gr dataset internally in R from data directory:
## load(file="data/cod.rda") or
## data(cod)
##
usethis::use_data("cod")

##
## Downloaded human GRCh37 gtf file from:
## \url{{https://www.gencodegenes.org/human/release_25lift37.html}
##
## Extracting gtf from gtf.gz:
##   use unzip, 7zip or gunzip tools.
##
## For extracting summary data from the downloaded gtf file:
## use the python code provided in "data-raw/python_extraction_code.R" directory.
##
## Loading gencode dataset internally in R from data directory:
## load(file="data/gencode.rda") or
## data(gencode)
##
usethis::use_data("gencode")

##
## ILEF dataset
##
## Loading ILEF dataset internally in R from data directory:
## load(file="data/ILEF.rda")
## data(ILEF)
##
usethis::use_data("ILEF")

##
## intron dataset
##
## Loading intron dataset internally in R from data directory:
## load(file="data/intron.rda")
## data(intron)
##
usethis::use_data("intron")

##
## GSE6136 matrix dataset
##
## Loading GSE6136 matrix dataset internally in R from data directory:
## load(file="data/GSE6136_matrix.rda") or
## data(GSE6136_matrix)
##
usethis::use_data("GSE6136_matrix")

##
## Summary of GSE6136 dataset = GSE6136_cli
##
## Loading GSE6136 cli dataset internally in R from data directory:
## load(file="data/GSE6136_cli.rda") or
## data(GSE6236_cli)
##
usethis::use_data("GSE6136_cli")
