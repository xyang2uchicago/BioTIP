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
## Loading gencode_gr dataset internally in R from data directory:
## load(file="data/gencode_gr.v19_chr21.rda")
##
devtools::use_data(gencode_gr, overwrite = TRUE)