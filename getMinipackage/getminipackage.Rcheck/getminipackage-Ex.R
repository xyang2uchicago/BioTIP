pkgname <- "getminipackage"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
base::assign(".ExTimings", "getminipackage-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('getminipackage')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("getBiotypes")
### * getBiotypes

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: getBiotypes
### Title: Classifying Transcripts Using getminipackage
### Aliases: getBiotypes getReadthrough

### ** Examples

#Input datasets locally
#load(file="data/gencode_gr.v19_chr21.rda")
#load(file="data/intron_gr.chr21.rda")
#load("data/ILEF_gr.chr21.rda")
#example 1 using getBiotypes
getBiotypes(ILEF_gr,gencode_gr)

#### example 2 using getReadthrough ####
#load(file="data/gencode_gr.v19_chr21.rda")
#load(file="data/ILEF_gr.chr21.rda")
cod_gr <- subset(gencode_gr, biotype == 'protein_coding')
getReadthrough(ILEF_gr,cod_gr)

## Not run: getBiotypes("intron_gr")




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("getBiotypes", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
