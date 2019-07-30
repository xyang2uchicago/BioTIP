## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ---- echo=FALSE, fig.cap="Fig 1. NPS workflow with five key analytic steps", fig.align='center', out.width = '60%'----
knitr::include_graphics("Fig1.jpg")

## ---- echo=FALSE, fig.align='center', out.width = '65%'------------------
knitr::include_graphics("Fig2.jpg")

## ------------------------------------------------------------------------
#library("GenomeInfoDb")
#library("GenomicRanges")
#library("IRanges")
#library("NPS")
#data("gencode")
#head(gencode)

## ----"python code", eval = FALSE-----------------------------------------
#  gtf = ("PATH_FILE")
#  outF = open("gtf_summary_transbiotype.txt","w")
#  
#  def getquote(str,f,target):
#      targetLen = len(target)+2
#      strInd = str.find(target)
#      st = strInd + len(target)+2
#      ed = st + str[st:].find("";")
#      #print(st,ed)
#      f.write("\t"+str[st:ed]) if strInd!= -1 else f.write("\t"+"NA.")
#  
#  with open(gtf, "r") as f:
#       for line in f:
#          if line[0] != "#":
#              chromosome = line.split("\t")[0]
#              st = line.split("\t")[3]
#              ed = line.split("\t")[4]
#              strand = line.split("\t")[6]
#              type = line.split("\t")[2]
#              outF.write(chromosome+"\t"+st+"\t"+ed+"\t"+strand+"\t"+type)
#              c = "transcript_id"
#              g = "gene_name"
#              t = "transcript_type"
#              getquote(line,outF,c)
#              getquote(line,outF,g)
#              getquote(line,outF,t)
#              outF.write("\n")
#  outF.close()

## ------------------------------------------------------------------------
#library("GenomeInfoDb")
library("GenomicRanges")
#library("IRanges")
#library("NPS")
data("gencode")
#head(gencode)

## ------------------------------------------------------------------------
query <- GRanges(c("chr1:2-10:+","chr1:6-10:-"), Row.names = c("trans1","trans2"),score = c(1,2))
head(query)


## ------------------------------------------------------------------------
#library(NPS)
#library("GenomeInfoDb")
#library("GenomicRanges")
#library("IRanges")
gr <- GRanges(c("chr1:1-5:+","chr1:2-3:+"),biotype = c("lincRNA","CPC"))
head(gr)

## ------------------------------------------------------------------------
#library(NPS)
#library("GenomeInfoDb")
#library("GenomicRanges")
#library("IRanges")

intron <- GRanges("chr1:6-8:+")
head(intron)

## ------------------------------------------------------------------------
#library("NPS")
#library("GenomeInfoDb")
#library("GenomicRanges")
#library("IRanges")

#intron_trncp <- getBiotypes(query, gr, intron)
#intron_trncp

## ------------------------------------------------------------------------
#library("NPS")
library("GenomeInfoDb")
library("GenomicRanges")
library("IRanges")
#data("intron")
#data("ILEF")
#data("gencode")
#gencode_gr = GRanges(gencode)
#ILEF_gr = GRanges(ILEF)
#cod_gr = GRanges(cod)
#intron_gr = GRanges(intron)

#non-coding <- getBiotypes(ILEF_gr, gencode_gr, intron_gr)
#non-coding
#
#coding <- getBiotypes(ILEF_gr, gencode_gr)
#coding

## ------------------------------------------------------------------------
#library("NPS")
#library("GenomicRanges")
#library("IRanges")
#library("GenomeInfoDb")

data("ILEF")
data("cod")
#ILEF_gr = GRanges(ILEF)
#cod_gr = GRanges(cod)

#rdthrough <- getReadthrough(ILEF_gr, cod_gr)
#rdthrough

## ------------------------------------------------------------------------
#install.packages("stringr")
library("stringr")
library("psych")
library("igraph")
library("cluster")

## ------------------------------------------------------------------------
source("C:/Users/benjk/Desktop/NPS/R/NPS.R")
GSE6136 = read.table("C:/Users/benjk/Desktop/GSE6136_matrix.txt", header = TRUE, comment = "!")

dim(GSE6136)               #[1] 22690rows and 27 columns
row.names(GSE6136) = GSE6136$ID_REF
GSE6136 = GSE6136[,-1]
dim(GSE6136)               #[1] 22690 rows and 26 columns

## ------------------------------------------------------------------------
cli = read.delim("C:/Users/benjk/Desktop/GSE6136_cli.txt", head = F)
dim(cli)
cli = t(cli)
colnames(cli) = str_split_fixed(cli[1,],'_',2)[,2]
cli = cli[-1,]
cli = data.frame(cli)
cli[,"cell-type:ch1"] = str_split_fixed(cli$characteristics_ch1.1,": ",2)[,2]
cli[,"Ig clonality:ch1"] = str_split_fixed(cli$characteristics_ch1.3,": ",2)[,2]

colnames(cli)[colnames(cli) == "cell-type:ch1"] = "group"
cli$Row.names = cli[,1]
head(cli)

## ------------------------------------------------------------------------
dat <- GSE6136
df <- log2(dat+1)
head(df)

## ------------------------------------------------------------------------
tmp <- names(table(cli$group))
samplesL <- split(cli[,1],f = cli$group)
head(samplesL)
test <- sd_selection(df, samplesL,0.01)
head(test[["activated"]])

## ------------------------------------------------------------------------
igraphL <- getNetwork(test, fdr = 0.05)
#summary(igraphL)
clstr <- getCluster_methods(igraphL)
#summary(clstr)
head(clstr)

## ------------------------------------------------------------------------
#par(mfrow <- c(1,5))
#membersL_noweight <- getMCI(cluster, test, adjust.size = F, ylim = c(0,8))
#head(getMCI)
#maxCIms <- getMaxCImember(membersL_noweight[[1]],membersL_noweight[[2]],min =3)
#head(maxCIms)

## ----SessionInfo---------------------------------------------------------
sessionInfo()

