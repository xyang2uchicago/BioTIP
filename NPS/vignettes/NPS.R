## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ------------------------------------------------------------------------
#library("DNS")
#library("GenomeInfoDb")
library("GenomicRanges")
library("IRanges")
data("gencode")
gencode

## ----"python code", eval = FALSE-----------------------------------------
#  gtf = ('PATH_FILE')
#  outF = open('gtf_summary_transbiotype.txt','w')
#  
#  def getquote(str,f,target):
#      targetLen = len(target)+2
#      strInd = str.find(target)
#      st = strInd + len(target)+2
#      ed = st + str[st:].find('";')
#      #print(st,ed)
#      f.write('\t'+str[st:ed]) if strInd!= -1 else f.write('\t'+'NA.')
#  
#  with open(gtf, 'r') as f:
#       for line in f:
#          if line[0] != '#':
#              chromosome = line.split('\t')[0]
#              st = line.split('\t')[3]
#              ed = line.split('\t')[4]
#              strand = line.split('\t')[6]
#              type = line.split('\t')[2]
#              outF.write(chromosome+'\t'+st+'\t'+ed+'\t'+strand+'\t'+type)
#              c = 'transcript_id'
#              g = 'gene_name'
#              t = 'transcript_type'
#              getquote(line,outF,c)
#              getquote(line,outF,g)
#              getquote(line,outF,t)
#              outF.write('\n')
#  outF.close()

## ----quickstart----------------------------------------------------------
#library("GenomeInfoDb")
library("GenomicRanges")
library("IRanges")
data("gencode")
head(gencode)

## ------------------------------------------------------------------------
query <- GRanges(c('chr1:2-10:+','chr1:6-10:-'),Row.names = c('trans1','trans2'),score = c(1,2))
head(query)


## ------------------------------------------------------------------------
#library(DNS)
#library("GenomeInfoDb")
library("GenomicRanges")
library("IRanges")
gr <- GRanges(c('chr1:1-5:+','chr1:2-3:+'),biotype = c('lincRNA','CPC'))
head(gr)

## ------------------------------------------------------------------------
#library(DNS)
#library("GenomeInfoDb")
library("GenomicRanges")
library("IRanges")

intron <- GRanges('chr1:6-8:+')
head(intron)

## ------------------------------------------------------------------------
library(DNS)
#library("GenomeInfoDb")
library("GenomicRanges")
library("IRanges")

coding_trncp <- getBiotypes(query, gr, intron)
head(coding_trncp)

## ------------------------------------------------------------------------
library("DNS")
#library("GenomeInfoDb")
library("GenomicRanges")
library("IRanges")
data("intron")
data("ILEF")
data("gencode")
gencode_gr = GRanges(gencode)
ILEF_gr = GRanges(ILEF)
cod_gr = GRanges(cod)
intron_gr = GRanges(intron)

biotyp <- getBiotypes(ILEF_gr, gencode_gr, intron_gr)
head(biotyp)

## ------------------------------------------------------------------------
#library("DNS")
library("GenomicRanges")
library("IRanges")
#library("GenomeInfoDb")

data("ILEF")
data("cod")
ILEF_gr = GRanges(ILEF)
cod_gr = GRanges(cod)

rdthrough <- getReadthrough(ILEF_gr, cod_gr)
head(rdthrough)

## ----SessionInfo---------------------------------------------------------
sessionInfo()

