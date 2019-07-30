---
title: "NPS -- Adopting tipping-point theory to transcriptome profiles to unravel disease regulatory trajectory"
author: "Zhezhen Wang, Biniam Feleke and Xinan H.Yang"
date: "07/26/2019"
abstract: >
  Network perturbation signature (NPS). The purpose of this R package NPS, is to
  analyze clinical transcriptome data in order to determine critical transition
  from phenotype-defined sample states. Fundamental of the tipping-point is
  applied to transcriptomic profiles where highly connected transcriptional
  configurations serve as attractors [1]. NPS improves on a method called
  Dynamic Network Biomarker (DNB) that was calculated within each state [2].
  When studying NB, expressional patterns are often the aggregates of multiple
  dynamic subsystems and thus a new method is needed to allow for cross-state
  comparisons. This makes it vague to use DNB for tipping-point detection
  because the autocorrelation of non-DNB-modules at other states could attribute
  to an even higher DNB-score at a predicted ‘tipping point’ (Figure 1). By
  incorporating an index of critical state transition (IC-score) [3] into the
  DNB model, NPS overcomes this pitfall by checking for  drastic deviations
  across samples around the tipping point. Thus, the NPS scheme acts as a hybrid
  model to join the advantages of both DNB and IC-scoring.
output:  
  html_document: 
    fig_caption: yes
    keep_md: yes
    toc: yes
  pdf_document: default
  word_document: default
vignette: >
  %\VignetteIndexEntry{title}
  %\VignetteEncoding{UTF-8}
  %\VignetteEngine{knitr::rmarkdown}
---

#### ####



<div class="figure" style="text-align: center">
<img src="Fig1.jpg" alt="Fig 1. NPS workflow with five key analytic steps" width="60%" />
<p class="caption">Fig 1. NPS workflow with five key analytic steps</p>
</div>

#### [Standard workflow](#Standard workflow)
   * ##### [Transcript Annotation and Biotype](#Transcript Annotation and Biotype)
     *  ##### [Quick Start](#Quick Start)
     *  ##### [Genomic Data Source](#Genomic Data Source)
     *  ##### [Extracting Summary Data](#Extracting Summary Data)
     *  ##### [Loading Data](#Loading Data)
     *  ##### [Filtering chr21](#Filtering chr21)
   * ##### [An Identification of Critical Tipping Point](#An Identification of Critical Tipping Point)
     *  ##### [Data Preprocessing](#Data Preprocessing)
     *  ##### [Pre-selection Transcript](#Pre-selection Transcript)
     *  ##### [Network Partition](#Network Partition)
     *  ##### [Identifying Dynamic Network Biomodule](#Identifying Dynamic Network Biomodule)
     *  ##### [Finding Tipping Point](#Finding Tipping Point)
   * ##### [Acknowledgements](#Acknowledgements)
   * ##### [SessionInfo](#SessionInfo)
   * ##### [References](#References)
<a name="Transcript Annotation and Biotype"></a>

#### __Transcript Annotation and Biotype__
<a name="Quick Start"></a>

#### __Quick Start__

  The R function getBiotype is used to group transcripts of interest into 11
  biotypes based on GENCODE annotation (Fig 2a). When a query transcript
  overlaps at least half with a GENCODE transcript on the same strand, this
  query transcript will inherit the biotype of the GENCODE transcript.

  In the previous study conducted, five out of the 11 biotypes showed high
  protein-coding potential while the others did not (Fig 2b) [4]. We thus
  concluded these seven other biotypes, including protein-coding antisense RNAs,
  to be lncRNAs. The remaining coding biotypes in this study included canonic
  protein coding (CPC), ‘PC_mixed’, and ‘PC_intron’.

  First start by loading the required libraries: “GenomeInfoDb,” “NPS,”
  “GenomicRanges,” “IRanges” and “NPS”. Next load the datasets: “gencode”,
  “ILEF”, “intron” and “cod”. Using these datasets, excute NPS functions
  getBiotypes and getReadthrough as follows. These steps assume you installed
  the “NPS” package. If you did not install the package, use the
  `install.packages("NPS")` to install in R.

<img src="Fig2.jpg" width="65%" style="display: block; margin: auto;" />
Fig 2. A getBiotypes workflow and protein-coding potential in real data analysis
[4]. (a) Workflow of an in-house R function (getBiotypes) to query transcripts
of interests and classify into biotypes. (b) Pie-chart of eleven types of
transcripts assembled from polyadenylated RNA(TARGET). (c) Empirical cumulative
distribution plot comparing the transcripts across all 11 biotypes. The
protein-coding potential was estimated with the Coding Potential Assessment Tool
(CPAT). Line color codes biotypes. The more a line towards the right-bottom
corner, the lower protein-coding potential it has.


```r
#library("GenomeInfoDb")
#library("GenomicRanges")
#library("IRanges")
#library("NPS")
#data("gencode")
#head(gencode)
```

  These illustrations above assumes you have installed "NPS" package. If you did
  not install the pacakge, use the `install.packages("NPS")` to install in R.

<a name="Genomic Data Source"></a>\

#### __Genomic Data Source__ 

  High quality human genome sequence data can be obtained from various sources.
  We obtained a comprehensive gene annotation of human GRCh37 from
  [GENCODE](https://www.gencodegenes.org/human). For our illustrations, human
  GRCh37 data will be used. A standard file structure, similar to gtf format, is
  required for this package. The general transfer format (gtf) organizes genomic
  data in rows and columns (fields). Each row contains information about
  specific samples. The columns are tab separated headers of the data
  frame.There are eight fixed columns with specific headers. An example of gtf
  format is shown below. For details of the gtf file format visit this
  [link](https://useast.ensembl.org/info/website/upload/gff.html#tracklines
  target="_blank").


   ![Chromosome 21 of human GRCh37 gtf](chr21.jpg)

  The table above shows head of chr21 dataset which was extracted from a full
  genome dataset GRCh37. The extraction of subset Chr21 is described below.

<a name="">[GoTop](#)</a> 
<a name="Extracting Summary Data"></a>\

##### __Extracting Summary Data__

  Before any further analysis, we need to summarize the content of the raw gtf
  data. There are two ways to get genome biotypes: a) "transctipt_type" b)
  "gene_type". Due to our interst in coding and noncoding regions, the
  `transcript_type` method was used to extract the regions of interest using
  python script shown below. __Note__ that the `"PATH_FILE"` refers to the path
  where the downloded gtf file is located. For instance, if the gtf file is
  located on your `desktop`, replace the `"PATH_FILE"` with
  `"Users/user/Desktop/gtf"`.

**Python codes:**


```r
gtf = ("PATH_FILE")
outF = open("gtf_summary_transbiotype.txt","w")

def getquote(str,f,target):
    targetLen = len(target)+2
    strInd = str.find(target)
    st = strInd + len(target)+2
    ed = st + str[st:].find("";")
    #print(st,ed)
    f.write("\t"+str[st:ed]) if strInd!= -1 else f.write("\t"+"NA.")

with open(gtf, "r") as f:
     for line in f:
        if line[0] != "#":
            chromosome = line.split("\t")[0]
            st = line.split("\t")[3]
            ed = line.split("\t")[4]
            strand = line.split("\t")[6]
            type = line.split("\t")[2]
            outF.write(chromosome+"\t"+st+"\t"+ed+"\t"+strand+"\t"+type)
            c = "transcript_id"
            g = "gene_name"
            t = "transcript_type"
            getquote(line,outF,c)
            getquote(line,outF,g)
            getquote(line,outF,t)
            outF.write("\n")
outF.close() 
```
***
<a name="Home">[GoTop](#)</a>
<a name="Loading Data"></a>\ 

#### __Loading Data__

  In order to load your data from a local drive, use the following format.
  __Note__ that the `"PATH_FILE"` refers to the location of the summary data
  from the above section. For more details on how to load datasets click
  [here](https://support.rstudio.com/hc/en-us/articles/218611977-Importing-Data-with-RStudio)

##### loading data from local drive
 > data <- read.delim("PATH_FILE", comment.char = "#")

  Internal NPS package data is included in the data folder. The data can be
  loaded into R working console using `data()` function. Here we show an example
  of how to load a dataset `gencode` from the data directory. A quick view of
  the data can be achieved using `head(gencode)`.


```r
#library("GenomeInfoDb")
library("GenomicRanges")
#library("IRanges")
#library("NPS")
data("gencode")
#head(gencode)
```
<a name="Home">[GoTop](#)</a>

<a name="Filtering chr21 gencode"></a>\ 

#### __Filtering chr21__

  Here we show an extraction of "gencode" dataset using R commands. Note
  to replace `PATH_FILE` with file direcotry path. `gtf` refers to the full 
  genome file. A subset function was used to filter chr21 dataset as follows.

`chr21 <- subset(gencode, seqnames == "chr21")` #"genecode" = whole genome gtf

    > gtf = read.table("PATH_FILE")
    > gtf = subset(gtf, biotype == "transcript")
    > colnames(gtf) = c("chr","start","end","strand","biotype")
    > gr = GRanges(gtf)

<a name="Home">[GoTop](#)</a> 
<a name="Examples"></a>\

#### __Examples__

##### Processing Query
***


```r
query <- GRanges(c("chr1:2-10:+","chr1:6-10:-"), Row.names = c("trans1","trans2"),score = c(1,2))
head(query)
```

```
## GRanges object with 2 ranges and 2 metadata columns:
##       seqnames    ranges strand |   Row.names     score
##          <Rle> <IRanges>  <Rle> | <character> <numeric>
##   [1]     chr1      2-10      + |      trans1         1
##   [2]     chr1      6-10      - |      trans2         2
##   -------
##   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```

##### Classifying Biotypes
***


```r
#library(NPS)
#library("GenomeInfoDb")
#library("GenomicRanges")
#library("IRanges")
gr <- GRanges(c("chr1:1-5:+","chr1:2-3:+"),biotype = c("lincRNA","CPC"))
head(gr)
```

```
## GRanges object with 2 ranges and 1 metadata column:
##       seqnames    ranges strand |     biotype
##          <Rle> <IRanges>  <Rle> | <character>
##   [1]     chr1       1-5      + |     lincRNA
##   [2]     chr1       2-3      + |         CPC
##   -------
##   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```

##### Extracting intron coordinates
*** 

      # Intron coordinates
      
       intron <- GRanges("chr1:6-8:+")
  

```r
#library(NPS)
#library("GenomeInfoDb")
#library("GenomicRanges")
#library("IRanges")

intron <- GRanges("chr1:6-8:+")
head(intron)
```

```
## GRanges object with 1 range and 0 metadata columns:
##       seqnames    ranges strand
##          <Rle> <IRanges>  <Rle>
##   [1]     chr1       6-8      +
##   -------
##   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```

##### Filtering coding transcripts
***

    # Filtering non-coding regions using products from example 1, 2 and 3


```r
#library("NPS")
#library("GenomeInfoDb")
#library("GenomicRanges")
#library("IRanges")

#intron_trncp <- getBiotypes(query, gr, intron)
#intron_trncp
```

    # Filtering Intron and Exons 

  Here we show how to obtain protien coding and non-coding from our datasets.
  The coding transcrptis are an expressed sections of the genome that is
  responsibe for protein fomation. Meanwhile the non-coding transcripts are
  vital in the fromation reguraltory elements such promoters, enhancers and
  silencers.


```r
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
```


##### Finding overlapping transcripts
***
    # Samples with overlapping coding regions.

```r
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
```

#### __An Identification of Critical Tipping Point__
#### __Data Preprocessing__
#### __Pre-selection Transcript__
#### __Network Partition__
#### __Identifying Dynamic Network Biomodule__
#### __Finding Tipping Point__

  Next we will filter clinical data and group the transcripts into five
  different stages of disease progression: `"activated"`,
  `"lymphoma_aggressive"`, `"lymphoma_marginal"`, `"lymphoma_transitional,"` and
  `"resting"` as disease stage names. These stages resemble the general disease
  progression patterns: normal(healthy) state, predisease state, tipping point,
  critical transition stage(CTS) and disease stage. Currently it is vey
  difficult to identify when a disease progress from predisease state to
  critical transition state. In most cases these transitions occurs very fast.
  Ultimately we want to identify this critical tipping point before it
  progresses into a early disease stage.

  An existing dataset, GSE6136, is used to demonestrate how our functions are
  applied. Start by loading packages "stringr", "psych", "igraph" and "cluster".
  If these packages are not installed, use the `install.packages("packages")`
  function to install them. Then load them using `library("package")`. Below are
  some examples.


```r
#install.packages("stringr")
library("stringr")
library("psych")
```

```
## 
## Attaching package: 'psych'
```

```
## The following object is masked from 'package:NPS':
## 
##     reflect
```

```
## The following object is masked from 'package:IRanges':
## 
##     reflect
```

```r
library("igraph")
```

```
## 
## Attaching package: 'igraph'
```

```
## The following objects are masked from 'package:NPS':
## 
##     decompose, spectrum, union
```

```
## The following object is masked from 'package:GenomicRanges':
## 
##     union
```

```
## The following object is masked from 'package:IRanges':
## 
##     union
```

```
## The following object is masked from 'package:S4Vectors':
## 
##     union
```

```
## The following objects are masked from 'package:BiocGenerics':
## 
##     normalize, path, union
```

```
## The following objects are masked from 'package:stats':
## 
##     decompose, spectrum
```

```
## The following object is masked from 'package:base':
## 
##     union
```

```r
library("cluster")
```

  Once all the required packages are installed, load using "read.table()"
  function as follows. Note to change the `read.table(PATH/TO/YOUR/FILE)` when
  running your datasets. To check the dimension of your data set use
  "dim(dataset)" function. Here we show a use of dim function "dim(GSE6136)".
  Notice that after editing the column and row, the dimension of our dataset
  changes from (22690,27) to (22690, 26) because we removed the downloaded first
  row after assigning it to be column-name of the final numeric data matrix.


```r
source("C:/Users/benjk/Desktop/NPS/R/NPS.R")
GSE6136 = read.table("C:/Users/benjk/Desktop/GSE6136_matrix.txt", header = TRUE, comment = "!")

dim(GSE6136)               #[1] 22690rows and 27 columns
```

```
## [1] 22690    27
```

```r
row.names(GSE6136) = GSE6136$ID_REF
GSE6136 = GSE6136[,-1]
dim(GSE6136)               #[1] 22690 rows and 26 columns
```

```
## [1] 22690    26
```

  The summary of GSE6136 dataset is shown below.


```r
cli = read.delim("C:/Users/benjk/Desktop/GSE6136_cli.txt", head = F)
dim(cli)
```

```
## [1] 36 27
```

```r
cli = t(cli)
colnames(cli) = str_split_fixed(cli[1,],'_',2)[,2]
cli = cli[-1,]
cli = data.frame(cli)
cli[,"cell-type:ch1"] = str_split_fixed(cli$characteristics_ch1.1,": ",2)[,2]
cli[,"Ig clonality:ch1"] = str_split_fixed(cli$characteristics_ch1.3,": ",2)[,2]

colnames(cli)[colnames(cli) == "cell-type:ch1"] = "group"
cli$Row.names = cli[,1]
head(cli)
```

```
##    geo_accession                status submission_date last_update_date
## V2     GSM142398 Public on Oct 28 2006     Oct 25 2006      Mar 18 2009
## V3     GSM142399 Public on Oct 28 2006     Oct 25 2006      Mar 18 2009
## V4     GSM142400 Public on Oct 28 2006     Oct 25 2006      Mar 18 2009
## V5     GSM142401 Public on Oct 28 2006     Oct 25 2006      Mar 18 2009
## V6     GSM142402 Public on Oct 28 2006     Oct 25 2006      Mar 18 2009
## V7     GSM142403 Public on Oct 28 2006     Oct 25 2006      Mar 18 2009
##    type channel_count                      source_name_ch1 organism_ch1
## V2  RNA             1             wildtype resting B-Cells Mus musculus
## V3  RNA             1           wildtype activated B-Cells Mus musculus
## V4  RNA             1            E-mu-BRD2 resting B-Cells Mus musculus
## V5  RNA             1            E-mu-BRD2 resting B-Cells Mus musculus
## V6  RNA             1 marginal E-mu-BRD2 mediated lymphoma Mus musculus
## V7  RNA             1 marginal E-mu-BRD2 mediated lymphoma Mus musculus
##    characteristics_ch1        characteristics_ch1.1  characteristics_ch1.2
## V2  genotype: wildtype           cell-type: resting     splenomegaly: none
## V3  genotype: wildtype         cell-type: activated     splenomegaly: none
## V4 genotype: E-mu-BRD2           cell-type: resting     splenomegaly: none
## V5 genotype: E-mu-BRD2           cell-type: resting     splenomegaly: none
## V6 genotype: E-mu-BRD2 cell-type: lymphoma_marginal splenomegaly: marginal
## V7 genotype: E-mu-BRD2 cell-type: lymphoma_marginal     splenomegaly: mild
##        characteristics_ch1.3 molecule_ch1
## V2  Ig clonality: polyclonal    total RNA
## V3  Ig clonality: polyclonal    total RNA
## V4  Ig clonality: polyclonal    total RNA
## V5  Ig clonality: polyclonal    total RNA
## V6 Ig clonality: oligoclonal    total RNA
## V7 Ig clonality: oligoclonal    total RNA
##                                                                                                            extract_protocol_ch1
## V2 Total RNA was isolated with the guanidinium thiocyanate method, extracted with acid phenol and precipitated from 2-propanol.
## V3 Total RNA was isolated with the guanidinium thiocyanate method, extracted with acid phenol and precipitated from 2-propanol.
## V4 Total RNA was isolated with the guanidinium thiocyanate method, extracted with acid phenol and precipitated from 2-propanol.
## V5 Total RNA was isolated with the guanidinium thiocyanate method, extracted with acid phenol and precipitated from 2-propanol.
## V6 Total RNA was isolated with the guanidinium thiocyanate method, extracted with acid phenol and precipitated from 2-propanol.
## V7 Total RNA was isolated with the guanidinium thiocyanate method, extracted with acid phenol and precipitated from 2-propanol.
##    label_ch1
## V2    biotin
## V3    biotin
## V4    biotin
## V5    biotin
## V6    biotin
## V7    biotin
##                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 label_protocol_ch1
## V2 Using a poly-dT primer incorporating a T7 promoter, double-stranded cDNA was synthesized from 10 micrograms total RNA using a cDNA synthesis kit (SuperScript double-stranded cDNA synthesis kit; Invitrogen, Carlsbad, CA).  Double-stranded cDNA was purified by phenol/chloroform extraction using a Phase-Lock Gel (PLG Heavy Brinkmann Instruments, Westbury, NY) followed by ethanol precipitation.  Biotin-labeled cRNA was generated from the double-stranded cDNA template though in-vitro transcription with T7 polymerase, and a nucleotide mix containing biotinylated CTP and UTP (Enzo RNA Transcript Labeling Kit; Enzo Diagnostics, Farmingdale, NY).  The biotinylated cRNA was purified using RNeasy affinity columns (Qiagen, Valencia, MD).
## V3 Using a poly-dT primer incorporating a T7 promoter, double-stranded cDNA was synthesized from 10 micrograms total RNA using a cDNA synthesis kit (SuperScript double-stranded cDNA synthesis kit; Invitrogen, Carlsbad, CA).  Double-stranded cDNA was purified by phenol/chloroform extraction using a Phase-Lock Gel (PLG Heavy Brinkmann Instruments, Westbury, NY) followed by ethanol precipitation.  Biotin-labeled cRNA was generated from the double-stranded cDNA template though in-vitro transcription with T7 polymerase, and a nucleotide mix containing biotinylated CTP and UTP (Enzo RNA Transcript Labeling Kit; Enzo Diagnostics, Farmingdale, NY).  The biotinylated cRNA was purified using RNeasy affinity columns (Qiagen, Valencia, MD).
## V4 Using a poly-dT primer incorporating a T7 promoter, double-stranded cDNA was synthesized from 10 micrograms total RNA using a cDNA synthesis kit (SuperScript double-stranded cDNA synthesis kit; Invitrogen, Carlsbad, CA).  Double-stranded cDNA was purified by phenol/chloroform extraction using a Phase-Lock Gel (PLG Heavy Brinkmann Instruments, Westbury, NY) followed by ethanol precipitation.  Biotin-labeled cRNA was generated from the double-stranded cDNA template though in-vitro transcription with T7 polymerase, and a nucleotide mix containing biotinylated CTP and UTP (Enzo RNA Transcript Labeling Kit; Enzo Diagnostics, Farmingdale, NY).  The biotinylated cRNA was purified using RNeasy affinity columns (Qiagen, Valencia, MD).
## V5 Using a poly-dT primer incorporating a T7 promoter, double-stranded cDNA was synthesized from 10 micrograms total RNA using a cDNA synthesis kit (SuperScript double-stranded cDNA synthesis kit; Invitrogen, Carlsbad, CA).  Double-stranded cDNA was purified by phenol/chloroform extraction using a Phase-Lock Gel (PLG Heavy Brinkmann Instruments, Westbury, NY) followed by ethanol precipitation.  Biotin-labeled cRNA was generated from the double-stranded cDNA template though in-vitro transcription with T7 polymerase, and a nucleotide mix containing biotinylated CTP and UTP (Enzo RNA Transcript Labeling Kit; Enzo Diagnostics, Farmingdale, NY).  The biotinylated cRNA was purified using RNeasy affinity columns (Qiagen, Valencia, MD).
## V6 Using a poly-dT primer incorporating a T7 promoter, double-stranded cDNA was synthesized from 10 micrograms total RNA using a cDNA synthesis kit (SuperScript double-stranded cDNA synthesis kit; Invitrogen, Carlsbad, CA).  Double-stranded cDNA was purified by phenol/chloroform extraction using a Phase-Lock Gel (PLG Heavy Brinkmann Instruments, Westbury, NY) followed by ethanol precipitation.  Biotin-labeled cRNA was generated from the double-stranded cDNA template though in-vitro transcription with T7 polymerase, and a nucleotide mix containing biotinylated CTP and UTP (Enzo RNA Transcript Labeling Kit; Enzo Diagnostics, Farmingdale, NY).  The biotinylated cRNA was purified using RNeasy affinity columns (Qiagen, Valencia, MD).
## V7 Using a poly-dT primer incorporating a T7 promoter, double-stranded cDNA was synthesized from 10 micrograms total RNA using a cDNA synthesis kit (SuperScript double-stranded cDNA synthesis kit; Invitrogen, Carlsbad, CA).  Double-stranded cDNA was purified by phenol/chloroform extraction using a Phase-Lock Gel (PLG Heavy Brinkmann Instruments, Westbury, NY) followed by ethanol precipitation.  Biotin-labeled cRNA was generated from the double-stranded cDNA template though in-vitro transcription with T7 polymerase, and a nucleotide mix containing biotinylated CTP and UTP (Enzo RNA Transcript Labeling Kit; Enzo Diagnostics, Farmingdale, NY).  The biotinylated cRNA was purified using RNeasy affinity columns (Qiagen, Valencia, MD).
##    taxid_ch1
## V2     10090
## V3     10090
## V4     10090
## V5     10090
## V6     10090
## V7     10090
##                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 hyb_protocol
## V2 For each GeneChip, 20 micrograms of the labeled product was fragmented in 40 mM Tris-acetate, pH 8.1, 100mM KOAc, 30mM MgOAc, for 35 minutes at 94 degrees-Celsius, to an average size of 35 to 200 bases.  15 micrograms of this fragmented, biotinylated cRNA, along with hybridization controls supplied by the manufacturer (Affymetrix), were hybridized to the arrays for 16 hours at 45 degrees-Celsius and 60 rpm.  Arrays were washed and stained according to the standard Antibody Amplification for Eukaryotic Targets protocol (Affymetrix).
## V3 For each GeneChip, 20 micrograms of the labeled product was fragmented in 40 mM Tris-acetate, pH 8.1, 100mM KOAc, 30mM MgOAc, for 35 minutes at 94 degrees-Celsius, to an average size of 35 to 200 bases.  15 micrograms of this fragmented, biotinylated cRNA, along with hybridization controls supplied by the manufacturer (Affymetrix), were hybridized to the arrays for 16 hours at 45 degrees-Celsius and 60 rpm.  Arrays were washed and stained according to the standard Antibody Amplification for Eukaryotic Targets protocol (Affymetrix).
## V4 For each GeneChip, 20 micrograms of the labeled product was fragmented in 40 mM Tris-acetate, pH 8.1, 100mM KOAc, 30mM MgOAc, for 35 minutes at 94 degrees-Celsius, to an average size of 35 to 200 bases.  15 micrograms of this fragmented, biotinylated cRNA, along with hybridization controls supplied by the manufacturer (Affymetrix), were hybridized to the arrays for 16 hours at 45 degrees-Celsius and 60 rpm.  Arrays were washed and stained according to the standard Antibody Amplification for Eukaryotic Targets protocol (Affymetrix).
## V5 For each GeneChip, 20 micrograms of the labeled product was fragmented in 40 mM Tris-acetate, pH 8.1, 100mM KOAc, 30mM MgOAc, for 35 minutes at 94 degrees-Celsius, to an average size of 35 to 200 bases.  15 micrograms of this fragmented, biotinylated cRNA, along with hybridization controls supplied by the manufacturer (Affymetrix), were hybridized to the arrays for 16 hours at 45 degrees-Celsius and 60 rpm.  Arrays were washed and stained according to the standard Antibody Amplification for Eukaryotic Targets protocol (Affymetrix).
## V6 For each GeneChip, 20 micrograms of the labeled product was fragmented in 40 mM Tris-acetate, pH 8.1, 100mM KOAc, 30mM MgOAc, for 35 minutes at 94 degrees-Celsius, to an average size of 35 to 200 bases.  15 micrograms of this fragmented, biotinylated cRNA, along with hybridization controls supplied by the manufacturer (Affymetrix), were hybridized to the arrays for 16 hours at 45 degrees-Celsius and 60 rpm.  Arrays were washed and stained according to the standard Antibody Amplification for Eukaryotic Targets protocol (Affymetrix).
## V7 For each GeneChip, 20 micrograms of the labeled product was fragmented in 40 mM Tris-acetate, pH 8.1, 100mM KOAc, 30mM MgOAc, for 35 minutes at 94 degrees-Celsius, to an average size of 35 to 200 bases.  15 micrograms of this fragmented, biotinylated cRNA, along with hybridization controls supplied by the manufacturer (Affymetrix), were hybridized to the arrays for 16 hours at 45 degrees-Celsius and 60 rpm.  Arrays were washed and stained according to the standard Antibody Amplification for Eukaryotic Targets protocol (Affymetrix).
##                                                                                                                   scan_protocol
## V2 The stained GeneChip arrays were scanned at 488 nm using an Affymetrix Gene Chip Scanner 3000 (Affymetrix, Santa Clara, CA).
## V3 The stained GeneChip arrays were scanned at 488 nm using an Affymetrix Gene Chip Scanner 3000 (Affymetrix, Santa Clara, CA).
## V4 The stained GeneChip arrays were scanned at 488 nm using an Affymetrix Gene Chip Scanner 3000 (Affymetrix, Santa Clara, CA).
## V5 The stained GeneChip arrays were scanned at 488 nm using an Affymetrix Gene Chip Scanner 3000 (Affymetrix, Santa Clara, CA).
## V6 The stained GeneChip arrays were scanned at 488 nm using an Affymetrix Gene Chip Scanner 3000 (Affymetrix, Santa Clara, CA).
## V7 The stained GeneChip arrays were scanned at 488 nm using an Affymetrix Gene Chip Scanner 3000 (Affymetrix, Santa Clara, CA).
##                                                                                                    description
## V2 flow-cytometry: B220+, CD5â\200“, CD19+, CD25â\200“, CD69â\200“, B7-1â\200“, B7-2â\200“, IgD+, sIgM+; clinical-signs: none
## V3       flow-cytometry: B220+, CD5â\200“, CD19+, CD25+, CD69+, B7-1+, B7-2+, IgDlo, sIgMlo; clinical-signs: none
## V4 flow-cytometry: B220+, CD5â\200“, CD19+, CD25â\200“, CD69â\200“, B7-1â\200“, B7-2â\200“, IgD+, sIgM+; clinical-signs: none
## V5 flow-cytometry: B220+, CD5â\200“, CD19+, CD25â\200“, CD69â\200“, B7-1â\200“, B7-2â\200“, IgD+, sIgM+; clinical-signs: none
## V6                                                              flow-cytometry: ND; clinical-signs: neck tumor
## V7                                                         flow-cytometry: ND; clinical-signs: failure to nest
##                   data_processing platform_id    contact_name
## V2 MAS5.0, target intensity = 500     GPL8321 Marc,E.,Lenburg
## V3 MAS5.0, target intensity = 500     GPL8321 Marc,E.,Lenburg
## V4 MAS5.0, target intensity = 500     GPL8321 Marc,E.,Lenburg
## V5 MAS5.0, target intensity = 500     GPL8321 Marc,E.,Lenburg
## V6 MAS5.0, target intensity = 500     GPL8321 Marc,E.,Lenburg
## V7 MAS5.0, target intensity = 500     GPL8321 Marc,E.,Lenburg
##      contact_email contact_phone  contact_fax    contact_department
## V2 mlenburg@bu.edu  617-414-1375 617-414-1646 Genetics and Genomics
## V3 mlenburg@bu.edu  617-414-1375 617-414-1646 Genetics and Genomics
## V4 mlenburg@bu.edu  617-414-1375 617-414-1646 Genetics and Genomics
## V5 mlenburg@bu.edu  617-414-1375 617-414-1646 Genetics and Genomics
## V6 mlenburg@bu.edu  617-414-1375 617-414-1646 Genetics and Genomics
## V7 mlenburg@bu.edu  617-414-1375 617-414-1646 Genetics and Genomics
##                       contact_institute          contact_address
## V2 Boston University School of Medicine 715 Albany Street, E613B
## V3 Boston University School of Medicine 715 Albany Street, E613B
## V4 Boston University School of Medicine 715 Albany Street, E613B
## V5 Boston University School of Medicine 715 Albany Street, E613B
## V6 Boston University School of Medicine 715 Albany Street, E613B
## V7 Boston University School of Medicine 715 Albany Street, E613B
##    contact_city contact_state contact_zip.postal_code contact_country
## V2       Boston            MA                    2130             USA
## V3       Boston            MA                    2130             USA
## V4       Boston            MA                    2130             USA
## V5       Boston            MA                    2130             USA
## V6       Boston            MA                    2130             USA
## V7       Boston            MA                    2130             USA
##    contact_web_link
## V2 http://gg.bu.edu
## V3 http://gg.bu.edu
## V4 http://gg.bu.edu
## V5 http://gg.bu.edu
## V6 http://gg.bu.edu
## V7 http://gg.bu.edu
##                                                                   supplementary_file
## V2 ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM142nnn/GSM142398/suppl/GSM142398.CEL.gz
## V3 ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM142nnn/GSM142399/suppl/GSM142399.CEL.gz
## V4 ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM142nnn/GSM142400/suppl/GSM142400.CEL.gz
## V5 ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM142nnn/GSM142401/suppl/GSM142401.CEL.gz
## V6 ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM142nnn/GSM142402/suppl/GSM142402.CEL.gz
## V7 ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM142nnn/GSM142403/suppl/GSM142403.CEL.gz
##    data_row_count             group Ig clonality:ch1 Row.names
## V2          22690           resting       polyclonal GSM142398
## V3          22690         activated       polyclonal GSM142399
## V4          22690           resting       polyclonal GSM142400
## V5          22690           resting       polyclonal GSM142401
## V6          22690 lymphoma_marginal      oligoclonal GSM142402
## V7          22690 lymphoma_marginal      oligoclonal GSM142403
```

  We normalized the expression of genes using log2() scale. This normalization 
  will ensure a more accurate comparison of the variance between the expression
  groups (clusters).


```r
dat <- GSE6136
df <- log2(dat+1)
head(df)
```

```
##              GSM142398 GSM142399 GSM142400 GSM142401 GSM142402 GSM142403
## 1415670_at   10.219169 10.290480  9.987548 10.076816  9.827343  9.871289
## 1415671_at   10.903581 11.159934 10.776186 10.929998 11.468268 11.200408
## 1415672_at   11.115304 10.892087 11.091303 11.040290 11.109504 11.325305
## 1415673_at    8.990388 10.265615  8.742141  8.422065  8.963474  8.874674
## 1415674_a_at  9.911692 10.665869  9.942661  9.793766 10.243650 10.147078
## 1415675_at    9.524933  9.896332  9.590774  9.375474  9.392747  9.422065
##              GSM142404 GSM142405 GSM142406 GSM142407 GSM142408 GSM142409
## 1415670_at    9.675428  9.950702  9.848153 10.103419 10.040838  9.823367
## 1415671_at   10.654726 10.875135 10.854245 10.898450 10.736909 10.000000
## 1415672_at   11.639702 10.989253 11.021605 11.255088 11.278449 11.232960
## 1415673_at    8.923327  9.214562  8.965496  9.212861  9.121793  9.634993
## 1415674_a_at 10.133271 10.344962 10.392210 10.551131  9.750205 10.642864
## 1415675_at    9.396605  9.894666  9.818103  9.794741  9.857981 10.175924
##              GSM142410 GSM142411 GSM142412 GSM142413 GSM142414 GSM142415
## 1415670_at    9.974271 10.226412  9.471878  9.483413  9.743320 10.126317
## 1415671_at   10.341185 10.515897 10.486835 11.247394 10.813861 10.714589
## 1415672_at   11.227135 10.945444 10.904334 11.208173 11.316168 11.251542
## 1415673_at    9.696794  9.112440  9.050121  9.026523  8.806711  9.242221
## 1415674_a_at 10.239360 10.351381 10.469744  9.948075  9.760886 10.089980
## 1415675_at    9.603997  9.392532  9.250062  9.563387  9.597680  9.594325
##              GSM142416 GSM142417 GSM142418 GSM142419 GSM142420 GSM142421
## 1415670_at    9.654457  9.976134 10.033423  9.712699 10.093286  9.939138
## 1415671_at   10.429407 10.590681 10.663825 10.401733 10.812819 10.574026
## 1415672_at   10.418696 11.120108 11.484219 11.799808 11.337064 11.169048
## 1415673_at    8.240791 10.588715 10.380136 10.484521 10.616273  9.051209
## 1415674_a_at  9.813941 10.572984 10.664181 10.667023 10.865347  9.748696
## 1415675_at    8.961739  9.798634  9.638436  9.446670  9.690522  9.389739
##              GSM142422 GSM142423
## 1415670_at   10.152792  9.838258
## 1415671_at   10.615905 10.375908
## 1415672_at   11.469845 11.542500
## 1415673_at    9.899055 10.382732
## 1415674_a_at  9.434003  9.690696
## 1415675_at    9.111657  9.116084
```

  After normailization, we now show how to classify diffrernt stages. The
  tipping point which is within the "activated" state in this case. Here we see
  the number of samples that are clasified into states "activated",
  "lymphoma_aggressive", "lymphoma_marginal", "lymphoma_transitional" and
  "resting". For instance, states "activated" and "resting" contain three and
  four samples, respectively. All the contents of the data set "test" can be 
  viewed using `View(test)`. Each clinical state's content can be viewed using
  `head(test["stage_name"])`. For instance, head(test["activated"]) shows
  contents of the activated state.


```r
tmp <- names(table(cli$group))
samplesL <- split(cli[,1],f = cli$group)
head(samplesL)
```

```
## $activated
##        V3       V26       V27 
## GSM142399 GSM142422 GSM142423 
## 26 Levels: GSM142398 GSM142399 GSM142400 GSM142401 GSM142402 ... GSM142423
## 
## $lymphoma_aggressive
##        V9       V10       V20       V21       V22       V23       V24 
## GSM142405 GSM142406 GSM142416 GSM142417 GSM142418 GSM142419 GSM142420 
## 26 Levels: GSM142398 GSM142399 GSM142400 GSM142401 GSM142402 ... GSM142423
## 
## $lymphoma_marginal
##        V6        V7        V8       V11       V12       V13 
## GSM142402 GSM142403 GSM142404 GSM142407 GSM142408 GSM142409 
## 26 Levels: GSM142398 GSM142399 GSM142400 GSM142401 GSM142402 ... GSM142423
## 
## $lymphoma_transitional
##       V14       V15       V16       V17       V18 
## GSM142410 GSM142411 GSM142412 GSM142413 GSM142414 
## 26 Levels: GSM142398 GSM142399 GSM142400 GSM142401 GSM142402 ... GSM142423
## 
## $resting
##        V2        V4        V5       V19       V25 
## GSM142398 GSM142400 GSM142401 GSM142415 GSM142421 
## 26 Levels: GSM142398 GSM142399 GSM142400 GSM142401 GSM142402 ... GSM142423
```

```r
test <- sd_selection(df, samplesL,0.01)
head(test[["activated"]])
```

```
##              GSM142399 GSM142422 GSM142423
## 1415766_at    7.600656  8.778077  9.036723
## 1415827_a_at 11.002252 13.079218 13.205503
## 1415904_at   12.229810  5.885086  4.217231
## 1415985_at   11.786106 10.736656 10.044940
## 1416034_at   11.417114 13.474682 13.760252
## 1416071_at   11.194388 10.362273  9.993646
```

  A graphical represetation of genes of interest can be achieved using the
  functions shown below. The `getNetwork` function will obtain an igraph object
  based on a pearson correlation of `test`. This `igraphL` object is then run
  using the `getCluster_methods` function classify nodes.


```r
igraphL <- getNetwork(test, fdr = 0.05)
```

```
## [1] "activated:7 nodes"
## [1] "lymphoma_aggressive:218 nodes"
## [1] "lymphoma_marginal:0 nodes"
## [1] "lymphoma_transitional:181 nodes"
## [1] "resting:0 nodes"
```

```r
#summary(igraphL)
clstr <- getCluster_methods(igraphL)
#summary(clstr)
head(clstr)
```

```
## $activated
## IGRAPH clustering walktrap, groups: 3, mod: 0.56
## + groups:
##   $`1`
##   [1] "1418580_at" "1419357_at"
##   
##   $`2`
##   [1] "1416632_at"   "1428476_a_at" "1449553_at"  
##   
##   $`3`
##   [1] "1416930_at"   "1422032_a_at"
##   
## 
## $lymphoma_aggressive
## IGRAPH clustering walktrap, groups: 6, mod: 0.33
## + groups:
##   $`1`
##    [1] "1416001_a_at" "1416002_x_at" "1416527_at"   "1416754_at"  
##    [5] "1416771_at"   "1417133_at"   "1417237_at"   "1417435_at"  
##    [9] "1417570_at"   "1417604_at"   "1417634_at"   "1417720_at"  
##   [13] "1417730_at"   "1417826_at"   "1417839_at"   "1418563_at"  
##   [17] "1418718_at"   "1418902_at"   "1419149_at"   "1419259_at"  
##   [21] "1419469_at"   "1419657_a_at" "1419997_at"   "1420553_x_at"
##   [25] "1420897_at"   "1421066_at"   "1422181_at"   "1422466_at"  
##   [29] "1422823_at"   "1422824_s_at" "1423266_at"   "1423438_at"  
##   [33] "1423841_at"   "1423909_at"   "1424033_at"   "1424072_at"  
##   + ... omitted several groups/vertices
## 
## $lymphoma_marginal
## [1] NA
## 
## $lymphoma_transitional
## IGRAPH clustering walktrap, groups: 30, mod: 0.34
## + groups:
##   $`1`
##   [1] "1417866_at" "1424822_at" "1426957_at"
##   
##   $`2`
##    [1] "1416246_a_at" "1416950_at"   "1416992_at"   "1417111_at"  
##    [5] "1417640_at"   "1417930_at"   "1418842_at"   "1419104_at"  
##    [9] "1421922_at"   "1423226_at"   "1423526_at"   "1423829_at"  
##   [13] "1423946_at"   "1425553_s_at" "1429563_x_at" "1434883_at"  
##   [17] "1448511_at"   "1448904_at"   "1450355_a_at" "1452132_at"  
##   [21] "1452430_s_at" "1452431_s_at" "1460188_at"   "1460407_at"  
##   + ... omitted several groups/vertices
## 
## $resting
## [1] NA
```

  Finally, we can view a graph of classified clustered samples for five
  different stages. Use `"head(maxCIms)"` to view the CI scores calculated.


```r
#par(mfrow <- c(1,5))
#membersL_noweight <- getMCI(cluster, test, adjust.size = F, ylim = c(0,8))
#head(getMCI)
#maxCIms <- getMaxCImember(membersL_noweight[[1]],membersL_noweight[[2]],min =3)
#head(maxCIms)
```

<a name="">[GoTop](#)</a>
<a name="Acknowledgements"></a>\

#### __Acknowledgements__

  The development of this package would not be possible without continous help
  and feedback from individuals and institutions including: The Bioconductor
  Core Team, Dr. Xianan H Yang, Dr. Tzintzuni Garcia, and National Institutes of
  Health R21LM012619.

<a name="">[GoTop](#)</a>
<a name="SessionInfo"></a>\


```r
sessionInfo()
```

```
## R version 3.6.1 (2019-07-05)
## Platform: x86_64-w64-mingw32/x64 (64-bit)
## Running under: Windows 10 x64 (build 18362)
## 
## Matrix products: default
## 
## locale:
## [1] LC_COLLATE=C                          
## [2] LC_CTYPE=English_United States.1252   
## [3] LC_MONETARY=English_United States.1252
## [4] LC_NUMERIC=C                          
## [5] LC_TIME=English_United States.1252    
## 
## attached base packages:
## [1] parallel  stats4    stats     graphics  grDevices utils     datasets 
## [8] methods   base     
## 
## other attached packages:
##  [1] cluster_2.1.0        igraph_1.2.4.1       psych_1.8.12        
##  [4] stringr_1.4.0        NPS_0.99.0           GenomicRanges_1.37.7
##  [7] GenomeInfoDb_1.21.1  IRanges_2.19.5       S4Vectors_0.23.3    
## [10] BiocGenerics_0.31.2 
## 
## loaded via a namespace (and not attached):
##  [1] xfun_0.8               remotes_2.1.0          purrr_0.3.2           
##  [4] lattice_0.20-38        testthat_2.1.1         usethis_1.5.1         
##  [7] htmltools_0.3.6        yaml_2.2.0             rlang_0.4.0           
## [10] pkgbuild_1.0.3         foreign_0.8-71         glue_1.3.1            
## [13] withr_2.1.2            sessioninfo_1.1.1      GenomeInfoDbData_1.2.1
## [16] zlibbioc_1.31.0        commonmark_1.7         devtools_2.1.0        
## [19] memoise_1.1.0          evaluate_0.14          knitr_1.23            
## [22] callr_3.3.0            ps_1.3.0               highr_0.8             
## [25] Rcpp_1.0.1             backports_1.1.4        desc_1.2.0            
## [28] pkgload_1.0.2          XVector_0.25.0         fs_1.3.1              
## [31] mnormt_1.5-5           digest_0.6.20          stringi_1.4.3         
## [34] processx_3.4.0         rprojroot_1.3-2        grid_3.6.1            
## [37] cli_1.1.0              tools_3.6.1            bitops_1.0-6          
## [40] magrittr_1.5           RCurl_1.95-4.12        pkgconfig_2.0.2       
## [43] crayon_1.3.4           xml2_1.2.0             prettyunits_1.0.2     
## [46] assertthat_0.2.1       rmarkdown_1.13         roxygen2_6.1.1        
## [49] rstudioapi_0.10        R6_2.4.0               nlme_3.1-140          
## [52] compiler_3.6.1
```
<a name="">[GoTop](#)</a>
<a name="References"></a>\

#### __References__ 
1. Scheffer M, Carpenter SR, Lenton TM, Bascompte J, Brock W, Dakos V, et al. Anticipating critical transitions. Science. 2012;338(6105):344-8. doi: 10.1126/science.1225244. PubMed PMID: 23087241.
2. Chen L, Liu R, Liu ZP, Li M, Aihara K. Detecting early-warning signals for sudden deterioration of complex diseases by dynamical network biomarkers. Sci Rep. 2012;2:342. Epub 2012/03/31. doi: 10.1038/srep00342. PubMed PMID: 22461973; PubMed Central PMCID: PMC3314989.
3. Mojtahedi M, Skupin A, Zhou J, Castano IG, Leong-Quong RY, Chang H, et al. Cell Fate Decision as High-Dimensional Critical State Transition. PLoS Biol. 2016;14(12):e2000640. doi: 10.1371/journal.pbio.2000640. PubMed PMID: 28027308; PubMed Central PMCID: PMCPMC5189937.
4. Wang ZZ, Cunningham JM, Yang XH. CisPi: a transcriptomic score for disclosing cis-acting disease-associated lincRNAs. Bioinformatics. 2018;34(17):664-70. doi: 10.1093/bioinformatics/bty574. PubMed PMID: WOS:000444317200009.



