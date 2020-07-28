---
title: "BioTIP- an R package for characterization of Biological Tipping-Points"
author: "Zhezhen Wang, Biniam Feleke, Qier An, Antonio Feliciano and Xinan H Yang"
date: '07/28/2020'
output:
  html_document:
    fig_caption: yes
    keep_md: yes
    toc: yes
  pdf_document: default
  word_document: default
vignette: >
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteIndexEntry{BioTIP- an R package for characterization of Biological Tipping-Point}
  %\VignetteEncoding{UTF-8}
abstract: "Heterogeneous cells within a stable state are controlled by gene regulatory networks. Cell-state state transitions are not just gradually and linearly changing phenomena, but rather can be described as nonlinear critical transitions or tipping points (Moris 2016). This abrupt transition between two stable states has been well described using a tipping-point model (Scheffer 2012). When a complex system arrives at a tipping point, there exists fundamental instability in the system. Although this brings a potential risks of unwanted collapse, it also presents an opportunity for positive change through reversal of trajectory. The idea of tipping-point prediction has been previously pursued and implimented, however existing methods each have their limitations. The 'early-warning' package in R (http://www.early-warning-signals.org/) is designed for longitudinal data analysis and is applicable to ecosystems, ecology, space, and other fields. However, we found  it hard to apply to high-throughput 'omics' data. Another published method for biological tipping-point characterization is Dynamic Network Biomarker (DNB) (Chen 2012). This method relies on indexing gene-oscillation per state, and thus is not always feasible for cross-state comparisons and indication of tipping points. Another published method, Index of critical transition (Ic-score), predicts tipping point is limited by its reliance on pre-defined features (Mojtahedi 2016). The purpose of this R package BioTIP is to characterize Biological Tipping-Points from a high-throughput data matrix where rows represent molecular features such as transcripts and columns represent samples or cells in a data matrix. BioTIP implicates one new method and the two aforementioned published methods(DNB and Ic-score) to allow for flexible analysis of not only longitudinal, but also cross-sectional datasets. Besides overcoming the limitations of previous methods, BioTIP provides functions to optimize feature pre-selection and evaluate empirical significance. Additionally, BioTIP defines transcript biotypes according the GENCODE annotation, allowing for study of noncoding transcription (Wang 2018). These improvements allow for an application to cross-sectional datasets. The BioTIP scheme acts as a hybrid model to join the advantages of both DNB and Ic-scoring, and providing enhanced features for high-throughput 'omics' data analysis."
---

#### ####



<div class="figure" style="text-align: center">
<img src="FigS1_BioTIP_pipeline_detailed_v7.png" alt="Fig 1. BioTIP workflow with five key analytic steps. RTF: relative transcript fluctuation; PCC: Pearson correlation coefficient; MCI: Module-Criticality Index; Ic: index of critical transition." width="60%" />
<p class="caption">Fig 1. BioTIP workflow with five key analytic steps. RTF: relative transcript fluctuation; PCC: Pearson correlation coefficient; MCI: Module-Criticality Index; Ic: index of critical transition.</p>
</div>

#### [Standard workflow](#Standard workflow)
   * [An Identification of Critical Tipping Point](#An Identification of Critical Tipping Point)
     *  ##### [Data Preprocessing](#Data Preprocessing)
     *  ##### [Pre-selection Transcript](#Pre-selection Transcript)
     *  ##### [Network Partition](#Network Partition)
     *  ##### [Identifying Dynamic Network Biomodule](#Identifying Dynamic Network Biomodule)
     *  ##### [Finding Tipping Point](#Finding Tipping Point)
   * [An Application to scRNA-seq Data](#An Application to scRNA-seq Data)
     *  ##### [Coming soom ...](#Coming soon)  
   * [Transcript Annotation and Biotype](#Transcript Annotation and Biotype)
     *  ##### [Quick Start](#Quick Start)
     *  ##### [Genomic Data Source](#Genomic Data Source)
     *  ##### [Extracting Summary Data](#Extracting Summary Data)
     *  ##### [Loading Data](#Loading Data)
     *  ##### [Prepare GRanges Object](#Prepare GRanges Object)
   * [Acknowledgements](#Acknowledgements)
   * [SessionInfo](#SessionInfo)
   * [References](#References)

<a name="An Identification of Critical Tipping Point"></a>

######################################################################
### Standard workflow: An Identification of Critical Tipping Point ###
######################################################################

<a name="Data Preprocessing"></a>
 __Data Preprocessing__

An existing dataset, GSE6136, is used to demonstrate how our functions are applied. Samples were collected from transgenic mouse lymphomas and divided into five groups based on clinical presentation, pathology and flow cytometry (Lenburg 2007), thus belonging to cross-sectional profiles. Noticing these five group represent a control stage similar to non-transgenic B cells and four major periods of B-cell lymphomagenesis, Dr. Chen and coauthors used the DNB method to identify the pre-disease state exists around the normal activated period (P2), i.e., the system transitions to the disease state around the transitional lymphoma period (Figure S12 in publication (Chen 2012)). 
Start by installing the package 'BioTIP' and other dependent packages such as stringr, psych, and igraph if necessary. Below are some examples.


```r
# load 
source("/Users/jennifersun98/Desktop/BioTIP/R/BioTIP_update_3.4_07012020.R")
library(psych)
library(stringr)
library(BioTIP)
#> 
#> Attaching package: 'BioTIP'
#> The following objects are masked _by_ '.GlobalEnv':
#> 
#>     avg.cor.shrink, cor.shrink, getBiotypes, getCluster_methods,
#>     getCTS, getIc, getIc.new, getMaxMCImember, getMaxStats, getMCI,
#>     getNetwork, getReadthrough, optimize.sd_selection,
#>     plot_Ic_Simulation, plot_MCI_Simulation, plot_SS_Simulation,
#>     plotBar_MCI, plotIc, plotMaxMCI, sd_selection, simulation_Ic,
#>     simulation_Ic_sample, simulationMCI
library(igraph)
#> 
#> Attaching package: 'igraph'
#> The following objects are masked from 'package:stats':
#> 
#>     decompose, spectrum
#> The following object is masked from 'package:base':
#> 
#>     union
```

Once all the required packages are installed, load using `“read.table()”` function as follows. Note to change the read.table(“PATH”) when running your datasets. To check the dimension of your data set use “dim(dataset)” function. Here we show a use of dim function “dim(GSE6136)”. Notice that after editing the column and row, the dimension of our dataset changes from (22690,27) to (22690, 26) because we removed the downloaded first row after assigning it to be column-name of the final numeric data matrix.


```r
data(GSE6136_matrix)

dim(GSE6136_matrix)               #[1] 22690 genes and 27 samples (cells) 
#> [1] 22690    27
row.names(GSE6136_matrix) = GSE6136_matrix$ID_REF
GSE6136 = GSE6136_matrix[,-1]
dim(GSE6136)               #[1] 22690 genes and 26 samples (cells) 
#> [1] 22690    26
```

  The summary of GSE6136_matrix is GSE6136_cli shown below. These two kind of data files can be downloaded from GSE database


```r
#requires library(stringr)
library(BioTIP)
data(GSE6136_cli)

#dim(GSE6136_cli) #check dimension
cli = t(GSE6136_cli)

library(stringr)
colnames(cli) = str_split_fixed(cli[1,],'_',2)[,2]
cli = cli[-1,]
cli = data.frame(cli)
cli[,"cell-type:ch1"] = str_split_fixed(cli$characteristics_ch1.1,": ",2)[,2]
cli[,"Ig clonality:ch1"] = str_split_fixed(cli$characteristics_ch1.3,": ",2)[,2]

colnames(cli)[colnames(cli) == "cell-type:ch1"] = "group"
cli$Row.names = cli[,1]
head(cli[,1:3])
#>    geo_accession                status submission_date
#> V2     GSM142398 Public on Oct 28 2006     Oct 25 2006
#> V3     GSM142399 Public on Oct 28 2006     Oct 25 2006
#> V4     GSM142400 Public on Oct 28 2006     Oct 25 2006
#> V5     GSM142401 Public on Oct 28 2006     Oct 25 2006
#> V6     GSM142402 Public on Oct 28 2006     Oct 25 2006
#> V7     GSM142403 Public on Oct 28 2006     Oct 25 2006
```

  We normalized the expression of genes using log2() scale. This normalization 
  will ensure a more accurate comparison of the variance between the expression
  groups (clusters).


```r
dat <- GSE6136
df <- log2(dat+1)
head(df[,1:6])
#>              GSM142398 GSM142399 GSM142400 GSM142401 GSM142402 GSM142403
#> 1415670_at   10.219169 10.290480  9.987548 10.076816  9.827343  9.871289
#> 1415671_at   10.903581 11.159934 10.776186 10.929998 11.468268 11.200408
#> 1415672_at   11.115304 10.892087 11.091303 11.040290 11.109504 11.325305
#> 1415673_at    8.990388 10.265615  8.742141  8.422065  8.963474  8.874674
#> 1415674_a_at  9.911692 10.665869  9.942661  9.793766 10.243650 10.147078
#> 1415675_at    9.524933  9.896332  9.590774  9.375474  9.392747  9.422065
```

<a name=>[Go to Top](#)</a> 
<a name="Pre-selection Transcript"></a>

 __Pre-selection Transcript__

  Once normalized, we can now classify different stages. The tipping point
  is within the "activated" state in this case. Here we see the number of
  samples that are classified into states "activated", "lymphoma_aggressive",
  "lymphoma_marginal", "lymphoma_transitional" and "resting". For instance,
  states "activated" and "resting" contain three and four samples, respectively.
  All the contents of the data set "test" can be viewed using `View(test)`. Each
  clinical state's content can be viewed using `head(test["stage_name"])`. For
  instance, head(test["activated"]) shows contents of the activated state.


```r
cli$group = factor(cli$group,
                   levels = c('resting','activated','lymphoma_marginal','lymphoma_transitional','lymphoma_aggressive'))
samplesL <- split(cli[,"geo_accession"],f = cli$group)
lapply(samplesL, length)
#> $resting
#> [1] 5
#> 
#> $activated
#> [1] 3
#> 
#> $lymphoma_marginal
#> [1] 6
#> 
#> $lymphoma_transitional
#> [1] 5
#> 
#> $lymphoma_aggressive
#> [1] 7
test <- sd_selection(df, samplesL,0.01)
head(test[["activated"]])
#>              GSM142399 GSM142422 GSM142423
#> 1415766_at    7.600656  8.778077  9.036723
#> 1415827_a_at 11.002252 13.079218 13.205503
#> 1415904_at   12.229810  5.885086  4.217231
#> 1415985_at   11.786106 10.736656 10.044940
#> 1416034_at   11.417114 13.474682 13.760252
#> 1416071_at   11.194388 10.362273  9.993646
```

<a name=">[Go to Top](#)</a> 
<a name="Network Partition"></a>

 __Network Partition__

  A graphical represetation of genes of interest can be achieved using the
  functions shown below. The `getNetwork` function will obtain an igraph object
  based on a pearson correlation of `test`. This `igraphL` object is then run
  using the `getCluster_methods` function classify nodes.


```r
igraphL <- getNetwork(test, fdr = 1)
#> resting:226 nodes
#> activated:226 nodes
#> lymphoma_marginal:226 nodes
#> lymphoma_transitional:226 nodes
#> lymphoma_aggressive:226 nodes
cluster <- getCluster_methods(igraphL)
```

```r
names(cluster)
#> [1] "resting"               "activated"             "lymphoma_marginal"    
#> [4] "lymphoma_transitional" "lymphoma_aggressive"
head(cluster[[1]])
#> $`1`
#>   [1] "1416080_at"   "1416375_at"   "1416651_at"   "1416960_at"   "1417115_at"  
#>   [6] "1417310_at"   "1417415_at"   "1417518_at"   "1417799_at"   "1417988_at"  
#>  [11] "1418063_at"   "1418215_at"   "1418612_at"   "1418753_at"   "1418925_at"  
#>  [16] "1419077_at"   "1419121_at"   "1419429_at"   "1419708_at"   "1419763_at"  
#>  [21] "1419974_at"   "1419983_at"   "1420186_at"   "1420448_at"   "1420770_at"  
#>  [26] "1420940_x_at" "1421735_a_at" "1421785_at"   "1422162_at"   "1422202_at"  
#>  [31] "1422286_a_at" "1422298_at"   "1422333_at"   "1422894_at"   "1423126_at"  
#>  [36] "1423935_x_at" "1424123_at"   "1424217_at"   "1424366_at"   "1424380_at"  
#>  [41] "1424980_s_at" "1425035_s_at" "1425221_at"   "1425322_at"   "1425361_at"  
#>  [46] "1425475_at"   "1425649_at"   "1425651_at"   "1425659_at"   "1425699_a_at"
#>  [51] "1425892_a_at" "1426022_a_at" "1426075_at"   "1426816_at"   "1426909_at"  
#>  [56] "1427001_s_at" "1427360_at"   "1427623_at"   "1427653_at"   "1427735_a_at"
#>  [61] "1427843_at"   "1428766_at"   "1431706_at"   "1431821_a_at" "1433210_at"  
#>  [66] "1433909_at"   "1436639_at"   "1436810_x_at" "1437328_x_at" "1437531_at"  
#>  [71] "1438144_x_at" "1438252_at"   "1438777_a_at" "1447632_at"   "1448359_a_at"
#>  [76] "1449369_at"   "1449374_at"   "1449470_at"   "1449651_x_at" "1449701_at"  
#>  [81] "1449820_at"   "1449987_at"   "1450098_at"   "1450125_at"   "1450174_at"  
#>  [86] "1450249_s_at" "1450319_at"   "1450328_at"   "1450527_at"   "1450553_at"  
#>  [91] "1450578_at"   "1451671_at"   "1451766_at"   "1451807_at"   "1452023_at"  
#>  [96] "1452149_at"   "1452370_s_at" "1452480_at"   "1452524_a_at" "1452537_at"  
#> [101] "1452562_at"   "1452645_x_at" "1455905_at"   "1456188_at"   "1457899_at"  
#> 
#> $`2`
#>  [1] "1416745_x_at" "1417323_at"   "1418055_at"   "1418108_at"   "1418725_at"  
#>  [6] "1418735_at"   "1419146_a_at" "1419365_at"   "1419705_at"   "1419894_at"  
#> [11] "1419946_s_at" "1419981_at"   "1420096_at"   "1420169_at"   "1420200_at"  
#> [16] "1420580_at"   "1420658_at"   "1421350_a_at" "1421468_at"   "1421635_at"  
#> [21] "1421672_at"   "1421718_at"   "1421886_at"   "1422427_a_at" "1422477_at"  
#> [26] "1422667_at"   "1422986_at"   "1423844_s_at" "1424260_at"   "1424742_at"  
#> [31] "1424867_a_at" "1425039_at"   "1425357_a_at" "1425457_a_at" "1425807_at"  
#> [36] "1426157_a_at" "1426185_at"   "1427056_at"   "1427157_at"   "1427635_at"  
#> [41] "1427647_at"   "1427791_a_at" "1428016_a_at" "1429888_a_at" "1431382_a_at"
#> [46] "1434340_at"   "1435856_x_at" "1435962_at"   "1436012_s_at" "1438550_x_at"
#> [51] "1440499_at"   "1448001_x_at" "1448194_a_at" "1449688_at"   "1449902_at"  
#> [56] "1450147_at"   "1450214_at"   "1450482_a_at" "1451305_at"   "1451562_at"  
#> [61] "1451579_at"   "1451788_at"   "1452619_a_at" "1455282_x_at" "1456173_at"  
#> [66] "1459990_at"  
#> 
#> $`3`
#>  [1] "1416174_at"   "1417356_at"   "1417359_at"   "1417529_at"   "1417748_x_at"
#>  [6] "1418303_at"   "1418329_at"   "1418378_at"   "1418485_at"   "1419353_at"  
#> [11] "1419539_at"   "1419673_at"   "1419885_at"   "1420521_at"   "1420578_at"  
#> [16] "1420582_at"   "1421117_at"   "1421451_at"   "1421455_at"   "1421482_at"  
#> [21] "1421530_a_at" "1421790_a_at" "1421973_at"   "1422692_at"   "1422970_at"  
#> [26] "1423300_at"   "1424881_at"   "1425746_at"   "1425937_a_at" "1426039_a_at"
#> [31] "1426560_a_at" "1427018_at"   "1428740_a_at" "1430483_a_at" "1432835_at"  
#> [36] "1433810_x_at" "1434204_x_at" "1435989_x_at" "1437367_at"   "1437836_x_at"
#> [41] "1438374_x_at" "1438551_at"   "1449100_at"   "1450242_at"   "1450447_at"  
#> [46] "1451578_at"   "1452212_at"   "1452575_at"   "1453560_at"   "1453614_a_at"
#> [51] "1454041_at"   "1455098_a_at" "1456612_at"   "1456637_at"   "1456644_at"
```

<a name=">[Go to Top](#)</a> 
<a name="Identifying Dynamic Network Biomodule"></a>

 __Identifying Dynamic Network Biomodule__

Here, ‘module’ refers to a cluster of network nodes (e.g. transcripts) highly linked (e.g. by correlation). “Biomodule” refers to the module resenting a highest score called “Module-Criticality Index (MCI)” per state.

  The following step shows a graph of classified clustered samples for five
  different stages. MCI score is calculated for each module using the `getMCI`
  function. The `getMaxMCImember` function will obtain a list of modules with highest
  MCI at each stage. Use `"head(maxCIms)"` to view the MCI scores calculated. Use
  `plotMaxMCI` function to view the line plot of highest MCI score at each stage.


```r
membersL_noweight <- getMCI(cluster,test,adjust.size = FALSE)
plotBar_MCI(membersL_noweight,ylim = c(0,6))
```

![](BioTIP_files/figure-html/unnamed-chunk-10-1.png)<!-- -->

```r
maxMCIms <- getMaxMCImember(membersL_noweight[[1]],membersL_noweight[[2]],min =10)
names(maxMCIms)
#> [1] "idx"     "members"
names(maxMCIms[[1]])
#> [1] "resting"               "activated"             "lymphoma_marginal"    
#> [4] "lymphoma_transitional" "lymphoma_aggressive"
names(maxMCIms[[2]])
#> [1] "resting"               "activated"             "lymphoma_marginal"    
#> [4] "lymphoma_transitional" "lymphoma_aggressive"
```

```r
head(maxMCIms[['idx']])
#>               resting             activated     lymphoma_marginal 
#>                     3                     3                     1 
#> lymphoma_transitional   lymphoma_aggressive 
#>                     5                     1
head(maxMCIms[['members']][['lymphoma_aggressive']])
#> [1] "1416001_a_at" "1416002_x_at" "1416527_at"   "1416754_at"   "1416771_at"  
#> [6] "1417133_at"
```


  To get the selected statistics of biomodules (the module that has the highest MCI score) of each state, please run the following 


```r
biomodules = getMaxStats(membersL_noweight[['members']],maxMCIms[[1]])
maxMCI = getMaxStats(membersL_noweight[['MCI']],maxMCIms[[1]])
maxMCI = maxMCI[order(maxMCI,decreasing=TRUE)]
head(maxMCI)
#>             activated   lymphoma_aggressive               resting 
#>              3.522785              3.334630              2.886293 
#> lymphoma_transitional     lymphoma_marginal 
#>              2.448548              2.387766
topMCI = getTopMCI(membersL_noweight[[1]],membersL_noweight[[2]],membersL_noweight[['MCI']],min =10)
head(topMCI)
#> activated 
#>  3.522785
```

```r
maxSD = getMaxStats(membersL_noweight[['sd']],maxMCIms[[1]])
head(maxSD)
#>               resting             activated     lymphoma_marginal 
#>              1.469989              2.139577              1.295115 
#> lymphoma_transitional   lymphoma_aggressive 
#>              1.790035              1.885470
```

  To get the biomodule with the highest MCI score among all states, as we call it CTS (Critical Transition Signals), please run the following.
  

```r
CTS = getCTS(topMCI, maxMCIms[[2]])
#> Length: 35
```

 Run the following to visualize the trendence of every state represented by the cluster with the highest MCI scores.

```r
par(mar = c(10,5,0,2))
plotMaxMCI(maxMCIms,membersL_noweight[[2]],states = names(samplesL),las = 2)
```

![](BioTIP_files/figure-html/unnamed-chunk-16-1.png)<!-- -->

  We then perform simulation for MCI scores based on identified signature size
  (length(CTS) ) using the `simulationMCI` function.Use `plot_MCI_simulation`
  function to visualize the result. This step usually takes 20-30mins, so here
  to save the time, we picked a small number 3 as the length of the CTS.
  

```r
simuMCI <- simulationMCI(3,samplesL,df, B=100)
#>   |                                                                              |                                                                      |   0%  |                                                                              |=                                                                     |   1%  |                                                                              |=                                                                     |   2%  |                                                                              |==                                                                    |   3%  |                                                                              |===                                                                   |   4%  |                                                                              |====                                                                  |   5%  |                                                                              |====                                                                  |   6%  |                                                                              |=====                                                                 |   7%  |                                                                              |======                                                                |   8%  |                                                                              |======                                                                |   9%  |                                                                              |=======                                                               |  10%  |                                                                              |========                                                              |  11%  |                                                                              |========                                                              |  12%  |                                                                              |=========                                                             |  13%  |                                                                              |==========                                                            |  14%  |                                                                              |==========                                                            |  15%  |                                                                              |===========                                                           |  16%  |                                                                              |============                                                          |  17%  |                                                                              |=============                                                         |  18%  |                                                                              |=============                                                         |  19%  |                                                                              |==============                                                        |  20%  |                                                                              |===============                                                       |  21%  |                                                                              |===============                                                       |  22%  |                                                                              |================                                                      |  23%  |                                                                              |=================                                                     |  24%  |                                                                              |==================                                                    |  25%  |                                                                              |==================                                                    |  26%  |                                                                              |===================                                                   |  27%  |                                                                              |====================                                                  |  28%  |                                                                              |====================                                                  |  29%  |                                                                              |=====================                                                 |  30%  |                                                                              |======================                                                |  31%  |                                                                              |======================                                                |  32%  |                                                                              |=======================                                               |  33%  |                                                                              |========================                                              |  34%  |                                                                              |========================                                              |  35%  |                                                                              |=========================                                             |  36%  |                                                                              |==========================                                            |  37%  |                                                                              |===========================                                           |  38%  |                                                                              |===========================                                           |  39%  |                                                                              |============================                                          |  40%  |                                                                              |=============================                                         |  41%  |                                                                              |=============================                                         |  42%  |                                                                              |==============================                                        |  43%  |                                                                              |===============================                                       |  44%  |                                                                              |================================                                      |  45%  |                                                                              |================================                                      |  46%  |                                                                              |=================================                                     |  47%  |                                                                              |==================================                                    |  48%  |                                                                              |==================================                                    |  49%  |                                                                              |===================================                                   |  50%  |                                                                              |====================================                                  |  51%  |                                                                              |====================================                                  |  52%  |                                                                              |=====================================                                 |  53%  |                                                                              |======================================                                |  54%  |                                                                              |======================================                                |  55%  |                                                                              |=======================================                               |  56%  |                                                                              |========================================                              |  57%  |                                                                              |=========================================                             |  58%  |                                                                              |=========================================                             |  59%  |                                                                              |==========================================                            |  60%  |                                                                              |===========================================                           |  61%  |                                                                              |===========================================                           |  62%  |                                                                              |============================================                          |  63%  |                                                                              |=============================================                         |  64%  |                                                                              |==============================================                        |  65%  |                                                                              |==============================================                        |  66%  |                                                                              |===============================================                       |  67%  |                                                                              |================================================                      |  68%  |                                                                              |================================================                      |  69%  |                                                                              |=================================================                     |  70%  |                                                                              |==================================================                    |  71%  |                                                                              |==================================================                    |  72%  |                                                                              |===================================================                   |  73%  |                                                                              |====================================================                  |  74%  |                                                                              |====================================================                  |  75%  |                                                                              |=====================================================                 |  76%  |                                                                              |======================================================                |  77%  |                                                                              |=======================================================               |  78%  |                                                                              |=======================================================               |  79%  |                                                                              |========================================================              |  80%  |                                                                              |=========================================================             |  81%  |                                                                              |=========================================================             |  82%  |                                                                              |==========================================================            |  83%  |                                                                              |===========================================================           |  84%  |                                                                              |============================================================          |  85%  |                                                                              |============================================================          |  86%  |                                                                              |=============================================================         |  87%  |                                                                              |==============================================================        |  88%  |                                                                              |==============================================================        |  89%  |                                                                              |===============================================================       |  90%  |                                                                              |================================================================      |  91%  |                                                                              |================================================================      |  92%  |                                                                              |=================================================================     |  93%  |                                                                              |==================================================================    |  94%  |                                                                              |==================================================================    |  95%  |                                                                              |===================================================================   |  96%  |                                                                              |====================================================================  |  97%  |                                                                              |===================================================================== |  98%  |                                                                              |===================================================================== |  99%  |                                                                              |======================================================================| 100%Done!
plot_MCI_Simulation(topMCI[1],simuMCI,las=2)
```

![](BioTIP_files/figure-html/unnamed-chunk-17-1.png)<!-- -->
<a name="">[Go to Top](#)</a> 


<a name="Finding Tipping Point"></a>
 __Finding Tipping Point__

The next step is to calculate an Index of Critical transition (Ic score) of the dataset. First, use the getIc function to calculate the Ic score based on the biomodule previously identified. We use the plotIc function to draw a line plot of the Ic score. 
  

```r
IC <- getIc(df,samplesL,CTS[[1]])
par(mar = c(10,5,0,2))
plotIc(IC,las = 2)
```

![](BioTIP_files/figure-html/unnamed-chunk-18-1.png)<!-- -->

Then use the two functions to evaluate two types of empirical significance,
respectively. The first function simulation_Ic calculates random Ic-scores by
shuffling features (transcripts). Showing in the plot is Ic-score of the
identification (red) against its corresponding size-controlled random scores
(grey).


```r
simuIC <- simulation_Ic(length(CTS[[1]]),samplesL,df,B=100)
#>   |                                                                              |                                                                      |   0%  |                                                                              |=                                                                     |   1%  |                                                                              |=                                                                     |   2%  |                                                                              |==                                                                    |   3%  |                                                                              |===                                                                   |   4%  |                                                                              |====                                                                  |   5%  |                                                                              |====                                                                  |   6%  |                                                                              |=====                                                                 |   7%  |                                                                              |======                                                                |   8%  |                                                                              |======                                                                |   9%  |                                                                              |=======                                                               |  10%  |                                                                              |========                                                              |  11%  |                                                                              |========                                                              |  12%  |                                                                              |=========                                                             |  13%  |                                                                              |==========                                                            |  14%  |                                                                              |==========                                                            |  15%  |                                                                              |===========                                                           |  16%  |                                                                              |============                                                          |  17%  |                                                                              |=============                                                         |  18%  |                                                                              |=============                                                         |  19%  |                                                                              |==============                                                        |  20%  |                                                                              |===============                                                       |  21%  |                                                                              |===============                                                       |  22%  |                                                                              |================                                                      |  23%  |                                                                              |=================                                                     |  24%  |                                                                              |==================                                                    |  25%  |                                                                              |==================                                                    |  26%  |                                                                              |===================                                                   |  27%  |                                                                              |====================                                                  |  28%  |                                                                              |====================                                                  |  29%  |                                                                              |=====================                                                 |  30%  |                                                                              |======================                                                |  31%  |                                                                              |======================                                                |  32%  |                                                                              |=======================                                               |  33%  |                                                                              |========================                                              |  34%  |                                                                              |========================                                              |  35%  |                                                                              |=========================                                             |  36%  |                                                                              |==========================                                            |  37%  |                                                                              |===========================                                           |  38%  |                                                                              |===========================                                           |  39%  |                                                                              |============================                                          |  40%  |                                                                              |=============================                                         |  41%  |                                                                              |=============================                                         |  42%  |                                                                              |==============================                                        |  43%  |                                                                              |===============================                                       |  44%  |                                                                              |================================                                      |  45%  |                                                                              |================================                                      |  46%  |                                                                              |=================================                                     |  47%  |                                                                              |==================================                                    |  48%  |                                                                              |==================================                                    |  49%  |                                                                              |===================================                                   |  50%  |                                                                              |====================================                                  |  51%  |                                                                              |====================================                                  |  52%  |                                                                              |=====================================                                 |  53%  |                                                                              |======================================                                |  54%  |                                                                              |======================================                                |  55%  |                                                                              |=======================================                               |  56%  |                                                                              |========================================                              |  57%  |                                                                              |=========================================                             |  58%  |                                                                              |=========================================                             |  59%  |                                                                              |==========================================                            |  60%  |                                                                              |===========================================                           |  61%  |                                                                              |===========================================                           |  62%  |                                                                              |============================================                          |  63%  |                                                                              |=============================================                         |  64%  |                                                                              |==============================================                        |  65%  |                                                                              |==============================================                        |  66%  |                                                                              |===============================================                       |  67%  |                                                                              |================================================                      |  68%  |                                                                              |================================================                      |  69%  |                                                                              |=================================================                     |  70%  |                                                                              |==================================================                    |  71%  |                                                                              |==================================================                    |  72%  |                                                                              |===================================================                   |  73%  |                                                                              |====================================================                  |  74%  |                                                                              |====================================================                  |  75%  |                                                                              |=====================================================                 |  76%  |                                                                              |======================================================                |  77%  |                                                                              |=======================================================               |  78%  |                                                                              |=======================================================               |  79%  |                                                                              |========================================================              |  80%  |                                                                              |=========================================================             |  81%  |                                                                              |=========================================================             |  82%  |                                                                              |==========================================================            |  83%  |                                                                              |===========================================================           |  84%  |                                                                              |============================================================          |  85%  |                                                                              |============================================================          |  86%  |                                                                              |=============================================================         |  87%  |                                                                              |==============================================================        |  88%  |                                                                              |==============================================================        |  89%  |                                                                              |===============================================================       |  90%  |                                                                              |================================================================      |  91%  |                                                                              |================================================================      |  92%  |                                                                              |=================================================================     |  93%  |                                                                              |==================================================================    |  94%  |                                                                              |==================================================================    |  95%  |                                                                              |===================================================================   |  96%  |                                                                              |====================================================================  |  97%  |                                                                              |===================================================================== |  98%  |                                                                              |===================================================================== |  99%  |                                                                              |======================================================================| 100%Done!
par(mar = c(10,5,0,2))
plot_Ic_Simulation(IC,simuIC,las = 2)
```

![](BioTIP_files/figure-html/unnamed-chunk-19-1.png)<!-- -->
  
Another function plot_simulation_sample calculates random Ic-scores by shuffling
samples and visualizes the results. Showing in the plot is observed Ic-score
(red vertical line) comparing to the density of random scores (grey), at the
tipping point identified.


```r
sample_Ic = simulation_Ic_sample(df, sampleNo=3, genes=CTS[[1]], plot=TRUE)
```

![](BioTIP_files/figure-html/unnamed-chunk-20-1.png)<!-- -->

```r
#simulated_Ic = plot_simulation_sample(df,length(samplesL[['lymphoma_aggressive']]),IC[['lym#phoma_aggressive']],CTS) 
```

<a name=">[Go to top](#)</a>
<a name="An Application to scRNA-seq Data"></a>


############################################
###   An Application to scRNA-seq Data   ###
############################################

<a name=">[Go to top](#)</a>
<a name="Coming soon"></a>
Coming soon ...
                     

#############################################
###   Transcript Annotation and Biotype   ###
#############################################

<a name="Quick Start"></a>

 __Quick Start__

  The R function getBiotype is used to group transcripts of interest into 11
  biotypes based on GENCODE annotation (Fig 2a). When a query transcript
  overlaps at least half with a GENCODE transcript on the same strand, this
  query transcript will inherit the biotype of the GENCODE transcript.

  In the previous study conducted, five out of the 11 biotypes showed high
  protein-coding potential while the others did not (Fig 2b) [4]. We thus
  concluded these seven other biotypes, including protein-coding antisense RNAs,
  to be lncRNAs. The remaining coding biotypes in this study included canonic
  protein coding (CPC), ‘PC_mixed’, and ‘PC_intron’.

  First start by loading the required libraries: “GenomeInfoDb,” “BioTIP,”
  “GenomicRanges,” “IRanges” and “BioTIP”. Next load the datasets: “gencode”,
  “ILEF”, “intron” and “cod”. Using these datasets, excute BioTIP functions
  getBiotypes and getReadthrough as follows. These steps assume you installed
  the “BioTIP” package. If you did not install the package, use the
  `install.packages("BioTIP")` to install in R.

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
library(BioTIP)

data(gencode)
head(gencode)
#> GRanges object with 6 ranges and 1 metadata column:
#>       seqnames          ranges strand |  biotype
#>          <Rle>       <IRanges>  <Rle> | <factor>
#>   [1]    chr21 9683191-9683272      + |    miRNA
#>   [2]    chr21 9683191-9683272      + |    miRNA
#>   [3]    chr21 9683191-9683272      + |    miRNA
#>   [4]    chr21 9825832-9826011      + |    miRNA
#>   [5]    chr21 9825832-9826011      + |    miRNA
#>   [6]    chr21 9825832-9826011      + |    miRNA
#>   -------
#>   seqinfo: 25 sequences from an unspecified genome; no seqlengths
```

  These illustrations above assumes you have installed "BioTIP" package. If you did
  not install the package already, use the `install.packages("BioTIP")` to install in R.

<a name="Genomic Data Source"></a>\

 __Genomic Data Source__ 

High quality human genome sequence data can be obtained from various sources. To
demonstrate this package, we obtained a comprehensive gene annotation of human
GRCh37 from [GENCODE](https://www.gencodegenes.org/human). For our
illustrations, human GRCh37 data will be used. A standard file structure,
similar to general transfer format (gtf) format, is required for this package.
This gtf file organizes genomic data in rows and columns (fields). Each row
contains information about specific samples. The columns are tab separated
headers of the data frame.  There are eight fixed columns with specific headers.
An example of gtf format is shown below. For details of the gtf file format
visit this
[link](https://useast.ensembl.org/info/website/upload/gff.html#tracklines
target="_blank").


   ![Chromosome 21 of human GRCh37 gtf](chr21.JPG)

The table above contains a chr21 dataset which was extracted from a full genome
dataset. An extraction method for filtering chr21 from `gencode` file is
described below.


<a name=>[Go to Top](#)</a> 
<a name="Extracting Summary Data"></a>\

 __Extracting Summary Data__

  Before any further analysis, we need to summarize the content of the raw gtf
  data. There are two ways to get genome biotypes: a) "transcript_type" b)
  "gene_type". Due to our interst in coding and noncoding regions, the
  `transcript_type` method was used to extract the regions of interest using
  python script shown below. __Note__ that the `"PATH_FILE"` refers to the path
  where the downloded gtf file is located. For instance, if the gtf file is
  located on your `desktop`, replace the `"PATH_FILE"` Cc
  `"Users/user/Desktop/gtf"`.

**Python codes:**


```r
gtf = ("Your/PATH/TO/THE/FILE")
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

<a name="Home">[Go to Top](#)</a>
<a name="Loading Data"></a>\ 

 __Loading Data__

In order to load your data from a local drive, use the following format.
  __Note__ that the `"PATH_FILE"` refers to the location of the summary data
  from the above section. For more details on how to load datasets click
  [here](https://support.rstudio.com/hc/en-us/articles/218611977-Importing-Data-with-RStudio).


##### loading data from local drive
 > data <- read.delim("PATH_FILE", comment.char = "#")

  Internal BioTIP package data is included in the data folder. The data can be
  loaded into R working console using `data()`function. Here we show an example
  of how to load a dataset `gencode` from the data directory. A quick view of
  the data can be achieved using `head(gencode)`.


```r
library(BioTIP)
library(GenomicRanges)
#> Loading required package: stats4
#> Loading required package: BiocGenerics
#> Loading required package: parallel
#> 
#> Attaching package: 'BiocGenerics'
#> The following objects are masked from 'package:parallel':
#> 
#>     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
#>     clusterExport, clusterMap, parApply, parCapply, parLapply,
#>     parLapplyLB, parRapply, parSapply, parSapplyLB
#> The following objects are masked from 'package:igraph':
#> 
#>     normalize, path, union
#> The following objects are masked from 'package:stats':
#> 
#>     IQR, mad, sd, var, xtabs
#> The following objects are masked from 'package:base':
#> 
#>     anyDuplicated, append, as.data.frame, basename, cbind, colnames,
#>     dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
#>     grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
#>     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
#>     rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
#>     union, unique, unsplit, which, which.max, which.min
#> Loading required package: S4Vectors
#> 
#> Attaching package: 'S4Vectors'
#> The following object is masked from 'package:base':
#> 
#>     expand.grid
#> Loading required package: IRanges
#> 
#> Attaching package: 'IRanges'
#> The following object is masked from 'package:psych':
#> 
#>     reflect
#> Loading required package: GenomeInfoDb
data(gencode)
head(gencode)
#> GRanges object with 6 ranges and 1 metadata column:
#>       seqnames          ranges strand |  biotype
#>          <Rle>       <IRanges>  <Rle> | <factor>
#>   [1]    chr21 9683191-9683272      + |    miRNA
#>   [2]    chr21 9683191-9683272      + |    miRNA
#>   [3]    chr21 9683191-9683272      + |    miRNA
#>   [4]    chr21 9825832-9826011      + |    miRNA
#>   [5]    chr21 9825832-9826011      + |    miRNA
#>   [6]    chr21 9825832-9826011      + |    miRNA
#>   -------
#>   seqinfo: 25 sequences from an unspecified genome; no seqlengths
```
<a name="Home">[Go to Top](#)</a>

<a name="Prepare GRanges Object"></a>\ 

 __Prepare GRanges Object__

  Here we show an extraction of "gencode" dataset using R commands. Note
  to replace `PATH_FILE` with file direcotry path. `gtf` refers to the full 
  genome file. A subset function was used to filter chr21 dataset as follows.

`chr21 <- subset(gencode, seqnames == "chr21")` #"genecode" = whole genome gtf

    > gtf = read.table("PATH_FILE")
    > gtf = subset(gtf, biotype == "transcript")
    > colnames(gtf) = c("chr","start","end","strand","biotype")
    > gr = GRanges(gtf)

<a name="Home">[Go to Top](#)</a> 


##### Processing Query
***


```r
query <- GRanges(c("chr1:2-10:+","chr1:6-10:-"), Row.names = c("trans1","trans2"), score = c(1,2))
head(query)
#> GRanges object with 2 ranges and 2 metadata columns:
#>       seqnames    ranges strand |   Row.names     score
#>          <Rle> <IRanges>  <Rle> | <character> <numeric>
#>   [1]     chr1      2-10      + |      trans1         1
#>   [2]     chr1      6-10      - |      trans2         2
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```

##### Classifying Biotypes
***


```r
library(BioTIP)
gr <- GRanges(c("chr1:1-5:+","chr1:2-3:+"),biotype = c("lincRNA","CPC"))
head(gr)
#> GRanges object with 2 ranges and 1 metadata column:
#>       seqnames    ranges strand |     biotype
#>          <Rle> <IRanges>  <Rle> | <character>
#>   [1]     chr1       1-5      + |     lincRNA
#>   [2]     chr1       2-3      + |         CPC
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```

##### Extracting intron coordinates
*** 

      # Intron coordinates
      
       intron <- GRanges("chr1:6-8:+")
  

```r

intron <- GRanges("chr1:6-8:+")
head(intron)
#> GRanges object with 1 range and 0 metadata columns:
#>       seqnames    ranges strand
#>          <Rle> <IRanges>  <Rle>
#>   [1]     chr1       6-8      +
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```

##### Filtering coding transcripts
***

    # Filtering non-coding regions using products from example 1, 2 and 3


```r

intron_trncp <- getBiotypes(query, gr, intron)
intron_trncp
#>   seqnames start end width strand Row.names score type.fullOverlap
#> 1     chr1     2  10     9      +    trans1     1          de novo
#> 2     chr1     6  10     5      -    trans2     2          de novo
#>   type.partialOverlap type.50Overlap hasIntron type.toPlot
#> 1       lincRNA,  CPC  lincRNA,  CPC       yes     lincRNA
#> 2             de novo        de novo        no     de novo
```

    # Filtering Intron and Exons 

Here we show how to obtain protein coding and non-coding from our datasets. The
coding transcripts are an expressed section of the genome that is responsible
for protein formation. Meanwhile the non-coding transcripts are vital in the
formation regulatory elements such promoters, enhancers and silencers.


```r
library(BioTIP)
data("intron")
data("ILEF")
data("gencode")

gencode_gr = GRanges(gencode)
ILEF_gr = GRanges(ILEF)
cod_gr = GRanges(cod)
intron_gr = GRanges(intron)

non_coding <- getBiotypes(ILEF_gr, gencode_gr, intron_gr)
#> Warning in .Seqinfo.mergexy(x, y): The 2 combined objects have no sequence levels in common. (Use
#>   suppressWarnings() to suppress this warning.)
dim(non_coding)
#> [1] 300  11
head(non_coding[,1:3])
#>   seqnames    start      end
#> 1    chr21 15608524 15710335
#> 2    chr21 43619799 43717938
#> 3    chr21 28208595 28217692
#> 4    chr21 28279059 28339668
#> 5    chr21 46493768 46646483
#> 6    chr21 45285030 45407475
```

```r
coding <- getBiotypes(ILEF_gr, gencode_gr)
dim(coding)
#> [1] 300  11
head(coding[,1:3])
#>   seqnames    start      end
#> 1    chr21 15608524 15710335
#> 2    chr21 43619799 43717938
#> 3    chr21 28208595 28217692
#> 4    chr21 28279059 28339668
#> 5    chr21 46493768 46646483
#> 6    chr21 45285030 45407475
```


##### Finding overlapping transcripts
***
    # Samples with overlapping coding regions.

```r
library(BioTIP)

data(ILEF)
data(cod)
ILEF_gr = GRanges(ILEF)
cod_gr = GRanges(cod)

rdthrough <- getReadthrough(ILEF_gr, cod_gr)
head(rdthrough)
#>   seqnames    start      end  width strand Row.names readthrough
#> 1    chr21 15608524 15710335 101812      +    ABCC13           0
#> 2    chr21 43619799 43717938  98140      +     ABCG1           1
#> 3    chr21 28208595 28217692   9098      -   ADAMTS1           1
#> 4    chr21 28279059 28339668  60610      -   ADAMTS5           1
#> 5    chr21 46493768 46646483 152716      +    ADARB1           1
#> 6    chr21 45285030 45407475 122446      +    AGPAT3           1
```

<a name="Acknowledgements"></a>\

 __Acknowledgements__

  The development of this package would not be possible without continous help
  and feedback from individuals and institutions including: The Bioconductor
  Core Team, Dr. Xianan H Yang, Dr. Tzintzuni Garcia, and National Institutes of
  Health R21LM012619.

<a name="SessionInfo"></a>\


```r
sessionInfo()
#> R version 4.0.0 (2020-04-24)
#> Platform: x86_64-apple-darwin17.0 (64-bit)
#> Running under: macOS Mojave 10.14.6
#> 
#> Matrix products: default
#> BLAS:   /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRblas.dylib
#> LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib
#> 
#> locale:
#> [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
#> 
#> attached base packages:
#> [1] parallel  stats4    stats     graphics  grDevices utils     datasets 
#> [8] methods   base     
#> 
#> other attached packages:
#> [1] GenomicRanges_1.40.0 GenomeInfoDb_1.24.0  IRanges_2.22.1      
#> [4] S4Vectors_0.26.1     BiocGenerics_0.34.0  igraph_1.2.5        
#> [7] BioTIP_1.2.0         stringr_1.4.0        psych_1.9.12.31     
#> 
#> loaded via a namespace (and not attached):
#>  [1] Rcpp_1.0.4.6           XVector_0.28.0         knitr_1.28            
#>  [4] cluster_2.1.0          magrittr_1.5           zlibbioc_1.34.0       
#>  [7] MASS_7.3-51.6          mnormt_1.5-7           lattice_0.20-41       
#> [10] rlang_0.4.6            highr_0.8              tools_4.0.0           
#> [13] grid_4.0.0             nlme_3.1-147           xfun_0.14             
#> [16] htmltools_0.4.0        yaml_2.2.1             digest_0.6.25         
#> [19] GenomeInfoDbData_1.2.3 bitops_1.0-6           RCurl_1.98-1.2        
#> [22] evaluate_0.14          rmarkdown_2.1          stringi_1.4.6         
#> [25] compiler_4.0.0         pkgconfig_2.0.3
```
<a name=">[Go to Top](#)</a>
<a name="References"></a>\

 __References__ 

* Scheffer M, Carpenter SR, Lenton TM, Bascompte J, Brock W, Dakos V, et al. Anticipating critical transitions. Science. 2012;338(6105):344-8. doi: 10.1126/science.1225244. PubMed PMID: 23087241.
* Chen L, Liu R, Liu ZP, Li M, Aihara K. Detecting early-warning signals for sudden deterioration of complex diseases by dynamical network biomarkers. Sci Rep. 2012;2:342. Epub 2012/03/31. doi: 10.1038/srep00342. PubMed PMID: 22461973; PubMed Central PMCID: PMC3314989.
* Lenburg, M. E., A. Sinha, D. V. Faller and G. V. Denis (2007). "Tumor-specific and proliferation-specific gene expression typifies murine transgenic B cell lymphomagenesis." J Biol Chem 282(7): 4803-4811.PMC2819333
* Moris, N., C. Pina and A. M. Arias (2016). “Transition states and cell fate decisions in epigenetic landscapes.” Nat Rev Genet 17(11): 693-703. PMID: 27616569.
* Mojtahedi M, Skupin A, Zhou J, Castano IG, Leong-Quong RY, Chang H, et al. Cell Fate Decision as High-Dimensional Critical State Transition. PLoS Biol. 2016;14(12):e2000640. doi: 10.1371/journal.pbio.2000640. PubMed PMID: 28027308; PubMed Central PMCID: PMCPMC5189937.
* Wang, Z. Z., J. M. Cunningham and X. H. Yang (2018). “CisPi: a transcriptomic score for disclosing cis-acting disease-associated lincRNAs.” Bioinformatics34(17): 664-670"



