### Introduction
Abrupt and irreversible changes (or tipping points) are decisive in the progression of biological processes. BioTIP is an R-package for characterization of biological tipping-points. BioTIP is the first toolset that amalgamates two computational impediments: detection of tipping-points accurately, and identification of non-stochastic critical transition signals (CTSs). 

Compared to existing tipping-point methods, BioTIP addresses three analytical challenges below:  

1. Large variability often exists in the sizes of statistical ensembles (populations);  
2. Not only multiple tipping points may coexist during an observed progression but also multiple CTSs may coexist in the same critical transition state; and  
3. There are multiple trajectories when cells are exposed to stimulus.

To adopt tipping-point theory to transcriptomic analysis, there are two common premises:  

1. The system of the individual cell has a dissipative structure (e.g., having discrete states 
2. Including the one showing multi-stability);  
each state has a characteristic gene expression profile and thus presents a distinct molecular phenotype.  

BioTIP apply to data meeting these premises, including both single-cell and bulk-cell transcriptomes. The impending-transition states are called ‘critical transition states’ or ‘tipping points.’

We have applied BioTIP to identified temporal features of molecular-network dynamics from gene expression profiles successfully. Importantly, the CTS identifications have helped infer the underlying gene-regulatory network and the involving key transcription factors.


![](imgs/Fig1_scRNARNAseq_pipeline_2021_xy.jpg)

Overview of BioTIP.

![](imgs/FigS1_BioTIP_pipeline_detailed_v7.jpg)

### Installation
To use the newest BioTIP package, either clone/download this repository, or you can install BioTIP with:

```r
library("devtools")
devtools::install_github("xyang2uchicago/BioTIP")
```

You can install the released version of BioTIP from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("BioTIP")
```
or even better
``` r
source('http://bioconductor.org/biocLite.R')
biocLite("BioTIP")
```

### Example
This is a basic example which shows you how to solve a common problem:

```{r example}
library(BioTIP)
data(GSE6136_matrix)
row.names(GSE6136_matrix) = GSE6136_matrix$ID_REF
GSE6136 = GSE6136_matrix[,-1]

data(GSE6136_cli)
cli = t(GSE6136_cli)

library(stringr)
colnames(cli) = str_split_fixed(cli[1,],'_',2)[,2]
cli = cli[-1,]
cli = data.frame(cli)
cli[,"cell-type:ch1"] = str_split_fixed(cli$characteristics_ch1.1,": ",2)[,2]
cli[,"Ig clonality:ch1"] = str_split_fixed(cli$characteristics_ch1.3,": ",2)[,2]

colnames(cli)[colnames(cli) == "cell-type:ch1"] = "group"
cli$Row.names = cli[,1]

df <- log2(dat+1)

tmp <- names(table(cli$group))
samplesL <- split(cli[,1],f = cli$group)
test <- sd_selection(df, samplesL,0.01)
igraphL <- getNetwork(test, fdr = 0.05)
cluster <- getCluster_methods(igraphL)
membersL_noweight <- getMCI(cluster,test,adjust.size = F)
plotBar_MCI(membersL_noweight,ylim = c(0,6))
```

Here are example analysis with bulk and single-cell data sets:
[BioTIP Applications](https://github.com/xyang2uchicago/BioTIP_application).

## Enjoy!

