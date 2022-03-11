# BioTIP: an R-package for Characterization of Biological Tipping Points
### What is BioTIP?
Abrupt and often irreversible changes (or tipping points) are decisive in the progression of biological processes. We, therefore, developed this R-package for the characterization of biological tipping points using gene expression profiles. BioTIP is the first toolset that amalgamates two computational impediments in multivariate expression-data analysis: (1) detection of tipping points accurately, and (2) identification of significant critical transition signals (CTSs). 

### Why BioTIP?
BioTIP addresses and shows robust performances by addressing three analytical challenges of the existing tipping-point methods:

1. The sizes of the statistical ensembles (e.g., cellular populations) vary considerably.
2. Multiple tipping points may coexist during an observed progression and multiple CTSs may coexist in the same critical transition state. 
3. Under exposure to a stimulus, the same population of cells faces multiple trajectories.

We applied BioTIP to six datasets and compared BioTIP's performance with other existing tools.

<img src="https://github.com/xyang2uchicago/BioTIP/blob/master/6db_for_git.jpg"> 

BioTIP also demonstrated robustness with respect to different clustering methods, as shown in the figure below. 

<img src="https://github.com/xyang2uchicago/BioTIP/blob/master/FigS1_robustness_xy_v3.jpg"> 


### Where to apply BioTIP?
To adopt tipping-point theory to transcriptomic analysis, there are two commonly accepted premises:  

1. The system of the individual cell has a dissipative structure (e.g., having discrete states, including the one showing semi-stability).
2. Each state has a characteristic gene expression profile and thus presents a distinct molecular phenotype.  

BioTIP applies to data meeting these premises, including both single-cell and bulk-cell transcriptomes. The impending transition states are called ‘critical transition states’ or ‘tipping points.’

We have successfully applied BioTIP to identify temporal features of molecular-network dynamics from gene expression profiles. Importantly, the CTS identifications helped infer the underlying gene-regulatory network and the involving key transcription factors.

### How does BioTIP work? 

1. [BioTIP tutorial](https://htmlpreview.github.io/?https://github.com/xyang2uchicago/BioTIP/blob/master/Gastrulation.html): This is a detailed walkthrough of BioTIP on one of our key results (Mouse Gastrulation, GSE87038, [E8.25 2019](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE87038)). 

2. [Vignette](https://bioconductor.org/packages/release/bioc/vignettes/BioTIP/inst/doc/BioTIP.html): This documented exampled case studies on bulk (GSE6136) and single-cell ([Nestorowa 2016](https://pubmed.ncbi.nlm.nih.gov/27365425/)) datasets.

### How to install?
To use the newest BioTIP package, either clone/download this repository, or you can install BioTIP with:

```r
library("devtools")
devtools::install_github("xyang2uchicago/BioTIP")
```

You can install the released version of BioTIP from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("BioTIP")
```
or even better:
``` r
source('http://bioconductor.org/biocLite.R')
biocLite("BioTIP")
```

### Acknowledgements
BioTIP is made possible by contributions from the following authors: Xinan H Yang, Zhezhen Wang, Andrew Goldstein, Yuxi Sun, Dannie Griggs, Antonio Feliciano, Yanqiu Wang, Biniam Feleke, Qier An, Ieva Tolkaciovaite, and John M Cunningham. 
