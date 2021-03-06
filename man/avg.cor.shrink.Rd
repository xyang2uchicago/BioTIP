% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BioTIP_update_3.3_02282020.R
\name{avg.cor.shrink}
\alias{avg.cor.shrink}
\title{Estimation of average values of correlation}
\usage{
avg.cor.shrink(
  X,
  Y = NULL,
  MARGIN = c(1, 2),
  shrink = TRUE,
  abs = FALSE,
  target = c("zero", "average", "half")
)
}
\arguments{
\item{X}{A G1 x S matrix of counts. Rows correspond to genes,
columns correspond to samples.}

\item{Y}{A G2 x S matrix of counts. Rows correspond to genes,
columns correspond to samples. By default is NULL.}

\item{MARGIN}{An integer indicate whether the rows (1,  genes)
or the columns (2,  samples) to be calculated for pairwise correlation.}

\item{shrink}{A flag specifying whether to shrink the correlation or not.}

\item{abs}{A flag specifying whether to take the absolute value
before taking the average (used for gene-gene correlations,
not sample-sample correlations)}

\item{target}{A character choose among ('zero', 'average',  'half'),
indicating whether to shrink towards zero (for gene-gene correlations),
shrink towards their empirical common average,
one or 0.5 (for sample-sample correlations).}
}
\value{
The average pairwise correlation between genes or samples.
}
\description{
This function takes in ine (or two) matrix X
(rows are genes, columns are samples) (or Y).
It then calculates the average pairwise correlation between genes or samples.
This method uses the method outlined by Schafer and Strimmer in
"A Shrinkage Approach to Large-Scale Covariance Matrix Estimation
and Implications for Functional Genomics" (2005)
}
\examples{
## Generating a data X as coming from a multivariate normal distribution 
## with 10 highly correlated variables, roughly simulating correlated genes.
M = matrix(.9, nrow = 10, ncol = 10)
diag(M) = 1
mu = rnorm(10)
X = MASS::mvrnorm(1000, mu, M)
dim(X)  #1000 10  

## Mean value of standard pairwise correlation between 1000 genes
# cor_tX = cor(t(X))
# mean(abs(cor_tX[upper.tri(cor_tX, diag = FALSE)])) # 0.9150228

## Calculating estimated pairwise correlation between 1000 genes
avg.cor.shrink(X, MARGIN=1,shrink = TRUE, targe='zero') # 0.8287838

M = matrix(.9, nrow = 20, ncol = 10)
diag(M) = 1
Y = rbind(M,X)
dim(Y)  #1020 10 
avg.cor.shrink(X, Y, MARGIN=1,shrink = TRUE, targe='zero') #0.8197959


}
\references{
Schafer and Strimmer (2005) "A Shrinkage Approach to Large-Scale
Covariance Matrix Estimation and Implications for Functional Genomics"
}
\author{
Andrew Goldstein \email{andrewgoldstein@uchicago.edu}; Xinan H Yang \email{xyang2@uchicago.edu}
}
