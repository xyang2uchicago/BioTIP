#'
#' Classifying Transcripts Using getminipackage
#'
#' @description
#' The getminipackage contains two functions called getBiotypes and
#' getReadthrough. The purpose of the getBiotypes function is to class both
#' coding and noncoding transcripts into biotypes using the most recent GENCODE
#' annotations. This tool can also be used to define potential lncRNAs, given an
#' available genome transcriptome assembly (a gtf file) or any genomic loci of
#' interest.The getReadthrough function is used to find long transcripts that
#' covers more than two coding regions of a genome.
#'
#' @references
#' Wang, Z. Z., J. M. Cunningham and X. H. Yang (2018).'CisPi: a transcriptomic
#' score for disclosing cis-acting disease-associated lincRNAs.'
#' Bioinformatics34(17): 664-670'
#' @docType package
#' @aliases getminipackage package-getminipackage
#' @import GenomicRanges
#' @importFrom stats aggregate
#
data("gencode_gr.v19_chr21")
data("ILEF_gr.chr21")
data("cod_gr.chr21")
data("intron_gr.chr21")

a <- getBiotypes(ILEF_gr, gencode_gr, intron_gr)
a

b <- getReadthrough(ILEF_gr, cod_gr)
b