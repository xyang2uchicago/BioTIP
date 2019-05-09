#'
#' getReadthrough function
#'
#' @description
#' A getReadthrough functions that will The getReadthrough function is used to
#' find long transcripts that covers more than two coding regions of a genome.
#' The variables are as follows.
#'
#' @param gr A GRanges object that shows the start and end loci on genome.
#' @param cod_gr A GRanges object that contains coding regions. For details
#'   please visit /R/data.R.
#'
#' @details For details of findOverlaps, type.partialOverlap, type.50Overlap
#'   type.toPlot, queryhits, and subjecthits see,
#'   [GenomicRanges](https://www.bioconductor.org/packages/release/bioc/html/GenomicRanges.html),
#'   [IRanges](https://www.bioconductor.org/packages/release/bioc/html/IRanges.html),
#'    and [BiocManager](http://bioconductor.org/install/index.html).
#'
#' @return Returns the classified transcriptome biotypes.
#'
#' @source [Refrence GRCh37 genome](https://www.gencodegenes.org/human/release_25lift37.html)
#' for details on gtf format visit [ensemble](https://useast.ensembl.org/info/website/upload/gff.html)
#'
#' @import GenomicRanges
#' @importFrom stats aggregate
#'
#' @references
#'
#' Wang, Z. Z., J. M. Cunningham and X. H. Yang (2018).'CisPi: a transcriptomic score for disclosing cis-acting disease-associated lincRNAs.'
#' Bioinformatics34(17): 664-670'
#'
#' @examples
#' data('gencode_gr.v19_chr21.rda')
#' data('ILEF_gr.chr21.rda')
#' data('cod_gr.chr21.rda')
#' getReadthrough(ILEF_gr,cod_gr)
#'
#' \dontrun{getReadthrough(cod_gr)}
#'
#' @note
#' Replace the path_file when loading data locally to the data directory.
#'
#' @import GenomicRanges IRanges
#' @importFrom stats aggregate
#' @export

#gencode_gr <- data("gencode_gr.v19_chr21")
#cod_gr <- subset(gencode_gr, biotype = "protein_coding")
#cod_gr <- data("cod_gr.chr21")
#
getReadthrough <- function(gr, cod_gr) {
    full_table = data.frame(gr)
    overlapcount = countOverlaps(gr, cod_gr)
    completeoverlap = unique(subjectHits(findOverlaps(cod_gr, GRanges(full_table$ID), type = "within")))
    if (length(completeoverlap) == 0) {
        full_table$readthrough = ifelse(overlapcount > 2, 1, 0)
    } else {
        full_table$readthrough = ifelse(overlapcount > 2 & row.names(completeoverlap) %in% completeoverlap,
            1, 0)
    }
    gr = GRanges(subset(full_table, readthrough == 1))
    idx = subset(full_table, readthrough == 1)$ID
    overlaps = as.data.frame(findOverlaps(gr, cod_gr))
    splitoverlaps = split(overlaps, f = overlaps$queryHits)
    table(sapply(splitoverlaps, nrow) > 1)
    cod_grL = sapply(splitoverlaps, function(x) cod_gr[x$subjectHits])
    overlapL = sapply(cod_grL, function(x) findOverlaps(x))
    notoverlap = sapply(overlapL, function(x) identical(queryHits(x), subjectHits(x)))
    full_table$readthrough = ifelse(full_table$readthrough == 1 & !notoverlap, 1, 0)
    return(full_table)
}
