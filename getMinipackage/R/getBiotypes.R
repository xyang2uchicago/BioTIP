#'
#' getBiotypes function
#'
#' @aliases getBiotypes
#'
#' @description
#' The purpose of the getBiotypes function is to class both coding and noncoding
#' transcripts into biotypes using the most recent GENCODE annotations. This
#' tool can also be used to define potential lncRNAs, given an available genome
#' transcriptome assembly (a gtf file) or any genomic loci of interest.
#'
#' @param full_gr A GRanges object of coding and noncoding transctipts. Unique
#'   identifications for each column must be assigned.More details can be found
#'   in the GRanges package.
#' @param gencode_gr This GRanges object contains a GENCODE reference
#'   annotation.It must have a column of biotypes.
#' @param intron_gr A GRanges object containing the coordinates of introns.For
#'   details see GRanges package.
#' @param minoverlap Detects minimum overlap between two IRanges objects.
#'   Details Overlap arguments are included in the IRanges package.
#'
#' @details For details of findOverlaps, type.partialOverlap, type.50Overlap
#' type.toPlot, queryhits, and subjecthits see
#' [GenomicRanges](https://www.bioconductor.org/packages/release/bioc/html/GenomicRanges.html)
#' [IRanges](https://www.bioconductor.org/packages/release/bioc/html/IRanges.html),
#' and [BiocManager](http://bioconductor.org/install/index.html).
#'
#' @return Returns the classified transcriptome biotypes.
#'
#' @source [Refrence GRCh37 genome](https://www.gencodegenes.org/human/release_25lift37.html)
#' for details on gtf format visit [ensemble](https://useast.ensembl.org/info/website/upload/gff.html)
#' @import GenomicRanges
#'
#' @references
#'
#' Wang, Z. Z., J. M. Cunningham and X. H. Yang (2018).'CisPi: a transcriptomic score for disclosing cis-acting disease-associated lincRNAs.'
#' Bioinformatics34(17): 664-670'
#'
#' @examples
#' #Input datasets locally
#' data("gencode_gr.v19_chr21.rda")
#' data("intron_gr.chr21.rda")
#' data("data/ILEF_gr.chr21.rda")
#'
#' getBiotypes(ILEF_gr,gencode_gr)
#'
#' \dontrun{getBiotypes('intron_gr')}
#'
#' @note
#' Replace the path_file when loading data locally to the data directory.
#'
#' @import GenomicRanges IRanges
#' @importFrom stats aggregate
#' @export

getBiotypes <- function(full_gr, gencode_gr, intron_gr = NULL, minoverlap = 1L) {
    ## check input format ##########
    if (class(full_gr) != "GRanges")
        stop("please give full_gr as a \"GRanges\" object")
    if (class(gencode_gr) != "GRanges")
        stop("pealse give gencode_gr as a \"GRanges\" object")
    if (class(intron_gr) != "GRanges" & !is.null(intron_gr))
        stop("please give intron_gr as a \"GRanges\" object")
    ## find transcripts overlap with any GENCODe annotated transcripts ##########
    hits = findOverlaps(full_gr, gencode_gr, type = "within", minoverlap = minoverlap)
    ## derive the unique index from full_gr
    full = as.data.frame(full_gr)
    full$type.fullOverlap = "de novo"
    idx = as.data.frame(mcols(full_gr[queryHits(hits)]))
    if (nrow(idx) != 0) {
        idx$biotype = as.data.frame(mcols(gencode_gr[subjectHits(hits)]))[, 1]
        idx_collapse = aggregate(as.list(idx["biotype"]), idx["Row.names"], FUN = function(X) paste(unique(X),
            collapse = ", "))
        idx_full = match(idx_collapse$Row.names, full$Row.names)
        full[idx_full, ]$type.fullOverlap = idx_collapse$biotype
    }
    hits = findOverlaps(full_gr, gencode_gr, minoverlap = minoverlap)
    overlaps <- pintersect(full_gr[queryHits(hits)], gencode_gr[subjectHits(hits)])
    percentOverlap <- width(overlaps)/width(gencode_gr[subjectHits(hits)])

    idx = as.data.frame(mcols(full_gr[queryHits(hits)]))
    idx$biotype = as.data.frame(mcols(gencode_gr[subjectHits(hits)]))
    idx_collapse = aggregate(as.list(idx["biotype"]), idx["Row.names"], FUN = function(X) paste(unique(X),
        collapse = ", "))

    full$type.partialOverlap = "de novo"
    idx_partial = match(idx_collapse$Row.names, full$Row.names)
    full[idx_partial, ]$type.partialOverlap = idx_collapse$biotype

    idx$percentOverlap = percentOverlap
    idx_50 = subset(idx, percentOverlap >= 0.5)
    idx_50collapse = aggregate(as.list(idx_50["biotype"]), idx_50["Row.names"], FUN = function(X) paste(unique(X),
        collapse = ", "))

    full$type.50Overlap = "de novo"
    idx_50 = match(idx_50collapse$Row.names, full$Row.names)
    full[idx_50, ]$type.50Overlap = idx_50collapse$biotype

    if (!is.null(intron_gr)) {
        ## check introns ##########
        hits = findOverlaps(full_gr, intron_gr)
        idx = unique(as.data.frame(mcols(full_gr[queryHits(hits)])))
        full$hasIntron = "no"
        idx_intron = match(idx$Row.names, full$Row.names)
        if (length(idx_intron) != 0)
            full[idx_intron, ]$hasIntron = "yes"
    } else (full$hasIntron = NA)
    full$type.toPlot = ifelse(full$hasIntron == "yes" & full$type.50Overlap == "protein_coding", "protein_coding_intron",
        full$type.50Overlap)
    ## grouping into 11 biotypes ##########
    full$type.toPlot = sapply(full$type.toPlot, function(x) ifelse(grepl("protein_coding", x) & grepl("antisense",
        x), "protein_coding_antisense", x))

    full$type.toPlot = sapply(full$type.toPlot, function(x) ifelse(grepl("protein_coding,", x), "protein_coding_mixed",
        x))
    full$type.toPlot = sapply(full$type.toPlot, function(x) ifelse(grepl(", protein_coding", x), "protein_coding_mixed",
        x))
    full$type.toPlot = sapply(full$type.toPlot, function(x) ifelse(grepl("lincRNA", x), "lincRNA",
        x))
    full$type.toPlot = sapply(full$type.toPlot, function(x) ifelse(grepl("antisense,", x), "antisense",
        x))
    full$type.toPlot = sapply(full$type.toPlot, function(x) ifelse(grepl(", antisense", x), "antisense",
        x))

    label = c("protein_coding", "protein_coding_mixed", "lincRNA", "antisense", "pseudogene, processed_pseudogene",
        "pseudogene, unprocessed_pseudogene", "de novo", "protein_coding_antisense", "protein_coding_intron",
        "miRNA")
    full$type.toPlot = sapply(full$type.toPlot, function(x) ifelse(!x %in% label, "other_noncoding",
        x))
    return(full)
}
