#' snps2assays: a package for preparing SNPtype-like assay orders from ddRAD or RAD locus data
#'
#' This package provides one main function, \code{\link{assayize}}, that takes, as
#' input, data frames of:
#' \enumerate{
#'  \item all the known variation in different (typically short) contigs,
#'  \item the SNPs you want to turn into assays,
#'  \item the sequences of all the contigs (or at least the ones that you want to turn into
#'  assays).
#'}
#' It also provides two other functions, \code{\link{grab_vcf}} and
#' \code{\link{grab_fasta}} for reading (possibly gzipped) VCF and fasta
#' files into data frames appropriate for use in \code{\link{assayize}}.
#'
#' @docType package
#' @name snps2assays
NULL
#> NULL
