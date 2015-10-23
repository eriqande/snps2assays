#' Example of target SNPs data frame
#'
#' A data frame that corresponds to the VCF and FASTA files
#' in "textdata" in the package.
#'
#' @format A tbl_df data frame with 25 rows and 3 variables:
#' \describe{
#'   \item{CHROM}{the contig upon which the SNP is found,}
#'   \item{POS}{The position, starting from one, of the SNP in the contig,}
#'   \item{FST}{Sham Fst values for these SNPs. These are just in there to show that you can
#'   put extra columns of information you into your "targets" data frame.}
#' }
#' @source Thanks to Kristen Ruegg.
"example_target_snps"


# I can't document these like this, so I have commented them out

# #' Example VCF file in package directory "textdata"
# #'
# #' A small gzipped VCF file for examples.  Appears in package directory "textdata"
# #'
# #' @source Thanks to Kristen Ruegg.
# "vcf.txt.gz"
#
#
# #' Example FASTA file in package directory "textdata"
# #'
# #' A small gzipped FASTA file for examples.  Appears in package directory "textdata"
# #'
# #' @source Thanks to Kristen Ruegg.
# "vcf.txt.gz"

