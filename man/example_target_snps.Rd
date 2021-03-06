% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{example_target_snps}
\alias{example_target_snps}
\title{Example of target SNPs data frame}
\format{
A tibble data frame with 25 rows and 3 variables:
\describe{
  \item{CHROM}{the contig upon which the SNP is found,}
  \item{POS}{The position, starting from one, of the SNP in the contig,}
  \item{FST}{Sham Fst values for these SNPs. These are just in there to show that you can
  put extra columns of information you into your "targets" data frame.}
}
}
\source{
Thanks to Kristen Ruegg.
}
\usage{
example_target_snps
}
\description{
A data frame that corresponds to the VCF and FASTA files
in "textdata" in the package.
}
\keyword{datasets}
