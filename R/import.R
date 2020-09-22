#' @importFrom dplyr arrange everything filter group_by mutate left_join rename right_join select ungroup
#' @importFrom magrittr %>%
#' @importFrom stats setNames
NULL


# quiets concerns of R CMD check re: the . and other column names
# that appear in dplyr chains
if (getRversion() >= "2.15.1")  {
  utils::globalVariables(
    c(
      ".",
      "ALT",
      "BasesPresent",
      "CHROM",
      "Designable",
      "GC_content",
      "IUPAC",
      "LENGTH",
      "LeftFlank",
      "POS",
      "REF",
      "RightFlank",
      "Seq"
    )
  )
}
