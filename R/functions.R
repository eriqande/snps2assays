

#    \item{LENGTH}{The length (number of bases) of the DNA segment referred to in CHROM. This is not a standard part of
#      most VCF files, so you might have to put this on yourself. }



#' add summaries about distance to flanking snps to assay rows
#'
#' More later
#' @param V a data frame of all the variants detected in sequences from individuals.
#' Typically this will be from a VCF file (maybe just the first few columns). It can
#' have any number of columns, but it must have the following:
#' \describe{
#'    \item{CHROM}{The name of the contiguous piece of DNA that the variant is found in.
#'      This will typically be a RAD locus identifier, etc.}
#'    \item{POS}{The position of this particular variant within the CHROM DNA.  Starts
#'      counting from 1.}
#'    \item{REF}{The nucleotide base of the reference sequence at each of the positions.  Must be uppercase A, C, G, or T}
#'    \item{ALT}{The alternate base(s) at each of the positions.  Must be uppercase, A, C, G, or T, or if more than 1
#'    alternate base is present, they must be listed, and comma-separated, and there can be no more than 3 alternate bases.}
#' }
#' Any columns other than \code{CHROM}, \code{POS}, \code{REF} and \code{ALT} will be ignored.
#' If consSeq is NA then V must also have a LENGTH column that gives the number of bases in the
#' contig.
#' @param targets A data frame of the variants that are targets for assay development.  This
#' can also have whatever columns that are desired, but it must have CHROM and POS, concordant with \code{V},
#' and it cannot have a column LENGTH in it.
#'  The extra columns in this
#'  data frame will all be represented in the output.
#' @param consSeq A data frame of consensus sequences.  Must have one column \code{CHROM}, exactly concordant with the
#' naming conventions of \code{V} and \code{targets}, and one column \code{SEQS} which are the consensus sequences as
#' strings.  There can be more columns, but they will be ignored.  If this is NA, then the function returns just the
#' flanking information summaries, without the sequence from which to build assays.
#' @param reqDist The required minimum number of bases between a target SNP and the nearest flanking SNP
#' for a target to be designable (according to whomever's criterion).
#' @return Returns a data frame with rows as in \code{targets} with additional columns appended to
#' it, namely:
#' \itemize{
#'    \item One here
#'    \item One there
#' }
#' @export
assayize <- function(V, targets,  consSeq = NA, reqDist = 20) {

  if(!("CHROM" %in% names(V))) stop("No column CHROM in V")
  if(!("POS" %in% names(V))) stop("No column POS in V")
  if(is.na(consSeq) && !("LENGTH" %in% names(V))) stop("No column LENGTH in V while consSeq is NA")
  if(!("REF" %in% names(V))) stop("No column REF in V")
  if(!("ALT" %in% names(V))) stop("No column ALT in V")
  if(!("CHROM" %in% names(targets))) stop("No column CHROM in targets")
  if(!("POS" %in% names(targets))) stop("No column POS in targets")
  if(("LENGTH" %in% names(targets))) stop("targets has a column called LENGTH, but it should not. That should be in V")

  # check here to make sure that LENGTH and POS are numeric and coerce REF and ALT to as.character
  # in case they are already factors.
  if(is.na(consSeq) && !is.numeric(V$LENGTH)) stop("consSeq is NA so V has column LENGTH, but LENGTH is not numeric!")
  if(!is.numeric(V$POS)) stop("POS not numeric in V")
  if(!is.numeric(targets$POS)) stop("POS not numeric in targets")
  if(!is.character(V$REF)) {
    V$REF <- as.character(V$REF)
    message("Coercing V$REF to character")
  }
  if(!is.character(V$ALT)) {
    V$ALT <- as.character(V$ALT)
    message("Coercing V$ALT to character")
  }

  # if LENGTH REF or ALT are in targets, we need to remove them from there along with a message to the user
  # saying that values from V will be used.
  if(any(names(targets) %in% c("LENGTH", "REF", "ALT"))) {
    message("Removing LENGTH, REF, and/or ALT from parameter \"targets\".  Values from V and consSeq will be used.")
  }


  # check to make sure REF and ALT are strings and do not include any [^ACGT,]
  if(any(stringr::str_detect(V$REF, "[^ACGTU]"))) stop("Detected something other than ACGT or U in the REF column in V")
  if(any(stringr::str_detect(V$ALT, "[^ACGTU,]"))) stop("Detected something other than ACGTU or \",\" in the ALT column in V")

  # if LENGTH is missing but consSeq is not, then I could certainly get LENGTH from that
  # file....

  # make a named vector of IUPAC codes here:
  iupac <- c(`A` = "A",
             `C` = "C",
             `G` = "G",
             `T` = "T",
             `U` = "U",
             `AG` = "R",
             `CT` = "Y",
             `CG` = "S",
             `AT` = "W",
             `GT` = "K",
             `AC` = "M",
             `CGT` = "B",
             `AGT` = "D",
             `ACT` = "H",
             `ACG` = "V",
             `N` = "N",
             `.` = ".",
             `-` = "-"
  )
  # a few functions that will be useful within this function
  left_flank <- function(pos) {
    lefties <- c(0, pos[-length(pos)])
    pos - lefties - 1
  }
  right_flank <- function(pos, length) {
    righties <- c(pos[-1], length[1])
    righties - pos - 1
  }
  bases_present <- function(ref, alt) { # return a comma-separated string of bases present, sorted in alphabetical order
    a <- stringr::str_split(ref, ",")
    b <- stringr::str_split(alt, ",")
    lapply(1:length(a), function(x) { paste(sort(c(a[[x]], b[[x]])), collapse = "")}) %>%
      unlist
  }

  chrom2get <- unique(targets$CHROM)

  # get the left and right flank summaries
  summaries <- V %>%
    filter(CHROM %in% chrom2get) %>%
    arrange(CHROM, POS)  %>%  # this is critical!
    select(CHROM, POS, LENGTH, REF, ALT) %>%
    group_by(CHROM) %>%
    mutate(LeftFlank = left_flank(POS),
           RightFlank = right_flank(POS, LENGTH),
           Designable = LeftFlank >= reqDist & RightFlank >= reqDist) %>%
    ungroup %>%
    mutate(BasesPresent = bases_present(REF, ALT),
           IUPAC = iupac[as.character(BasesPresent)])

  if(is.na(consSeq)) {
    return(summaries)
  }

  # down here we will attach a column with strings showing where variation is (with iupac codes) and with
  # the SNP that is to be assayed as, for example, "[A/G]"
  # first, we just put iupac codes in everywhere...


  # and then we filter out just the SNPs that are to have assays designed, and we
  # put their designations into the sequences.

}



#' read in a fasta file and turn into a tbl_df data frame
#'
#' This reads in a fasta file and removes the leading ">" from the name of the contig
#' so that it will match stuff in a corresponding VCF or other variant file.
#' It then returns a tbl_df data frame with columns CHROM and consSeq.  Currently this
#' only woks on files that have two line for each sequence, not with files that break the
#' sequence over multiple lines.
#' @param path  The path to the fasta file. If it is gzipped, it will automatically
#' be gunzipped on the fly.
#' @export
#' @examples
#' grab_fasta(system.file("textdata", "fasta.txt.gz", package = "snps2assays"))
grab_fasta <- function(path) {
  con <- gzfile(path)

  ret <- readLines(con) %>%
    matrix(ncol = 2, byrow = TRUE) %>%
    as.data.frame(stringsAsFactors = FALSE) %>%
    setNames(c("CHROM", "Seq")) %>%
    tbl_df %>%
    mutate(CHROM = stringr::str_replace_all(CHROM, "^>", ""))

  close(con)
  ret
}


#' read in a (possibly gzipped) VCF file and return it as a data frame
#'
#' This scans a VCF file to find the header line (by searching for a line
#' starting with "#CHROM"), and then it reads it in and returns it as a
#' data frame.  It doesn' worry about the formatting of the genotype columns
#' because we are typically concerned here just with getting CHROM, POS,
#' REF, and ALT.  In fact, by default, those four columns are all that
#' this function returns
#' @param path The path to a VCF file (which may or may not be gzipped)
#' @param all_cols Logical indicating whether to return all columns (TRUE)
#' or only CHROM, POS, REF, and ALT.
#' @param top_test How many lines to scan from the top of the file to look
#' for the "#CHROM" showing that it is the header.
#' @export
#' @examples
#' grab_vcf(system.file("textdata", "vcf.txt.gz", package = "snps2assays"))
grab_vcf <- function(path, all_cols = FALSE, top_test = 200) {
  con <- gzfile(path)

  # first scan for the position of the header line
  tmp <- readLines(con, n = top_test) %>%
    stringr::str_detect("^#CHROM")
  #close(con)

  if(sum(tmp) == 0) stop("Didn't find pattern \"^#CHROM\" in the first ", top_test, " lines of file ", path)
  if(sum(tmp) > 1) stop("Found more than one line matching pattern \"^#CHROM\" in the first ", top_test, " lines of file ", path)

  head_line <- which(tmp)

  ret <- read.table(con,
                    comment = "",
                    skip = head_line - 1,
                    header = TRUE,
                    stringsAsFactors = FALSE) %>%
    tbl_df

  names(ret)[1] <- "CHROM"  # deal with the comment character at the beginning there


  if(all_cols == FALSE) {
    ret <- ret %>%
      select(CHROM, POS, REF, ALT)
  }

  # sort it as it should be
  ret %>% arrange(CHROM, POS)
}
