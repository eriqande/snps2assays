

#    \item{LENGTH}{The length (number of bases) of the DNA segment referred to in CHROM. This is not a standard part of
#      most VCF files, so you might have to put this on yourself. }



#' Prepare SNPs for assay orders (calculate SNP flanking distances, GC content, and write out sequence)
#'
#' This function takes, as
#' input, data frames of:
#' \enumerate{
#'  \item all the known variation in different (typically short) contigs,
#'  \item the SNPs you want to turn into assays,
#'  \item the sequences of all the contigs (or at least the ones that you want to turn into
#'  assays).
#'}
#' It returns information needed for screening variation and ordering SNP assays.
#' @param V a data frame of all the variants detected in sequences from individuals.
#' Typically this will be from a VCF file (maybe just the first few columns). It can
#' have any number of columns, but it must have the following:
#' \describe{
#'    \item{CHROM}{The name of the contiguous piece of DNA that the variant is found in.
#'      This will typically be a RAD locus identifier, etc.}
#'    \item{POS}{The position, starting from 1, of this particular variant within the CHROM DNA.}
#'    \item{REF}{The nucleotide base of the reference sequence at each of the positions.  Must be uppercase A, C, G, or T}
#'    \item{ALT}{The alternate base(s) at each of the positions.  Must be uppercase, A, C, G, or T, or, if more than 1
#'    alternate base is present, they must be listed, and comma-separated, and there can be no more
#'    than 3 alternate bases, (or perhaps a "-" or ".", though in those cases you probably don't want
#'    to be developing a SNP assay, anyway.)}
#' }
#' Any columns other than \code{CHROM}, \code{POS}, \code{REF} and \code{ALT} will be ignored.
#' If \code{consSeq} is NULL then V must also have a LENGTH column that gives the number of bases in the
#' contig.
#' @param targets A data frame of the variants that are targets for assay development.  This
#' can also have whatever columns that are desired, but it must have CHROM and POS, concordant with \code{V},
#' and it cannot have a column LENGTH in it.
#'  The extra columns in this
#'  data frame will all be represented in the output.
#' @param consSeq A data frame of consensus sequences.  Must have one column \code{CHROM}, exactly concordant with the
#' naming conventions of \code{V} and \code{targets}, and one column \code{Seq} which holds the consensus sequences as
#' strings. The sequence should consist of capital A, C, G, or T.
#' There can be more columns, but they will be ignored.  If this is NULL, then the function returns just the
#' flanking information summaries, without the sequence from which to build assays. However, in that
#' case, there must be a LENGTH column in V.
#' @param reqDist The required minimum number of bases between a target SNP and the nearest flanking SNP
#' for a target to be designable (according to whomever's criterion).
#' @param reqDistFromEnd The required number of bases between the SNP site and the start or end of
#' the contig for the SNP to be designable as an assay.
#' @param GCmax The maximum GC content of the reference contig allowable for an assay to still be
#' considered designable.
#' @param allVar Logical indicating whether all variable sites within the CHROMs containing the targets
#' should be returned, or just the targets themselves.
#' @return Returns a data frame with all the columns in \code{targets} along with additional columns
#' appended to it. Namely, you can expect the following columns.
#' \describe{
#'  \item{CHROM}{The contig identifier}
#'  \item{POS}{Position of the focal SNP within the contig}
#'  \item{LENGTH}{Length of the contig.}
#'  \item{REF}{The reference base at the site}
#'  \item{ALT}{The alternative base(s) at the site}
#'  \item{LeftFlank}{Number of SNP-free bases to the "left" of the focal SNP}
#'  \item{RightFlank}{Number of SNP-free bases to the "right" of the focal SNP}
#'  \item{GC_content}{Fraction of sites in the reference sequence that are G's or C's.  Does not appear if
#'  \code{consSeq} is NULL}
#'  \item{Designable}{Logical indicating whether this focal SNP is designable given the
#'  criteria specified in \code{reqDist}, \code{reqDistFromEnd}, and, if \code{consSeq} is not
#'  NULL, also given \code{GCmax}.}
#'  \item{BasesPresent}{All the bases (or other variation) observed at the SNP in sorted order}
#'  \item{IUPAC}{The IUPAC code that describes the variation observed at the focal SNP}
#'  \item{Seq}{The sequence of the contig as prepared for assay order.  Variation at non-focal SNPs
#'  is given by IUPAC codes while variation at the focal SNP is specified like so: [T/G]. This
#'  column is not present if \code{consSeq} is NULL.}
#'  \item{Any other columns that were present in \code{targets}}{Whatever extra columns that were
#'  in \code{targets} that were not duplicates of the above colums get stuck on in after CHROM
#'  and POS}

#' }
#'
#' If \code{allVar} is set to FALSE then only rows corresponding to SNPs in \link{targets} will
#' be returned.  If \code{allVar} is TRUE, then all the SNPs in V that fall in the contigs within
#' which the SNPs in \code{targets} fall will be returned (one row for each).  In that case, the
#' "\code{Any other columns that were present in targets}" will be NA.
#' @export
#' @examples
#' # get the vcf file:
#' vcf <- grab_vcf(system.file("textdata", "vcf.txt.gz", package = "snps2assays"))
#'
#' # get the fasta file
#' fasta <- grab_fasta(system.file("textdata", "fasta.txt.gz", package = "snps2assays"))
#'
#' # get our data frame of target SNPs
#' data(example_target_snps)
#'
#' # now, assayize them!
#' assays <- assayize(vcf, example_target_snps, fasta)
#' assays
#'
#' # here we run it without the consensus sequences.  The length of each contig is part
#' # of the RAD locus name (the 6th "_"-separated field), so we can use that:
#' vcf2 <- vcf %>%
#'  mutate(LENGTH = stringr::str_split(CHROM, "_") %>%
#'    lapply("[", 6) %>%
#'    unlist %>%
#'    as.numeric # Note that if LENGTH is not numeric assayize will throw an error!
#'    )
#'
#' # now that vcf2 has a LENGTH column, we can assayize it
#' # and for fun, let's require more distance around the SNP
#' # and only return the target SNPs
#' assays2 <- assayize(vcf2, example_target_snps, reqDist = 25, allVar = FALSE)
#' assays2
assayize <- function(V,
                     targets,
                     consSeq = NULL,
                     reqDist = 20,
                     reqDistFromEnd = 40,
                     GCmax = 0.65,
                     allVar = TRUE) {

  if(!("CHROM" %in% names(V))) stop("No column CHROM in V")
  if(!("POS" %in% names(V))) stop("No column POS in V")
  if(is.null(consSeq) && !("LENGTH" %in% names(V))) stop("No column LENGTH in V while consSeq is NULL")
  if(!("REF" %in% names(V))) stop("No column REF in V")
  if(!("ALT" %in% names(V))) stop("No column ALT in V")
  if(!("CHROM" %in% names(targets))) stop("No column CHROM in targets")
  if(!("POS" %in% names(targets))) stop("No column POS in targets")
  if(("LENGTH" %in% names(targets))) {
    message("targets has a column called LENGTH, which will be ignored. LENGTH is taken from the sequences in consSeq or should appear in V if consSeq is NULL")
    targets <- targets %>% select(-LENGTH)
  }

  # check here to make sure that LENGTH and POS are numeric and coerce REF and ALT to as.character
  # in case they are already factors.
  if(is.null(consSeq) && !is.numeric(V$LENGTH)) stop("consSeq is NULL so V has column LENGTH, but LENGTH is not numeric!")
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
  if(any(names(targets) %in% c("REF", "ALT"))) {
    message("Removing REF, and/or ALT from parameter \"targets\".  Values from V and consSeq will be used.")
  }

  # SHOULD CHECK TO MAKE SURE THERE ARE NO DUPLICATED CHROM+POSs

  # check to make sure REF and ALT are strings and do not include any [^ACGT,]
  if(any(stringr::str_detect(V$REF, "[^ACGTU.]"))) stop("Detected something other than ACGTU or . in the REF column in V")
  if(any(stringr::str_detect(V$ALT, "[^ACGTU,.-]"))) stop("Detected something other than ACGTU or , or . or - in the ALT column in V")

  # if we have consSeq then we add a column LENGTH to it, and we also add that to V
  if(!is.null(consSeq)) {
    consSeq <- consSeq %>%
      mutate(LENGTH = nchar(Seq)) %>%
      select(CHROM, Seq, LENGTH)  # strip it down to only these

    if("LENGTH" %in% names(V)) {
      message("Note: V has a column called LENGTH, but consSeq is not NULL.  We will take LENGTH from consSeq")
       V <- V %>%
        rename(LENGTH.V = LENGTH)
    }

    # here we add LENGTH from consSeq to V
    V <- V %>%
      left_join(., consSeq %>% select(CHROM, LENGTH))
  }

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
             `ACGT` = "N",
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
           Designable = (LeftFlank >= reqDist) &
             (RightFlank >= reqDist) &
             (POS > reqDistFromEnd) &
             (LENGTH - POS >= reqDistFromEnd)) %>%
    ungroup %>%
    mutate(BasesPresent = bases_present(REF, ALT),
           IUPAC = iupac[as.character(BasesPresent)])




  # if we *do* have those consensus sequences, this is the part we are all waiting for!
  # We put the IUPAC codes and
  # snp design specifications in the actual sequences (as long as we have them---i.e.
  # consSeq is not NULL)
  if(!is.null(consSeq)) {
    # quick functions.  seq is a sequence, pos is a vector of positions within that sequence
    # and iup is a vector of iupac codes for those positions at those positions.
    gc_content <- function(seq) { # in this case, seq is a whole vector of sequencies
      ret <- stringr::str_split(seq, "") %>%
        lapply(function(x) mean(x %in% c("G", "C"))) %>%
        unlist
      ret
    }
    insert_iupac <- function(seq, pos, iup) {
      tmp <- stringr::str_split(seq, "")[[1]]
      tmp[pos] <- iup
      paste(tmp, collapse = "")
    }
    insert_target <- function(seq, pos, ref, alt) {
      tmp <- stringr::str_split(seq, "")[[1]]
      tmp[pos] <- paste("[", ref[1], "/", alt[1], "]", sep = "")  # shouldn't need the [1]'s in there but put them in just in case...
      paste(tmp, collapse = "")
    }

    # down here we will attach a column with strings showing where variation is (with iupac codes) and with
    # the SNP that is to be assayed as, for example, "[A/G]"
    # first, we just put iupac codes in everywhere...
    summaries <- summaries %>%
      left_join(., consSeq %>% select(CHROM, Seq)) %>%
      mutate(GC_content = gc_content(Seq)) %>%
      group_by(CHROM) %>%
      mutate(Seq = insert_iupac(Seq, POS, IUPAC)) %>%
      group_by(CHROM, POS) %>%
      mutate(Seq = insert_target(Seq, POS, REF, ALT)) %>%
      ungroup %>%
      mutate(Designable = Designable & (GC_content <= GCmax))
  }



  # prep to return just the stuff that is desired
  if(allVar == TRUE) {
    ret <-left_join(summaries, targets)
  } else {
    ret <- right_join(summaries, targets)
  }

  if(!is.null(consSeq)) {
    ret <- ret %>%
      select(CHROM, POS, LENGTH, REF, ALT, LeftFlank, RightFlank, GC_content, Designable, BasesPresent, IUPAC, Seq, everything())
  } else {
    ret <- ret %>%
      select(CHROM, POS, LENGTH, REF, ALT, LeftFlank, RightFlank, Designable, BasesPresent, IUPAC, everything())
  }
  ret
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
