
#### make some example VCF data: ####
allSnps <- read.table("/Volumes/PerseDrive3/SWTH_WIWA_BIOINFORMATICS/Berkeley_KRuegg_201110/PopGen/UCLA_birds_popgen_075_w55_k75_r75_BC03N007vsSWTH.txt",
                      comment = "",
                      skip = 17,
                      header = TRUE,
                      stringsAsFactors = FALSE)

allSnps <-  allSnps %>% tbl_df
tmp <- allSnps %>%
  arrange(X.CHROM, POS)

vcf <- tmp[1:2000,] %>%
  mutate(CHROM = X.CHROM) %>%
  select(CHROM, everything()) %>%
  select(-X.CHROM)

# now write that stuff out
readLines("/Volumes/PerseDrive3/SWTH_WIWA_BIOINFORMATICS/Berkeley_KRuegg_201110/PopGen/UCLA_birds_popgen_075_w55_k75_r75_BC03N007vsSWTH_annotated_SNPs.txt",
          n = 18) %>%
  cat(., file = "inst/textdata/vcf.txt", sep = "\n")

write.table(vcf,
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE,
            sep = "\t",
            append = "TRUE",
            file = "inst/textdata/vcf.txt")

system("gzip inst/textdata/vcf.txt")


# note, we can read it like this:
readLines(gzfile("inst/textdata/vcf.txt.gz"), n = 18)


#### Take a subset of those and make sham SNPs you want to develop ####
set.seed(5)
example_target_snps <- vcf %>%
  sample_n(25) %>%
  select(CHROM, POS) %>%
  mutate(FST = runif(length(CHROM)))

save(example_target_snps, file = "data/example_target_snps.rda")





#### Now get some fasta output for all those ####
fasta1 <- readLines("/Volumes/PerseDrive3/SWTH_WIWA_BIOINFORMATICS/Berkeley_KRuegg_201110/Genome/UCLA_Birds_SbfI_BC03N0007_v1_filtered_assembly.txt")

fasta2  <- fasta1 %>%
  matrix(ncol = 2, byrow = T)


fastadf <- data.frame(CHROM = stringr::str_replace(fasta2[,1], ">", ""), SEQ = fasta2[,2], stringsAsFactors = FALSE) %>%
  tbl_df

# now get the ones that correspond to vcf
ours <- fastadf %>%
  filter(CHROM %in% unique(vcf$CHROM))

# now, write it out to a fasta file:
paste(">", ours$CHROM, "\n", ours$SEQ, "\n", sep = "") %>%
  cat(sep = "", file = "inst/textdata/fasta.txt")

system("gzip inst/textdata/fasta.txt")


