pacman::p_load(
  tidyverse,
  stringr,
  purrr,
  assertthat,
  GenomicRanges,
  VariantAnnotation,
  rtracklayer
  )

vcf <- readVcf('input_data/all.vcf', genome = 'GRCh37')
snp <- rowRanges(vcf)
# change from NCBI to UCSC style (e.g. 1 -> chr1)
seqlevelsStyle(snp) <- 'UCSC'
# code below expects SNP name to be in mcols$name
mcols(snp)$name <- names(snp)

read_bedfiles <- function(bedfiles) {
  xs <- list()
  for (x in bedfiles) {
    name <- tools::file_path_sans_ext(basename(x))
    xs[[name]] <- import(x)
  }
  xs
}


# Prepare data for further processing in other scripts.
# Create a mapping of SNPs overlapping DHS sites from various files.
#
find_dhs_snps <- function(snp, dhs) {
  assert_that(snp %>% mcols %has_name% 'name')
  not_NA <- function(x) !is.na(x)
  is_overlapping <- findOverlaps(snp, dhs, select = 'first') %>% not_NA
  snp[is_overlapping]$name
}

dhs_snp_df <- function(snp, dhs_name, dhs) {
  snps <- find_dhs_snps(snp, dhs)
  tibble(
    SNP = snps,
    DHS = dhs_name
  )
}

#
# Read DHS sites
#
bluebed_paths <- list.files('output_data/bluebed_shortname', full.names = TRUE, pattern = 'bed$')
gsTCC_paths <- list.files('input_data/gsTCC_dhs', full.names = TRUE, pattern = 'bed$') %>%
  purrr::discard(~ str_detect(.x, 'merged'))
dhs <- read_bedfiles(c(bluebed_paths, gsTCC_paths))

#
# SNPs in strong LD with a CeD tag SNP
#
ced_snps <- import('input_data/LD_SNPs/CeD_LD_SNPs_Iris_maf-.001_r2-9.bed')

#
# Create Mapping between CeD LD SNPs and DHS sites
#
names(dhs) %>%
  map(~ dhs_snp_df(ced_snps, .x, dhs[[.x]])) %>%
  bind_rows() %>%
  write_csv('output_data/ced_snp_dhs_map.csv')

#
# Create Mapping between all immunochip SNPs and DHS sites
#
names(dhs) %>%
  map(~ dhs_snp_df(snp, .x, dhs[[.x]])) %>%
  bind_rows() %>%
  write_csv('output_data/all_snp_dhs_map.csv')

find_dhs_snps(ced_snps, dhs[[1]])
