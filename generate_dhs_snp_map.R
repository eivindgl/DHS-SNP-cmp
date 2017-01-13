pacman::p_load(
  tidyverse,
  stringr,
  purrr,
  assertthat,
  GenomicRanges,
  VariantAnnotation,
  rtracklayer,
  biomaRt
  )

vcf <- readVcf('input_data/all.vcf', genome = 'GRCh37')
snp <- rowRanges(vcf)
# change from NCBI to UCSC style (e.g. 1 -> chr1)
seqlevelsStyle(snp) <- 'UCSC'
# code below expects SNP name to be in mcols$name
mcols(snp)$name <- names(snp)

ensembl_version = 'feb2014.archive.ensembl.org'
snp_mart <- useMart(
  'ENSEMBL_MART_SNP',
  host = ensembl_version,
  dataset = 'hsapiens_snp')

ced_top <- list(path = 'input_data/dubois_2010_top1K_SNPs.tsv')
ced_top$orig <- read_tsv(ced_top$path)
ced_top$bm <- getBM(
    attributes = c('refsnp_id', 'chr_name', 'chrom_start'),
    filters = 'snp_filter',
    values = ced_top$orig$SNP,
    mart = snp_mart) %>%
  as_tibble()

ced_top$bm_clean <- ced_top$bm %>%
  dplyr::rename(name = refsnp_id) %>%
  group_by(name) %>%
  summarise(
    chrom = first(chr_name),
    start = first(chrom_start)) %>%
  mutate(end = start + 1) %>%
  dplyr::select(chrom, start, end, name)

ced_top$gr <- makeGRangesFromDataFrame(ced_top$bm_clean, keep.extra.columns = TRUE)
seqlevelsStyle(ced_top$gr) <- 'UCSC'

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
  not_NA <- compose(`!`, is.na)
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
# Create mapping between Top 1000 SNPs from Dubois 2010 CeD GWAS paper and DHS sites
#
names(dhs) %>%
  map(~ dhs_snp_df(ced_top$gr, .x, dhs[[.x]])) %>%
  bind_rows() %>%
  write_csv('output_data/ced_dubois_snp_dhs_map.csv')


#
# Create Mapping between all immunochip SNPs and DHS sites
#
names(dhs) %>%
  map(~ dhs_snp_df(snp, .x, dhs[[.x]])) %>%
  bind_rows() %>%
  write_csv('output_data/all_snp_dhs_map.csv')

find_dhs_snps(ced_snps, dhs[[1]])

#
# Create Mapping between eQTL SNPs and DHS sites within 1MB
#
eqtl_snp_names <- unique(read_tsv('output_data/rasqual_fulll_ow_pvalue.tsv')$snps)
eqtl_snps <- snp[eqtl_snp_names]
eqtl_snps_1Mb <- flank(eqtl_snps, 5e5, both = TRUE)

find_dhs_snps <- function(snp, dhs) {
}

f <- function(snp, dhs_name, dhs) {
  assert_that(snp %>% mcols %has_name% 'name')
  hits <- findOverlaps(snp, dhs)
  snps <- snp[from(hits)]
  d <- dhs[to(hits)]
  tibble(
    SNP = snps$name,
    DHS = dhs_name,
    bp = width(d)
  )
}


x <- names(dhs) %>%
  map(~ f(eqtl_snps, .x, dhs[[.x]])) %>%
  bind_rows() %>%
  write_tsv('output_data/eqtl_snp_1Mb_dhs_map.tsv')
