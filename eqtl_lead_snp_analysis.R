#
# Currently full SNP set (i.e. all SNPs gene pairs with p-value < 1e-3)
# This scripts is currently a little messy.
#
# Counts eQTL SNPs per DHS type
# Counts eQTL SNPs per CeD loci (only two; 2p16.1 and 21q22.3 have multiple eQTLs)
pacman::p_load(
  tidyverse,
  rtracklayer,
  VariantAnnotation,
  biomaRt
)

dir.create('tmp', showWarnings = FALSE)
ced_regions <- import('output_data/hg19_Immunobase_Celiac_Disease.bed')
eqtl <- read_tsv('output_data/rasqual_low_pvalue.tsv') %>%
  dplyr::rename(SNP = snps)
feqtl <- read_tsv('output_data/rasqual_full_low_pvalue.tsv') %>%
  dplyr::rename(SNP = snps)

snp_dhs <- read_csv('output_data/all_snp_dhs_map.csv')

#
# Prints number of SNPs per DHS type (should be normalized by coverage?)
#
dhs_snp_stats <- feqtl %>%
  group_by(SNP) %>%
  summarise(num_genes = n()) %>%
  inner_join(snp_dhs) %>%
  group_by(DHS) %>%
  summarise(eqtl_SNPs = n()) %>%
  arrange(desc(eqtl_SNPs))

ced_dhs_summary <- read_csv('output_data/DHS_CeD_LD_summary.csv') %>%
  dplyr::rename(DHS = name)

df <- inner_join(dhs_snp_stats, ced_dhs_summary)
df %>%
  mutate(eqtl_SNPs_per_mb = eqtl_SNPs / (tot_cov / 1000.) ) %>%
  dplyr::select(DHS, eqtl_SNPs, eqtl_SNPs_per_mb) %>%
  arrange(desc(eqtl_SNPs_per_mb)) %>%
  print(n = Inf)

vcf <- readVcf('input_data/all.vcf', genome = 'GRCh37')
snp <- rowRanges(vcf)[unique(feqtl$SNP)]
seqlevelsStyle(snp) <- 'UCSC'
export(snp, '~/Downloads/rasqual_full_snps.bed')

idx <- findOverlaps(snp, ced_regions, select = 'first')
idx <- which(!is.na(idx))
ced_lead_snps <- names(snp[idx]) # SNPs in immunobase regions
feqtl %>%
  filter(SNP %in% ced_lead_snps)

ensembl_version = 'feb2014.archive.ensembl.org'
mart <- useMart(
  'ENSEMBL_MART_ENSEMBL',
  host = ensembl_version,
  dataset = 'hsapiens_gene_ensembl')

hgnc <- getBM(
  attributes = c('chromosome_name', 'start_position', 'end_position', 'external_gene_id'),
  filters = 'ensembl_gene_id',
  values = feqtl$gene,
  mart = mart) %>%
  as_tibble %>%
  mutate(chromosome_name = paste('chr', chromosome_name, sep = ''))

hgnc %>%
  arrange(chromosome_name, start_position) %>%
  dplyr::distinct() %>%
  write_tsv('~/Downloads/feQTL_genes.bed', col_names = FALSE)
