pacman::p_load(
  tidyverse,
  rtracklayer
)

dir.create('tmp', showWarnings = FALSE)
ced_regions <- import('output_data/hg19_Immunobase_Celiac_Disease.bed')
eqtl <- read_tsv('output_data/rasqual_low_pvalue.tsv') %>%
  dplyr::rename(SNP = snps)
snp_dhs <- read_csv('output_data/all_snp_dhs_map.csv')

#
# Prints number of SNPs per DHS type (should be normalized by coverage?)
#
dhs_snp_stats <- eqtl %>%
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
  select(DHS, eqtl_SNPs, eqtl_SNPs_per_mb) %>%
  arrange(desc(eqtl_SNPs_per_mb)) %>%
  print(n = Inf)
