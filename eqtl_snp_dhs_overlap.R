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
eqtl <- read_tsv('output_data/rasqual_all_sig_FDR.tsv') %>%
  dplyr::rename(SNP = snps)
feqtl <- read_tsv('output_data/rasqual_full_all_sig_FDR.tsv') %>%
  dplyr::rename(SNP = snps)

snp_dhs <- read_csv('output_data/all_snp_dhs_map.csv')
ced_snp_dhs <- read_csv('output_data/ced_snp_dhs_map.csv')
eqtl_snp_dhs <- read_tsv('output_data/eqtl_snp_1Mb_dhs_map.tsv')

#
# Prints number of SNPs per DHS type (should be normalized by coverage?)
#
dhs_snp_stats <- feqtl %>%
  group_by(SNP) %>%
  summarise(num_genes = n()) %>%
  inner_join(eqtl_snp_dhs) %>%
  group_by(DHS) %>%
  summarise(eqtl_SNPs = n(), dhs_1mb_cov = sum(bp)) %>%
  arrange(desc(eqtl_SNPs)) %>%
  mutate(eqtl_SNPs_per_kb_dhs = eqtl_SNPs / (dhs_1mb_cov / 1000.))

ced_dhs_summary <- read_csv('output_data/DHS_CeD_LD_summary.csv') %>%
  dplyr::rename(DHS = name)

dhs_snp_stats %>%
  inner_join(ced_dhs_summary) %>%
  ggplot(aes(dhs_1mb_cov, eqtl_SNPs, color = group)) + geom_point() + scale_x_log10()

dhs_snp_stats %>%
  inner_join(ced_dhs_summary) %>%
  dplyr::select(group, DHS, eqtl_SNPs, dhs_1mb_cov, eqtl_SNPs_per_kb_dhs) %>%
  ggplot(aes(group, eqtl_SNPs_per_kb_dhs, color = group)) + geom_boxplot()



vcf <- readVcf('input_data/all.vcf', genome = 'GRCh37')
snp <- rowRanges(vcf)[unique(feqtl$SNP)]
seqlevelsStyle(snp) <- 'UCSC'
export(snp, '~/Downloads/rasqual_full_snps.bed')
