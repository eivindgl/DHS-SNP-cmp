pacman::p_load(
  tidyverse,
  stringr,
  forcats,
  VariantAnnotation,
  biomaRt
)

all_snps <- readVcf('input_data/all.vcf', genome = 'GRCh37') %>% rowRanges()
# change from NCBI to UCSC style (e.g. 1 -> chr1)
seqlevelsStyle(all_snps) <- 'UCSC'

ced_top <- read_tsv('input_data/dubois_2010_top1K_SNPs.tsv')

ced_top <- ced_top %>%  mutate(
  logp = -log10(P_GWAS),
  pcat = cut(P_GWAS, breaks = c(0, 1e-5, 1e-4, 2.5e-4, 5e-4, 1e-3), labels = c(
  '<1e-5', '1e-5<1e-4', '1e-4<2.5e-4', '2.5e-4<5e-4', '5e-4<1e-3')),
  pcat = fct_reorder(pcat, logp))

ced_top %>% count(pcat)

# ced_top <- ced_top %>%  mutate(
#     logp = -log10(P_GWAS),
#     pcat = cut_number(logp, n = 3))



dhs_snp_map <- read_csv('output_data/all_snp_dhs_map.csv')
dhs_snp_map <- read_csv('output_data/ced_dubois_snp_dhs_map.csv') %>%
  bind_rows(dhs_snp_map) %>%
  dplyr::distinct()


snp_glob_summary <- dhs_snp_map %>%
  group_by(DHS) %>%
  summarise(
    tot_snps = length(all_snps),
    glob_dhs_snps = n(),
    glob_frac = glob_dhs_snps / tot_snps)

ced_top_summary <- ced_top %>%
  inner_join(dhs_snp_map) %>%
  group_by(pcat, DHS) %>%
  summarise(
    top_dhs_snps = n()
  ) %>% ungroup()

snp_summary <- ced_top %>%
  group_by(pcat) %>%
  summarise(SNPs_in_pcat = n()) %>%
  inner_join(ced_top_summary, by = 'pcat') %>%
  mutate(pcat_freq = top_dhs_snps / SNPs_in_pcat) %>%
  inner_join(snp_glob_summary, by = 'DHS') %>%
  mutate(fold_change = pcat_freq / glob_frac)

tmp <- snp_summary %>%
  filter(!str_detect(DHS, '^TCC')) %>%
  mutate(group = str_replace(DHS, '_\\d$', ''))

snp_summary <- snp_summary %>%
  filter(str_detect(DHS, '^TCC')) %>%
  mutate(group = paste('TCC', str_extract(DHS, '\\d+$'), sep = '-')) %>%
  bind_rows(tmp) %>%
  mutate(group = fct_relevel(group,
                             'TCC-0', 'TCC-10', 'TCC-30', 'TCC-180',
                             'T-naive', 'Th-1', 'Th-2', 'Th-17', 'T-reg', 'Jurkat'),
         broad_group = ifelse(str_detect(group, '^TCC'), 'gsTCC', 'Encode'))



snp_summary %>%
  ggplot(aes(fold_change, fill = broad_group)) +
  geom_density(alpha = 0.75)

snp_df <- snp_summary %>%
  mutate(
    group = ifelse(str_detect(group, 'TCC'), 'gsTCC', as.character(group)))

snp_df <- snp_df %>%
  filter(pcat == '<1e-5') %>%
  group_by(group) %>%
  summarise(median_fchange = median(fold_change)) %>%
  inner_join(snp_df)

snp_df %>%
  write_csv('output_data/Maurano_2012_fig5_inspired_data.csv')

snp_df %>%
  mutate(
    group = fct_reorder(group, median_fchange)
  ) %>%
  ggplot(aes(group, fold_change, color = pcat)) +
  geom_boxplot() +
  scale_color_discrete(name = 'GWAS SNP P-value intervals') +
  labs(
    title = 'SNPs in DHS sites per cell type for top 1000 CeD GWAS SNPs ',
    subtitle = 'Fold change over genome average for cell type. Grouped by P-value intervals.',
    x = '',
    y = 'Fold Change',
    caption = 'P-values from Dubois et al 2010.\nInspired by figure 5 in Maurano et. al 2012'
  )
