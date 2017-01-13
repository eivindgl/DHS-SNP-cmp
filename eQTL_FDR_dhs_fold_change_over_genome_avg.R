pacman::p_load(
  tidyverse,
  stringr,
  forcats,
  VariantAnnotation,
  biomaRt
)
if (!file.exists('input_data/gene_list.tsv')) {
  ensembl_version = 'feb2014.archive.ensembl.org'
  mart <- useMart(
    'ENSEMBL_MART_ENSEMBL',
    host = ensembl_version,
    dataset = 'hsapiens_gene_ensembl')

    getBM(
      attributes = c('chromosome_name', 'start_position', 'end_position',
                     'ensembl_gene_id', 'strand', 'gene_biotype'),
    #filters = 'ensembl_gene_id',
    #values = ensembl_ids,
    mart = mart) %>%
    as_tibble %>%
    dplyr::select(
      chrom = chromosome_name,
      start = start_position,
      end = end_position,
      gene = ensembl_gene_id,
      strand,
      gene_biotype) %>%
    filter(chrom %in% c(as.character(1:22), 'X', 'Y')) %>%
    mutate(
      chrom = paste('chr', chrom, sep = ''),
      strand = ifelse(strand > 0, '+', '-')) %>%
    write_tsv('input_data/gene_list.tsv')
}
all_genes_df <- read_tsv('input_data/gene_list.tsv')
all_genes <- makeGRangesFromDataFrame(
  all_genes_df %>%
    filter(gene_biotype == 'protein_coding'))
all_snps <- readVcf('input_data/all.vcf', genome = 'GRCh37') %>% rowRanges()
# change from NCBI to UCSC style (e.g. 1 -> chr1)
seqlevelsStyle(all_snps) <- 'UCSC'

# test: using subset of SNPs for background
TSS_window <- promoters(all_genes, upstream = 1e5, downstream = 1e5)
TSS_window_snps <- !is.na(findOverlaps(all_snps, TSS_window, select = 'first'))
win_snps <- all_snps[TSS_window_snps]
cat('Selecting', length(win_snps), 'SNPs out of a possible', length(all_snps))

eqtl <- read_tsv('output_data/rasqual_all_sig_FDR.tsv') %>%
  dplyr::rename(SNP = snps) %>%
  mutate(pcat = cut_number(df$pvalue, n = 3))

eqtl %>%
  count(pcat)
dhs_snp_map <- read_csv('output_data/all_snp_dhs_map.csv') %>%
  filter(SNP %in% names(win_snps))

snp_glob_summary <- dhs_snp_map %>%
  group_by(DHS) %>%
  summarise(
    tot_snps = length(win_snps),
    glob_dhs_snps = n(),
    glob_frac = glob_dhs_snps / tot_snps)

eqtl_summary <- eqtl %>%
  inner_join(dhs_snp_map) %>%
  group_by(pcat, DHS) %>%
  summarise(
    top_dhs_snps = n()
  ) %>% ungroup()

snp_summary <- eqtl %>%
  group_by(pcat) %>%
  summarise(SNPs_in_pcat = n()) %>%
  inner_join(eqtl_summary, by = 'pcat') %>%
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
  group_by(group) %>%
  summarise(median_fchange = median(fold_change)) %>%
  inner_join(snp_df)

snp_df %>%
  write_csv('output_data/Maurano_2012_fig5_inspired_eqtl.csv')

snp_df %>%
  mutate(
    group = fct_reorder(group, median_fchange)
  ) %>%
  ggplot(aes(group, fold_change, color = pcat)) +
  geom_boxplot() +
  scale_color_discrete(name = 'GWAS SNP P-value intervals') +
  labs(
    title = 'SNPs in DHS sites per cell type for 1309 CeD eQTLs',
    subtitle = 'Fold change over genome average for cell type. Grouped by P-value intervals.',
    x = '',
    y = 'Fold Change',
    caption = 'P-values from gsTCC eQTL study.\nInspired by figure 5 in Maurano et. al 2012'
  )

snp_df %>%
  mutate(
    group = fct_reorder(group, median_fchange),
    fold_change = log(fold_change)

  ) %>%
  ggplot(aes(group, fold_change, color = group)) +
  geom_boxplot() +
  scale_color_discrete(name = 'GWAS SNP P-value intervals') +
  labs(
    title = 'SNPs in DHS sites per cell type for 1309 CeD eQTLs',
    subtitle = 'Fold change over genome average for cell type. Grouped by P-value intervals.',
    x = '',
    y = 'Log Fold Change',
    caption = 'P-values from gsTCC eQTL study.\nInspired by figure 5 in Maurano et. al 2012'
  )
