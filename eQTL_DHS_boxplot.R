pacman::p_load(
  tidyverse,
  forcats,
  stringr,
  magrittr
)

order_by_median_pvalue <- function(df) {
  cdf <- df %>%
    group_by(DHS) %>%
    summarize(medianp = median(pvalue), meanp = mean(pvalue)) %>%
    arrange(medianp)
  df %>% mutate(DHS = factor(DHS, levels = cdf$DHS))
}

gen_boxplot <- function(snp_score, bed_snp_map, timepoint) {
  df <- inner_join(bed_snp_map, snp_score)
  df <- df %>% filter(!is.na(pvalue))

  tmp <- df %>%
    filter(!str_detect(DHS, '^TCC')) %>%
    mutate(group = str_replace(DHS, '_\\d$', ''))

  df <- df %>%
    filter(str_detect(DHS, '^TCC')) %>%
    mutate(group = paste('TCC', str_extract(DHS, '\\d+$'), sep = '-')) %>%
    bind_rows(tmp)

  # Add ALL SNPs for comparison
  df <- snp_score %>%
    mutate(DHS = 'ALL', group = 'ALL') %>%
    filter(!is.na(pvalue)) %>%
    bind_rows(df)
  #
  df %<>% order_by_median_pvalue

  p <- df %>%
    ggplot(aes(DHS, pvalue, color = group)) +
    geom_boxplot() +
    scale_y_continuous(limits = c(0, 0.25)) +
    coord_flip()
  p +
    ggtitle('Boxplot of eQTL p-values for SNPs in DHS sites',
            subtitle = timepoint)
}

bed_snp_map <- read_csv('output_data/snp_dhs_map.csv')

for (score_path in list.files('input_data/processed-rasqual-snps', full.names = TRUE)) {
  timepoint <- str_extract(score_path, 'time_\\d+')
  outpath <- paste('output_data/eQTL_pvalues_per_DHS-dataset_boxplot_', timepoint, '.png', sep = '')
  snp_score <- read_delim(score_path, delim = '\t', col_names = TRUE)
  p <- gen_boxplot(snp_score, bed_snp_map, timepoint)
  print(p)
  ggsave(outpath, plot = p)
}
