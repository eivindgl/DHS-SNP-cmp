pacman::p_load(
  tidyverse,
  forcats,
  stringr
)

bed_snp_map <- read_csv('output_data/snp_dhs_map.csv')
snp_score <- read_delim('input_data/processed-rasqual-snps/time_180.tsv',
                        delim = '\t', col_names = TRUE)
df <- inner_join(bed_snp_map, snp_score)
df <- df %>% filter(!is.na(pvalue))

tmp <- df %>%
  filter(!str_detect(DHS, '^TCC')) %>%
  mutate(group = str_replace(DHS, '_\\d$', ''))

df <- df %>%
  filter(str_detect(DHS, '^TCC')) %>%
  mutate(group = paste('TCC', str_extract(DHS, '\\d+$'), sep = '-')) %>%
  bind_rows(tmp)

cdf <- df %>%
  group_by(DHS) %>%
  summarize(median = median(pvalue), mean = mean(pvalue)) %>%
  arrange(mean)

df <- df %>% mutate(DHS = factor(DHS, levels = cdf$DHS))

p <- df %>%
  ggplot(aes(DHS, pvalue, color = group)) +
  geom_boxplot() +
  scale_y_continuous(limits = c(0, 0.25)) +
  coord_flip()
p +
  ggtitle('Boxplot of eQTL p-values for SNPs in DHS sites')
ggsave('output_data/eQTL_pvalues_per_DHS-dataset_boxplot.png')
