
gced_df %>%
  mutate(source = ifelse(str_detect(name, 'TCC'),
                         'gsTCC', 'Encode')) %>%
  ggplot(aes(x = group, y = snp_per_kb, color = source)) +
  geom_boxplot() +
  xlab('T-cell type') +
  ylab('SNPs per kbp open chromtin at CeD loci') +
  ggtitle('Enrichment of Celiac Disease GWAS SNPs in DNase-sites') #+
  theme(legend.position = "none") # removes legend
ggsave('~/Desktop/simpe_CeD_DNase_overlap_comparison.png')
