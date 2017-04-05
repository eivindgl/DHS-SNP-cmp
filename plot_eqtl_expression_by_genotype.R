pacman::p_load(
  tidyverse,
  forcats,
  stringr,
  VariantAnnotation
)
vst_path = 'out/post_proc/normalized_gene_counts.tsv'
vst <- read_tsv(vst_path) %>%
  gather(experiment, vst, -gene_id) %>%
  mutate(timepoint = str_extract(experiment, '^[^_]+'),
         timepoint = factor(timepoint, levels = c('time0', 'time10', 'time30', 'time180')),
         sample = str_extract(experiment, '[^_]+$'))
eqtl <- read_tsv('out/post_proc/CeD_LD_subset/LD_subset_all_all.tsv') %>%
  dplyr::rename(snps = rs_id, hgnc = external_gene_id)
vcf <- readVcf('out/preprocessing/all.vcf')
snps <- unique(eqtl$snps)
gt <- geno(vcf)$GT %>%
  as.data.frame %>%
  rownames_to_column('SNP') %>%
  as_tibble %>%
  filter(SNP %in% snps) %>%
  gather(sample, genotype, -SNP) %>%
  mutate(genotype = fct_recode(genotype,
                               ref = '0|0',
                               heterozygous = '0|1',
                               heterozygous = '1|0',
                               alt = '1|1'))

gene_snp_map <- eqtl %>%
  dplyr::select(hgnc, gene_id, SNP = snps)

df <- gene_snp_map %>%
  inner_join(vst) %>%
  inner_join(gt)

# print highly expressed genes
df %>%
  group_by(hgnc, gene_id) %>%
  summarise(n = n(), median = median(vst), mean = mean(vst)) %>%
  arrange(desc(median))

dir.create('output_data/eqtl_expression', showWarnings = FALSE)
for (gene in unique(df$hgnc)) {
  p <- df %>%
    filter(hgnc == gene) %>%
    ggplot(aes(genotype, vst, color = genotype)) +
    geom_boxplot() +
    geom_jitter(width = 0.2) +
    facet_grid(hgnc ~ timepoint) +
    ggtitle(gene)
  outpath <- paste('output_data/eqtl_expression/', gene, '.png', sep = '')
  ggsave(outpath, plot = p)
}

df %>%
  filter(hgnc == 'CD28') %>%
  mutate(source = ifelse(str_detect(sample, 'TCC-'),
                         'Koning', 'Sollid')) %>%
  ggplot() +
  geom_boxplot(aes(genotype, vst, fill = genotype)) +
  geom_jitter(width = 0.2, aes(genotype, vst, color = source)) +
  facet_wrap(~ timepoint, nrow = 2) +
  ggtitle('CD28')


eqtl %>%
  group_by(hgnc) %>%
  summarise(pval = min(pval), FDR = min(FDR)) %>%
  filter(hgnc == 'CD28')
