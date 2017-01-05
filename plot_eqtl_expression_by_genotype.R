pacman::p_load(
  tidyverse,
  forcats,
  stringr,
  VariantAnnotation
)
vst_dir = 'output_data/expression'
vst <- list.files(vst_dir, full.names = TRUE) %>%
  map(function(p) {
    tp_lvls <- c('t0', 't10', 't30', 't180')
    vst_name <- str_extract(p, 't\\d+')
    read_tsv(p) %>%
      gather(sample, vst, -geneid) %>%
      mutate(timepoint = factor(vst_name, levels = tp_lvls))
  }) %>%
  bind_rows()
eqtl <- read_tsv('output_data/rasqual_low_pvalue.tsv')
vcf <- readVcf('input_data/all.vcf')
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

gene_snp_map <- eqtl %>% dplyr::select(hgnc, geneid = gene, SNP = snps)

df <- gene_snp_map %>%
  inner_join(vst) %>%
  inner_join(gt)

# print highly expressed genes
df %>%
  group_by(hgnc, geneid) %>%
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

