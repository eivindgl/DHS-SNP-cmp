pacman::p_load(
  tidyverse,
  magrittr,
  stringr,
  biomaRt
)

ensembl_version = 'feb2014.archive.ensembl.org'
mart <- useMart(
  'ENSEMBL_MART_ENSEMBL',
  host = ensembl_version,
  dataset = 'hsapiens_gene_ensembl')


df <- read_tsv('input_data/raw-rasqual-snps/gsTCC_lead_eigenfmt/lead_SNP_gsTCC_time_180.merged.eigenmt.tsv')

load_exp <- function(p, p_lim = 1e-3) {
  df <- read_tsv(p, col_types = 'ccdddd', progress = FALSE)
  df$FDR <- p.adjust(df$pvalue, method = 'BH')
  name <- str_extract(p, 'time_\\d+')
  cat(name, ", min p-value: ", min(df$pvalue), ', min fdr-adjusted pvalue: ', min(df$FDR), '\n')
  df %<>% filter(df$pvalue < p_lim)
  cat('Limiting to ', nrow(df), ' values with a p-value < ', p_lim, '.\n')
  df$experiment <- name
  return (df)
}

score_paths <- list.files('input_data/raw-rasqual-snps/gsTCC_lead_eigenfmt', full.names = TRUE)
df <- score_paths %>%
  map(load_exp) %>%
  bind_rows

hgnc <- getBM(
  attributes = c('ensembl_gene_id', 'external_gene_id', 'gene_biotype'),
  filters = 'ensembl_gene_id',
  values = df$gene,
  mart = mart) %>%
  rename(
    gene = ensembl_gene_id,
    hgnc = external_gene_id,
    type = gene_biotype) %>%
  as_tibble

df <- inner_join(df, hgnc) %>%
  dplyr::select(hgnc, experiment, pvalue, type, everything()) %>%
  arrange(pvalue)

df %>% write_tsv('output_data/rasqual_low_pvalue.tsv')
