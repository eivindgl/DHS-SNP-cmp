pacman::p_load(
  tidyverse,
  magrittr,
  stringr,
  assertthat,
  biomaRt
)

ensembl_version = 'feb2014.archive.ensembl.org'
mart <- useMart(
  'ENSEMBL_MART_ENSEMBL',
  host = ensembl_version,
  dataset = 'hsapiens_gene_ensembl')

exp_path_to_name <- function(p) {
  assert_that(is_scalar_character(p))
  name <- str_extract(p, 'time_\\d+')
  if (is.na(name)) {
    name <- str_match(p, '([^_]+).merged')[[2]]
  }
  name
}

load_exp <- function(p, p_lim = 1e-3) {
  df <- read_tsv(p, col_types = 'ccdddd', progress = FALSE)
  df$FDR <- p.adjust(df$pvalue, method = 'BH')
  name <- exp_path_to_name(p)
  assert_that(is.character(name), nchar(name) > 0)
  cat(name, ", min p-value: ", min(df$pvalue), ', min fdr-adjusted pvalue: ', min(df$FDR), '\n')
  df %<>% filter(df$pvalue < p_lim)
  cat('Limiting to ', nrow(df), ' values with a p-value < ', p_lim, '.\n')
  df %>%
    mutate(fdr_cat = cut(FDR, c(1e-15, 1e-3, 1e-2, 1e-1, 1))) %>%
    group_by(fdr_cat) %>%
    summarise(n_fdr = n()) %>%
    print()

  df$experiment <- name
  return(df)
}
#
# Operate on lead SNPs
#

load_rasqual_exp <- function(paths) {
  n <- paths %>% map_chr(exp_path_to_name)
  paths %>%
    map(load_exp) %>%
    setNames(nm = n)
}

combine_common_experiments <- function(exp_dfs, exp_id = 'time') {
  # by default, combines genes with low p-values from individual experiments
  exp_dfs %>%
    bind_rows %>%
    filter(str_detect(experiment, exp_id))
}

download_gene_metadata <- function(ensembl_ids) {
  getBM(
    attributes = c('ensembl_gene_id', 'external_gene_id', 'gene_biotype'),
    filters = 'ensembl_gene_id',
    values = ensembl_ids,
    mart = mart) %>%
    rename(
      gene = ensembl_gene_id,
      hgnc = external_gene_id,
      type = gene_biotype) %>%
    as_tibble
}

add_gene_meta <- function(df) {
  download_gene_metadata(df$gene) %>%
    inner_join(df, by = 'gene') %>%
    dplyr::select(hgnc, experiment, FDR, pvalue, type, everything())
}
#
# Operate on lead SNPs
#

score_paths <- list.files('input_data/raw-rasqual-snps/gsTCC_lead_eigenfmt', full.names = TRUE)
exp_dfs <- load_rasqual_exp(score_paths)
time_df <- exp_dfs %>%
  combine_common_experiments(exp_id = 'time') %>%
  add_gene_meta()
time_df %>%
  arrange(pvalue) %>%
  write_tsv('output_data/rasqual_low_pvalue.tsv')
all_df <- exp_dfs[['all']] %>%
  dplyr::filter(FDR < 0.05) %>%
  add_gene_meta() %>%
  arrange(FDR)
all_df %>%
  write_tsv('output_data/rasqual_all_sig_FDR.tsv')

# very poor overlap between timepoint df and all df...
# check very carefully that this does not indicate a bug somewhere...
# 1. what about mapping between vcf and counts?
all_df %>%
  semi_join(time_df, by = 'gene') %>%
  arrange(FDR) %>%
  print(n = Inf)
#
# Operate on all SNPs
#

full_score_paths <- list.files('input_data/raw-rasqual-snps/gsTCC_all_eigenfmt', full.names = TRUE)
full_exp_dfs <- load_rasqual_exp(full_score_paths)
full_time_df <- full_exp_dfs %>%
  combine_common_experiments(exp_id = 'time') %>%
  add_gene_meta()
full_time_df %>%
  arrange(pvalue) %>%
  write_tsv('output_data/rasqual_full_low_pvalue.tsv')
full_all_df <- full_exp_dfs[['all']] %>%
  dplyr::filter(FDR < 1e-2) %>%
  add_gene_meta() %>%
  arrange(FDR)
full_all_df %>%
  write_tsv('output_data/rasqual_full_all_sig_FDR.tsv')

full_all_df %>%
  semi_join(full_time_df, by = 'gene') %>%
  arrange(FDR) %>%
  print(n = Inf)


df <- full_score_paths %>%
  map(load_exp) %>%
  bind_rows

df %>% write_tsv('output_data/rasqual_full_low_pvalue.tsv')
