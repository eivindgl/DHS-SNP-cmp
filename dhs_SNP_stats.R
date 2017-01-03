pacman::p_load(
  GenomicRanges,
  rtracklayer,
  magrittr,
  purrr,
  ggplot2,
  stringr,
  tidyverse
  )

region <- rtracklayer::import('input_data/LD_SNPs/CeD_LD_SNPs_Iris_maf-.001_r2-9.regions.bed')
snp <- rtracklayer::import('input_data/LD_SNPs/CeD_LD_SNPs_Iris_maf-.001_r2-9.bed')

genome_coverage <- function(x, by = 1) {
  cov_bp = x %>% coverage %>% sum %>% sum
  cov_bp / by
}

collect_stats <- function(region, snp, dhs, name = 'unknown') {
  rdhs <- GenomicRanges::intersect(dhs, region)
  list(
    name = name,
    tot_cov = genome_coverage(dhs, by = 1e3),
    reg_cov = genome_coverage(rdhs, by = 1e3),
    reg_num = length(rdhs),
    num_snp = length(GenomicRanges::intersect(rdhs, snp)),
    tot_snp = length(GenomicRanges::intersect(region, snp))
  )
}

read_bedfiles <- function(bedfiles) {
  xs <- list()
  for (x in bedfiles) {
    name <- tools::file_path_sans_ext(basename(x))
    xs[[name]] <- import(x)
  }
  xs
}

bluebed_paths <- list.files('output_data/bluebed_shortname', full.names = TRUE, pattern = 'bed$')
gsTCC_paths <- list.files('input_data/gsTCC_dhs', full.names = TRUE, pattern = 'bed$')
beds <- read_bedfiles(c(bluebed_paths, gsTCC_paths))

#
# compute summaries
#
get_summaries <- function(xs, region, snp) {
  summaries <- list()
  for (name in names(xs)) {
    print(name)
    summaries[[name]] <- collect_stats(region, snp, xs[[name]], name)
  }
  summaries
}
summaries <- get_summaries(beds, region, snp)

df <- plyr::ldply(summaries, data.frame, .id = NULL) %>%
  as_tibble %>%
  mutate(snp_per_kb = num_snp / reg_cov)

tmp <- df %>%
  filter(!str_detect(df$name, '^(TCC|merge)')) %>%
  mutate(group = str_replace(name, '_\\d$', ''))

df %<>%
  filter(str_detect(df$name, '^(TCC|merge)')) %>%
  mutate(group = paste('TCC', str_extract(name, '\\d+$'), sep = '-')) %>%
  bind_rows(tmp)

df %<>%
  mutate(group = fct_relevel(group,
                             'TCC-0', 'TCC-10', 'TCC-30', 'TCC-180',
                             'T-naive', 'Th-1', 'Th-2', 'Th-17', 'T-reg', 'Jurkat'))


df %>%
  filter(!str_detect(name, 'merged')) %>%
  ggplot(aes(x = group, y = snp_per_kb, color = group)) +
  geom_boxplot() +
  xlab('T-cell type') +
  ylab('SNPs per kbp') +
  ggtitle('Average number of SNPs per kbp open chromatin in strong LD with CeD tag SNPs')

ggsave('output_data/DHS_CeD_LD-SNPs_boxplot.png')
