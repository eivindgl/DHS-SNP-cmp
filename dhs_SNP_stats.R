pacman::p_load(
  GenomicRanges,
  rtracklayer,
  magrittr,
  purrr,
  ggplot2,
  stringr,
  tidyverse
  )

region <- import('output_data/overlap/Immunobase_Celiac_Disease.bed')
dhs <- import('input_data/bluebed_dhs/s6997/DNaseI/ENCFF001CKO.bed')
snp <- import('output_data/overlap/CEL_SNPs.bed')

genome_coverage <- function(x, by = 1) {
  cov_bp = x %>% coverage %>% sum %>% sum
  cov_bp / by
}

collect_stats <- function(snp, dhs, name = 'unknown') {
  rdhs <- GenomicRanges::intersect(dhs, region)
  list(
    name = name,
    tot_cov = genome_coverage(dhs, by = 1e3),
    reg_cov = genome_coverage(rdhs, by = 1e3),
    reg_num = length(rdhs),
    num_snp = genome_coverage(GenomicRanges::intersect(rdhs, snp)),
    tot_snp = length(snp)
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
summaries <- list()
for (name in names(beds)) {
  print(name)
  summaries[[name]] <- collect_stats(snp, beds[[name]], name)
}
df <- plyr::ldply(summaries, data.frame, .id = NULL) %>% tibble::as_tibble()

df$group <- str_sub(df$name, 1, 3)
df %>%
  #dplyr::filter(reg_cov < 400) %>%
  ggplot() + geom_point(aes(x = reg_cov, y = num_snp, color = group))
#ggplot(df) + geom_point(aes(x = reg_cov, y = num_snp))

df %>% mutate(snp_cov_kb = num_snp / reg_cov) %>% ggplot + geom_histogram(aes(snp_cov_kb))
