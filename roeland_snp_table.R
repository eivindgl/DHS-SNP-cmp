# roeland asked about his (r^2 > 0.8) SNPs overlapping encode data sets.
pacman::p_load(
  tidyverse,
  stringr,
  purrr,
  GenomicRanges,
  rtracklayer
)
rdf <- read_csv('input_data/roeland_snps.csv') %>%
  dplyr::select(1:4)
roe <- rdf %>%
  GRanges()
seqlevelsStyle(roe) <- 'UCSC'

read_bedfiles <- function(bedfiles) {
  xs <- list()
  for (x in bedfiles) {
    name <- tools::file_path_sans_ext(basename(x))
    xs[[name]] <- import(x)
  }
  xs
}

bluebed_paths <- list.files('output_data/bluebed_shortname', full.names = TRUE, pattern = 'bed$')
gsTCC_paths <- list.files('input_data/gsTCC_dhs', full.names = TRUE, pattern = 'bed$') %>%
  purrr::discard(~ str_detect(.x, 'merged'))
beds <- read_bedfiles(c(bluebed_paths, gsTCC_paths))

has_overlap <- purrr::compose(`!`, is.na)
snp_overlap <- beds %>%
  map(function(x) {
    findOverlaps(roe, x, select = 'first') %>%
      has_overlap()
  }) %>%
  bind_cols

df <- bind_cols(rdf, snp_overlap)
df %>%
  write_tsv('output_data/roeland_snp_overlap.tsv')

x <- df %>%
  gather(experiment, DHS, -chr, -start, -end, -SNP) %>%
  mutate(group = ifelse(str_detect(experiment, '^TCC'),
                        paste('gsTCC_t', str_extract(experiment, '\\d+$'), sep = ''),
                        str_extract(experiment, '[^_]+')))

hits_cnt <- x %>%
  group_by(chr, start,end, SNP, group) %>%
  summarise(n_hits = sum(DHS))

tot_rep_cnt <- x %>%
  dplyr::distinct(experiment, group) %>%
  group_by(group) %>%
  summarize(total_rep = n())

hits_cnt %>%
  inner_join(tot_rep_cnt) %>%
  mutate(fraction = n_hits / total_rep) %>%
  dplyr::select(-n_hits, -total_rep) %>%
  spread(group, fraction) %>%
  write_tsv('output_data/roeland_snp_overlap_merged.tsv')
