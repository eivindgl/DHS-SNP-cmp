pacman::p_load(
  tidyverse,
  stringr,
  purrr
)
counts <- read_tsv('input_data/counts.tsv.gz')
names(counts)[1] <- 'geneid'
vst <- read_tsv('input_data/counts_VST.tsv.gz')
names(vst)[1] <- 'geneid'


# Compute median expression per gene per timepoint.
# Constructs a dataframe with one gene per row and median
# expression per timepoint in columns.
median_per_gene_timepoint <- function(df) {
  df %>%
    gather(sample, count, -geneid) %>%
    mutate(timepoint = str_extract(sample, 't\\d+')) %>%
    group_by(timepoint, geneid) %>%
    summarize(median_expression = median(count)) %>%
    ungroup %>%
    spread(timepoint, median_expression)
}
cdf <- median_per_gene_timepoint(counts)
vdf <- median_per_gene_timepoint(vst)

cdf %>% write_tsv('input_data/counts_median.tsv')
vdf %>% write_tsv('input_data/counts_VST_median.tsv')

#
# rename count table so it corresponds to vcf names
#
x <- vst %>%
  gather(sample, count, -geneid) %>%
  mutate(timepoint = str_extract(sample, 't\\d+'))
dfs <- split(x, x$timepoint)
for (timepoint in names(dfs)) {
  name_map_path <- paste('input_data/ASE_counts/time_',
                         str_sub(timepoint, start = 2),
                         '/sample_bam_map.tsv',
                         sep = '')
  if (!file.exists(name_map_path)) {
    cat('Missing file, skipping: ', name_map_path, '\n')
    next
  }
  extract_name <- compose(basename, tools::file_path_sans_ext)
  name_map <- read_tsv(name_map_path,
                       col_names = c('name', 'path')) %>%
    mutate(long_name = extract_name(path))

  dfs[[timepoint]] %>%
    inner_join(name_map, by = c('sample' = 'long_name')) %>%
    select(name, geneid, count) %>%
    spread(name, count) %>%
    write_tsv(paste('output_data/expression_vst_', timepoint, '_renamed.tsv', sep=''))

}

