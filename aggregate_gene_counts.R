pacman::p_load(
  tidyverse,
  stringr
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
