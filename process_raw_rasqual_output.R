pacman::p_load(
  tidyverse,
  stringr
)

if (!dir.exists('input_data/processed-rasqual-snps')) {
  dir.create('input_data/processed-rasqual-snps')
}

ldf <- read_delim('input_data/LD_SNPs/CeD_LD_SNPs_Iris_maf-.001_r2-9.bed', delim = '\t',
                  col_names = c('chrom', 'start', 'end', 'snps'))
# a few LD snps does not start with the rs prefix. I skip these for now.
# TODO figure out why these values are present!
ldf <- ldf %>% filter(str_detect(snps, '^rs'))

raw_paths <- list.files('input_data/raw-rasqual-snps/gsTCC_all_eigenfmt',
                        full.names = TRUE,
                        pattern = 'tsv$')

for (p in raw_paths) {
  timepoint <- str_extract(p, 'time_\\d+')
  outpath <- paste('input_data/processed-rasqual-snps/', timepoint, '.tsv', sep = '')
  rdf <- read_delim(p, delim = '\t', col_names = TRUE)
  # Each SNP may have multiple pvalues (for multiple genes). Use the min pvalue per SNP
  # Not sure how to handle SNPs distant to all genes
  # (so rasqual had too little information for any calculations)
  # I could replace them with e.g. 1. But it is perhaps better to skip them.
  #df %>% mutate(pvalue = if_else(is.na(pvalue), 1, pvalue))
  df <- right_join(rdf, ldf, by = 'snps') %>% rename(SNP = snps) %>%
    group_by(chrom, start, end, SNP) %>%
    summarise(pvalue = min(pvalue))
  df %>% write_delim(path = outpath,
                     delim = '\t', col_names = TRUE)
}


