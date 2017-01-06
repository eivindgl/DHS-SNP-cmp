pacman::p_load(
  tidyverse,
  rtracklayer,
  biomaRt
)

gwas_snps <- import('input_data/CeD_tag_SNPs.bed')
gwas_snps_table <- read_tsv('input_data/CeD_tag_SNPs.bed', col_names = c('chrom', 'snp_start', 'snp_end', 'SNP'))
eqtl_table <- read_tsv('output_data/rasqual_low_pvalue.tsv')

ensembl_version = 'feb2014.archive.ensembl.org'
mart <- useMart(
  'ENSEMBL_MART_ENSEMBL',
  host = ensembl_version,
  dataset = 'hsapiens_gene_ensembl')

eqtl_df <- getBM(
  attributes = c('chromosome_name', 'start_position', 'end_position', 'external_gene_id'),
  filters = 'ensembl_gene_id',
  values = eqtl_table$gene,
  mart = mart)  %>%
  dplyr::rename(
    chrom = chromosome_name,
    start = start_position,
    end = end_position,
    hgnc = external_gene_id
  ) %>%
  as_tibble

eqtl <- makeGRangesFromDataFrame(
  eqtl_df,
  keep.extra.columns = TRUE)
seqlevelsStyle(eqtl) <- 'UCSC'

gwas_snps_1mb <- gwas_snps %>% flank(5e5, both = TRUE)

hits <- findOverlaps(gwas_snps_1mb, eqtl)

eqtl_gwas <- bind_cols(
  eqtl_df[to(hits), ],
  gwas_snps_table[from(hits), ]
) %>%
  as_tibble %>%
  dplyr::select(
    chrom, start, end, hgnc, SNP, SNP_pos = snp_start)

eqtl_gwas %>%
  arrange(chrom, start) %>%
  write_tsv('output_data/eQTL_genes_within_1MB_of_CeD_tag_SNP.tsv')
