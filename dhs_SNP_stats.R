pacman::p_load(
  GenomicRanges,
  rtracklayer,
  magrittr
  )

region <- import('output_data/overlap/Immunobase_Celiac_Disease.bed')
dhs <- import('input_data/bluebed_dhs/s6997/DNaseI/ENCFF001CKO.bed')
snp <- import('output_data/overlap/CEL_SNPs.bed')

genome_coverage <- function(x, by = 1) {
  cov_bp = x %>% coverage %>% sum %>% sum
  cov_bp / by
}

collect_stats <- function(snp, dhs) {
  rdhs <- intersect(dhs, region)
  c(
    tot_cov = genome_coverage(dhs, by = 1e3),
    reg_cov = genome_coverage(rdhs, by = 1e3),
    reg_num = length(rdhs),
    num_snp = genome_coverage(intersect(rdhs, snp)),
    tot_snp = length(snp)
  )
}

collect_stats(snp, dhs)
