import csv
import sys
import collections

ld_snps = collections.defaultdict(list)
entries = csv.DictReader(open(sys.argv[1]))
chrom_map = {}

for entry in entries:
    lead_snp = entry['SNP1 Name'].split(':', 1)[0]
    lead_chrom = entry['SNP1 Chr']
    ld_snp_pos = int(entry['SNP2 Pos'])
    ld_snps[lead_snp].append(ld_snp_pos)
    if not lead_chrom.startswith('chr'):
        lead_chrom = 'chr' + lead_chrom
    chrom_map[lead_snp] = lead_chrom

for lead_snp, snps in ld_snps.items():
    print(
            chrom_map[lead_snp], 
            min(snps),
            max(snps) + 1,
            lead_snp, 
            sep='\t')



