# DHS-SNP-cmp
## Extract common SNPs in GWAS disease associated regions.
See extract\_SNPs.py script. Generates SNP frequency list
for each region. This is again merged to a common bed file of SNPs.
## Find overlap with DHS-sites from various T-cell subsets.
see count\_overlaps.bash
## Generate a summary.

## Input Data

### Immunobase summary for Celiac Disease
Easily available from https://www.immunobase.org/disease/CEL/

Downloaded file has UTF-8 BOM, which trips Python.
Can be removed with
```
tail --bytes=+4 input_data/Immunobase_Celiac_Disease.csv > tmp
mv tmp input_data/Immunobase_Celiac_Disease.csv
```

### 1000g Phase 3
I use this a the source for SNPs. I only extract common SNPs (maf 0.01)
in the immunobase regions
