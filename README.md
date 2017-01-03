# DHS-SNP-cmp
Extract common SNPs in GWAS disease associated regions.
Find overlap with DHS-sites from various T-cell subsets.
Generate a summary.

## Input Data

### Immunobase summary for Celiac Disease
Easily available from https://www.immunobase.org/disease/CEL/

Downloaded file has UTF-8 BOM, which trips Python.
Can be removed with
```
tail --bytes=+4 input_data/Immunobase_Celiac_Disease.csv > tmp
mv tmp input_data/Immunobase_Celiac_Disease.csv
```

### LD SNPs
CeD tag SNPs are located in `input_data/CeD_tag_SNPs.bed`.
I used the tag list with the web service http://raggr.usc.edu to find
SNPs in strong LD.

The paramters were:
* MAF 0.001
* R^2 > 0.9
* CEU european population (1000g)

