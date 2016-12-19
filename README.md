# DHS-SNP-cmp
See if gsTCC cells are overlapping more SNPs in GWAS associated disease regions
than other T-cell types.

## Status
A course analysis works and finds a very linear relationships between coverage in CeD regions
and SNP overlap. This might be correct, but it is possbile that I do something stupid when
selecting SNPs. I currently select all snps with MAF > 1% in disease associated regions.
This gives me a lot of SNPs that seems to be quite uniformly distributed across the regions.
Therefore, it seems likely that the linear relationship between DHS coverage and #SNPs is
given by the uniform distribution of SNPs.

# Notes / Thoughts / Possible extentions
Maybe it makes more sense to only look at very common MAFs (%5+).
Rationale is that GWAS studies can usually only pick these up.

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
