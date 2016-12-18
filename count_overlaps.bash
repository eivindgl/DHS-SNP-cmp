set -u
set -e
base_out=output_data/overlap
SNPs_orig=output_data/CEL_SNPs.bed
SNPs=$base_out/CEL_SNPs.bed
regions_orig=output_data/Immunobase_Celiac_Disease.bed
regions=$base_out/Immunobase_Celiac_Disease.bed

mkdir -p $base_out

set -x

[ -f $SNPs ] || awk '{print "chr" $0}' $SNPs_orig | sort-bed - > $SNPs
[ -f $regions ] || awk '{print "chr" $0}' $regions_orig | sort-bed - > $regions

dhs=input_data/bluebed_dhs/s6998/DNaseI/ENCFF001EJF.bed

#bedmap --skip-unmapped --echo --echo-map-id-uniq $dhs $SNPs
