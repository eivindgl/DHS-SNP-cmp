# csvcut is part of csvkit (pip3 install csvkit)
# csvcut -n <csvfile> # lists all column names
# awk does 3 things:
# 	1. adds 'chr' prefix
#	2. adds +1 to end column 
#	3. removes everything after : at the end.
# 	(so rs12312:12321:A:T becomes just rs12312)
for x in *.csv ; do
  csvcut -c 9,10,10,8 $x | csvformat -T | tail -n+2 | awk '{$3=$3+1; sub(/:.*/,""); print "chr" $0}' OFS='\t' > ${x%.csv}.bed
  python find_tag_regions.py $x | bedtools sort | bedtools merge > ${x%.csv}.regions.bed
done
