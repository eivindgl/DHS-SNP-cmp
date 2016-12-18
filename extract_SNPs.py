import os
import csv
import re
import glob

class Region:

    def __init__(self, chrom, start, end, name):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.name = name

    def __str__(self):
        return '{chrom}\t{start}\t{end}\t{name}'.format(**self.__dict__)

    def extract_SNPs_command(self, vcf_file, out_dir='.', maf=0.01):
        return '''\
vcftools --gzvcf {vcf_file} --out {out_prefix} --remove-indels \
--maf {maf} --freq2 --chr {chrom} --from-bp {start} --to-bp {end}'''.format(
        vcf_file=vcf_file, maf=maf,
        out_prefix=os.path.join(out_dir, 'region_' + self.name),
        **self.__dict__)


#
# Read Input
#
def read_immuobase_regions(ib_path):
    p = re.compile(r'chr(\d+):(\d+)-(\d+)')
    regions = []
    with open(ib_path) as f:
        for row in csv.DictReader(f):
            m = p.search(row['Position'])
            if m is not None:
                chrom, start, end = m.groups()
                region_name = row['Region'].split()[1] # contains a 'CEL ' prefix
                region = Region(chrom,start,end, region_name)
                regions.append(region)
    return regions

#
# Write output regions bed 
#
def write_ib_bed(ib_bed_path, regions, overwrite=True):
    if os.path.isfile(ib_bed_path) and not overwrite:
        print('File exists. skipping', ib_bed_path)
        return
    with open(ib_bed_path, 'w') as f:
        for x in regions:
            print(x, file=f)
    print('Wrote to', immunobase_bed)

#
# Write extarctSNP bash script (no execute)
#
def write_extractSNP_script(script_path, vcf_dir, vcf_name_pattern,
        overwrite=True):
    if os.path.isfile(script_path) and not overwrite:
        print('File exists. skipping', script_path)
        return
    with open(script_path, 'w') as f:
        print('set -u', file=f)
        print('set -e', file=f)
        for x in regions:
            vcf_path = os.path.join(vcf_dir,
                    vcf_name_pattern.format(chrom=x.chrom))
            print(x.extract_SNPs_command(vcf_path), file=f)
    print('Wrote to', gen_script_path)

#
# Merge freq counts to bed file with region names
#
def extract_region_name(path):
    '''Extract immunobase name from vcf count file path.

    >>> extract_region_name('dst/region_10p15.1.frq')
    '10p15.1'
    '''
    fname = os.path.basename(path)
    name = os.path.splitext(fname)[0]
    return name.split('_', 1)[1]

def print_bed_stream(p, name, fh):
    with open(p) as f:
        for row in csv.DictReader(f, dialect='excel-tab'):
            line = '{CHROM}\t{POS}\t{end}\t{name}'.format(
                    end=int(row['POS']) + 1,
                    name=name,
                    **row)
            print(line, file=fh)

def merge_SNP_counts(count_dir, dst_path, overwrite=True):
    '''Merges SNPs from all region count files into a single bed file.
    '''
    if os.path.isfile(dst_path) and not overwrite:
        print('File exists. Skipping', dst_path)
        return

    with open(dst_path, 'w') as f:
        for p in glob.glob('{}/*.frq'.format(count_dir)):
            name = extract_region_name(p)
            print_bed_stream(p, name, f)
    print('Wrote', dst_path)

#
# Exec Flow
#
def exec():
    out_dir = 'output_data'
    immunobase_regions = 'input_data/Immunobase_Celiac_Disease.csv'
    immunobase_bed = os.path.join(out_dir, 'Immunobase_Celiac_Disease.bed')
    gen_script_path= os.path.join(out_dir, 'vcf_region_count', 'extract_SNPs.bash')
    count_dir = os.path.dirname(gen_script_path)
    vcf_dir = '/apps/data/1000G/phase3/20130502'
    vcf_name_pattern = 'ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz'

    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)
    if not os.path.isdir(count_dir):
        os.makedirs(count_dir)
    regions = read_immuobase_regions(immunobase_regions)
    write_ib_bed(immunobase_bed, regions, overwrite=False)
    write_extractSNP_script(gen_script_path, vcf_dir, vcf_name_pattern,
        overwrite=False)
    merge_SNP_counts(count_dir, os.path.join(out_dir, 'CEL_SNPs.bed'), overwrite=False)

if __name__ == '__main__':
    exec()
