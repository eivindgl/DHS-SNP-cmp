import json
import glob
import os

short_name = {
        'naive thymus-derived CD4-positive, alpha-beta T cell': 'T-naive',
        'regulatory T-lymphocyte': 'T-reg',
        'T-helper 1 cell': 'Th-1',
        'T-helper 2 cell': 'Th-2'
}

def coverage(p):
    def length(line):
        start, end = line.split()[1:3]
        return int(end) - int(start)

    with open(p) as f:
        return '{0:.2f}MB'.format(sum(map(length, f)) / 1e6)


def biosample_name(m):
    return m['sample_info']['biosample_term_name']

def peak_type(m):
    return m['extra_metadata']['file_type'].split()[0]

def is_biological_replicate(m):
    return 'biological_replicate_number' in m['extra_metadata']

json_files = glob.glob('input_data/bluebed_dhs/*/DNaseI/*.json')

for p in json_files:
    with open(p) as f:
        m = json.load(f)
    long_name = biosample_name(m)
    peak = peak_type(m)
    name = short_name.get(long_name, long_name)
    if peak == 'bed' and not is_biological_replicate(m):
        bed_path = p[:(len(p) - len('_meta.json'))] + '.bed'
        print(bed_path, os.path.isfile(bed_path))
        print(name, peak, p, coverage(bed_path))
