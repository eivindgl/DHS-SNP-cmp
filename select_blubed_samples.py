import json
import glob
import os
import collections
import shutil

short_name = {
        'naive thymus-derived CD4-positive, alpha-beta T cell': 'T-naive',
        'regulatory T-lymphocyte': 'T-reg',
        'T-helper 1 cell': 'Th-1',
        'T-helper 2 cell': 'Th-2',
        'T-helper 17 cell': 'Th-17'
}

skip_list = {
        'ENCFF001UUF_meta.json': 'Extremely low read depth (see https://www.encodeproject.org/experiments/ENCSR000EIG/)'
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
out_dir = 'output_data/bluebed_shortname'
if not os.path.isdir(out_dir):
    os.makedirs(out_dir)

seen_times = collections.defaultdict(int)

log_file = open(os.path.join(out_dir, 'sample_map.tsv'), 'w')
print('short_name', 'src', 'src_desc', 'dst', sep='\t', file=log_file)

for p in json_files:
    if os.path.basename(p) in skip_list:
        print('Skipping {} because: {}'.format(p, skip_list[os.path.basename(p)]))
        continue
    with open(p) as f:
        m = json.load(f)
    long_name = biosample_name(m)
    peak = peak_type(m)
    name = short_name.get(long_name, long_name)
    if peak == 'bed':# and not is_biological_replicate(m):
        bed_path = p[:(len(p) - len('_meta.json'))] + '.bed'
        fname = '{}_{}.bed'.format(name, seen_times[name])
        dst = os.path.join(out_dir, fname)
        seen_times[name] += 1
        print('Copying {} to {}'.format(name, dst))
        shutil.copyfile(bed_path, dst)
        print(name, bed_path, p, dst, sep='\t', file=log_file)

