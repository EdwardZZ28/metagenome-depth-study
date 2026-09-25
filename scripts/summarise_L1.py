#!/usr/bin/env python3
"""Collect per-dataset metrics from the stage A pipeline output into one table.

Usage:  python3 scripts/summarise_L1.py [RESULTS_DIR] [OUT_TSV]

RESULTS_DIR must contain one directory per dataset, named
<sample>__t<tier>__s<seed>, holding:
    argoap/normalized_16S.subtype.txt
    size.txt                      (contigs=, assembly_bp=, bins=)
    genomad/<id>_summary/<id>_plasmid_summary.tsv
    genomad/<id>_summary/<id>_virus_summary.tsv
    rgi/<id>.txt

See METHODS.md section 5 for the definition of every column.
"""
import glob
import os
import re
import sys

RESULTS = sys.argv[1] if len(sys.argv) > 1 else '/QRISdata/Q6636/depth_study_results/stageA'
OUT = sys.argv[2] if len(sys.argv) > 2 else os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'data', 'L1_summary.tsv')

HIGH_RISK = re.compile(
    r'\b(NDM|KPC|OXA-48|OXA-23|IMP-|VIM-|MCR-|mcr-|VanA|vanA|VanB|vanB|CTX-M|TEM-|SHV-)',
    re.I)
ID_MIN = 80.0
COLS = ['sample', 'tier', 'arg_subtypes', 'copies_16S', 'contigs', 'assembly_bp',
        'bins', 'plasmid', 'conjugative', 'virus', 'rgi_ps', 'rgi_id80',
        'highrisk', 'highrisk_genes']


def to_int(x):
    try:
        return int(x)
    except (TypeError, ValueError):
        return 0


def read_argoap(path):
    if not os.path.exists(path):
        return None, None
    lines = [l for l in open(path).read().splitlines() if l.strip()]
    total = 0.0
    for l in lines[1:]:
        try:
            total += float(l.split('\t')[-1])
        except ValueError:
            pass
    return len(lines) - 1, round(total, 4)


def read_genomad(d):
    summaries = glob.glob(f'{d}/genomad/*_summary')
    if not summaries:
        return 0, 0, 0
    base = os.path.basename(summaries[0]).replace('_summary', '')
    plasmid = conjugative = virus = 0
    fp = f'{summaries[0]}/{base}_plasmid_summary.tsv'
    if os.path.exists(fp):
        lines = open(fp).read().splitlines()
        plasmid = max(0, len(lines) - 1)
        if lines:
            header = lines[0].split('\t')
            if 'conjugation_genes' in header:
                ci = header.index('conjugation_genes')
                for l in lines[1:]:
                    f = l.split('\t')
                    if len(f) > ci and f[ci].strip() not in ('', 'NA'):
                        conjugative += 1
    fv = f'{summaries[0]}/{base}_virus_summary.tsv'
    if os.path.exists(fv):
        virus = max(0, len(open(fv).read().splitlines()) - 1)
    return plasmid, conjugative, virus


def read_rgi(d):
    files = glob.glob(f'{d}/rgi/*.txt')
    if not files:
        return 0, 0, set()
    lines = open(files[0]).read().splitlines()
    if not lines:
        return 0, 0, set()
    header = lines[0].split('\t')
    need = ('Cut_Off', 'Best_Identities', 'Best_Hit_ARO')
    if not all(c in header for c in need):
        return 0, 0, set()
    c_cut, c_id, c_aro = (header.index(c) for c in need)
    n_ps = n_id80 = 0
    genes = set()
    for l in lines[1:]:
        f = l.split('\t')
        if len(f) <= max(c_cut, c_id, c_aro) or f[c_cut] not in ('Perfect', 'Strict'):
            continue
        n_ps += 1
        try:
            identity = float(f[c_id])
        except ValueError:
            continue
        if identity >= ID_MIN:
            n_id80 += 1
            if HIGH_RISK.search(f[c_aro]):
                genes.add(f[c_aro])
    return n_ps, n_id80, genes


rows = []
for d in sorted(glob.glob(f'{RESULTS}/*')):
    dataset_id = os.path.basename(d)
    m = re.match(r'(.+)__t(\d+)__s(\d+)$', dataset_id)
    if not m:
        continue
    r = {'sample': m.group(1), 'tier': int(m.group(2))}

    subtypes, copies = read_argoap(f'{d}/argoap/normalized_16S.subtype.txt')
    r['arg_subtypes'] = subtypes if subtypes is not None else ''
    r['copies_16S'] = copies if copies is not None else ''

    size = {}
    if os.path.exists(f'{d}/size.txt'):
        for token in open(f'{d}/size.txt').read().split():
            if '=' in token:
                k, v = token.split('=', 1)
                size[k] = v
    r['contigs'] = to_int(size.get('contigs'))
    r['assembly_bp'] = to_int(size.get('assembly_bp'))
    r['bins'] = to_int(size.get('bins'))

    r['plasmid'], r['conjugative'], r['virus'] = read_genomad(d)
    r['rgi_ps'], r['rgi_id80'], genes = read_rgi(d)
    r['highrisk'] = len(genes)
    r['highrisk_genes'] = ';'.join(sorted(genes))
    rows.append(r)

os.makedirs(os.path.dirname(OUT), exist_ok=True)
with open(OUT, 'w') as fh:
    fh.write('\t'.join(COLS) + '\n')
    for r in sorted(rows, key=lambda x: (x['sample'], x['tier'])):
        fh.write('\t'.join(str(r.get(c, '')) for c in COLS) + '\n')
print(f'{len(rows)} datasets written to {OUT}')
