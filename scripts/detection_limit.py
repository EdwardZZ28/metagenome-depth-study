#!/usr/bin/env python3
"""What abundance of ARG does a given depth fail to see?

For each sample, the 25 Gbp subtype profile is the reference. For each shallower
rung, two numbers: the highest abundance among subtypes that rung missed (an
effective detection floor), and the share of total abundance those misses carry.
"""
import glob, os, re, csv, collections, statistics as st

R = '/QRISdata/Q6636/depth_study_results/stageA'
PLAN = os.path.expanduser('~/depth_study/sample_plan.tsv')
TIERS = [2, 5, 10, 15, 20]
MATRICES = ['influent', 'effluent', 'preserved_wastewater', 'river_water', 'soil']
plan = {r['sample']: r['matrix'] for r in csv.DictReader(open(PLAN), delimiter='\t')}


def load(d):
    f = f'{d}/argoap/normalized_16S.subtype.txt'
    out = {}
    if not os.path.exists(f):
        return out
    for l in [x for x in open(f).read().splitlines() if x.strip()][1:]:
        p = l.split('\t')
        if len(p) >= 2:
            try: out[p[0]] = float(p[-1])
            except ValueError: pass
    return out


by = collections.defaultdict(dict)
for d in sorted(glob.glob(f'{R}/*')):
    m = re.match(r'(.+)__t(\d+)__s1$', os.path.basename(d))
    if m:
        by[m.group(1)][int(m.group(2))] = load(d)

floor = collections.defaultdict(lambda: collections.defaultdict(list))
share = collections.defaultdict(lambda: collections.defaultdict(list))
for s, d in by.items():
    if 25 not in d or len(d[25]) < 50:
        continue
    ref = d[25]
    mx = plan.get(s, '?')
    for t in TIERS:
        if t not in d or not d[t]:
            continue
        missed = [ref[k] for k in ref if k not in d[t]]
        found = [ref[k] for k in ref if k in d[t]]
        if missed and found:
            floor[mx][t].append(max(missed))
            share[mx][t].append(100 * sum(missed) / sum(ref.values()))


def show(title, tab, fmt):
    print(f'\n=== {title}')
    print(f"{'matrix':22s}" + ''.join(f'{t:>11d}' for t in TIERS))
    for mx in MATRICES:
        if mx not in tab:
            continue
        row = ''.join(f'{st.median(tab[mx][t]):{fmt}}' if tab[mx][t] else '          -'
                      for t in TIERS)
        print(f'{mx:22s}{row}')


show('effective detection floor: highest copies/16S among subtypes this depth missed',
     floor, '11.2e')
show('abundance missed, as % of total ARG abundance at 25 Gbp', share, '10.2f')
