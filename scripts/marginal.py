#!/usr/bin/env python3
"""Marginal return per Gbp within the 0-25 Gbp window, computed per sample."""
import csv, collections, os, statistics as st

R = [r for r in csv.DictReader(open(os.path.expanduser('~/depth_study/core_summary.tsv')),
                               delimiter='\t') if int(r['seed']) == 1]
INT = ['arg_subtypes','contigs','bins','plasmid','conjugative','virus','rgi_id80','bracken_species']
for r in R:
    r['tier'] = int(r['tier'])
    for k in INT:
        r[k] = int(r[k]) if r[k] else 0
    r['copies_16S'] = float(r['copies_16S']) if r['copies_16S'] else 0.0

by = collections.defaultdict(dict)
for r in R:
    by[(r['matrix'], r['sample'])][r['tier']] = r

STEPS = [(2,5),(5,10),(10,15),(15,20),(20,25)]
METRICS = ['copies_16S','arg_subtypes','rgi_id80','bracken_species',
           'contigs','plasmid','conjugative','virus','bins']

print('=== marginal gain per extra Gbp, median across samples (units per Gbp) ===')
print(f"{'matrix':22s}{'metric':18s}" + ''.join(f'{a}->{b:<6d}' for a,b in STEPS))
for mx in ['influent','effluent','preserved_wastewater','river_water','soil']:
    ss = {k:v for k,v in by.items() if k[0]==mx}
    if not ss: continue
    for m in METRICS:
        row = ''
        for a,b in STEPS:
            v = [(d[b][m]-d[a][m])/(b-a) for k,d in ss.items() if a in d and b in d]
            row += f'{st.median(v):9.1f}' if v else '        -'
        print(f'{mx:22s}{m:18s}{row}')
    print()

print('=== deceleration: last step (20->25) as a fraction of the first (2->5) ===')
print(f"{'matrix':22s}" + ''.join(f'{m[:11]:>13s}' for m in METRICS))
for mx in ['influent','effluent','preserved_wastewater','river_water','soil']:
    ss = {k:v for k,v in by.items() if k[0]==mx}
    if not ss: continue
    row = ''
    for m in METRICS:
        rat = []
        for k,d in ss.items():
            if all(t in d for t in (2,5,20,25)):
                first = (d[5][m]-d[2][m])/3
                last  = (d[25][m]-d[20][m])/5
                if first > 0: rat.append(last/first)
        row += f'{st.median(rat):12.2f} ' if rat else '           - '
    print(f'{mx:22s}{row}')
