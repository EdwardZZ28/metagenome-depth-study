#!/usr/bin/env python3
"""Do the biological conclusions change with depth, not just their completeness?

Three within-sample checks against the same sample at 25 Gbp:
  1. species profile      (Bracken)       - composition
  2. drug-class profile   (ARGs-OAP type) - composition
  3. ARG-carriage rate among plasmid contigs (geNomad + RGI) - a ratio
"""
import glob, os, re, csv, collections, statistics as st
from scipy.stats import spearmanr

RDM = '/QRISdata/Q6636/depth_study_results'
PLAN = os.path.expanduser('~/depth_study/sample_plan.tsv')
TIERS = [2, 5, 10, 15, 20]
plan = {r['sample']: r['matrix'] for r in csv.DictReader(open(PLAN), delimiter='\t')}
MATRICES = ['influent', 'effluent', 'preserved_wastewater', 'river_water', 'soil']


def header(title):
    print(f'\n=== {title}')
    print(f"{'matrix':22s}" + ''.join(f'{t:>9d}' for t in TIERS))


def show(table, fmt='9.3f'):
    for mx in MATRICES:
        if mx not in table: continue
        row = ''.join(f'{st.median(table[mx][t]):{fmt}}' if table[mx][t] else '        -'
                      for t in TIERS)
        print(f'{mx:22s}{row}')


def profile_check(loader, label, topn):
    by = collections.defaultdict(dict)
    for sample, tier, prof in loader():
        by[sample][tier] = prof
    rho = collections.defaultdict(lambda: collections.defaultdict(list))
    top = collections.defaultdict(lambda: collections.defaultdict(list))
    for s, d in by.items():
        if 25 not in d or not d[25]: continue
        ref = d[25]
        ref_top = [k for k, _ in sorted(ref.items(), key=lambda x: -x[1])[:topn]]
        for t in TIERS:
            if t not in d or not d[t]: continue
            common = set(ref) & set(d[t])
            if len(common) < 8: continue
            mx = plan.get(s, '?')
            rho[mx][t].append(spearmanr([ref[k] for k in common],
                                        [d[t][k] for k in common])[0])
            cur = [k for k, _ in sorted(d[t].items(), key=lambda x: -x[1])[:topn]]
            top[mx][t].append(len(set(cur) & set(ref_top)))
    header(f'Spearman correlation of {label} vs the same sample at 25 Gbp')
    show(rho)
    header(f'how many of the top {topn} {label} at 25 Gbp are also top {topn} here')
    show(top, '9.1f')


def bracken():
    for f in sorted(glob.glob(f'{RDM}/kraken/*.bracken_S.tsv')):
        m = re.match(r'(.+)__t(\d+)__s1\.bracken_S\.tsv$', os.path.basename(f))
        if not m: continue
        d = {}
        L = open(f).read().splitlines()
        if not L: continue
        h = L[0].split('\t')
        if 'name' not in h or 'fraction_total_reads' not in h: continue
        ni, fi = h.index('name'), h.index('fraction_total_reads')
        for l in L[1:]:
            p = l.split('\t')
            if len(p) > max(ni, fi):
                try: d[p[ni]] = float(p[fi])
                except ValueError: pass
        yield m.group(1), int(m.group(2)), d


def drugclass():
    for dd in sorted(glob.glob(f'{RDM}/stageA/*')):
        m = re.match(r'(.+)__t(\d+)__s1$', os.path.basename(dd))
        if not m: continue
        f = f'{dd}/argoap/normalized_16S.type.txt'
        d = {}
        if os.path.exists(f):
            for l in [x for x in open(f).read().splitlines() if x.strip()][1:]:
                p = l.split('\t')
                if len(p) >= 2:
                    try: d[p[0]] = float(p[-1])
                    except ValueError: pass
        yield m.group(1), int(m.group(2)), d


def plasmid_arg_rate():
    out = collections.defaultdict(dict)
    for dd in sorted(glob.glob(f'{RDM}/stageA/*')):
        m = re.match(r'(.+)__t(\d+)__s1$', os.path.basename(dd))
        if not m: continue
        g = glob.glob(f'{dd}/genomad/*_summary')
        if not g: continue
        b = os.path.basename(g[0]).replace('_summary', '')
        fp = f'{g[0]}/{b}_plasmid_summary.tsv'
        if not os.path.exists(fp): continue
        pl = {l.split('\t')[0] for l in open(fp).read().splitlines()[1:] if l.strip()}
        if not pl: continue
        args = set()
        rg = glob.glob(f'{dd}/rgi/*.txt')
        if rg:
            L = open(rg[0]).read().splitlines()
            if L:
                h = L[0].split('\t')
                if all(c in h for c in ('Contig', 'Cut_Off', 'Best_Identities')):
                    ci, cc, cid = (h.index(c) for c in ('Contig', 'Cut_Off', 'Best_Identities'))
                    for l in L[1:]:
                        p = l.split('\t')
                        if len(p) <= max(ci, cc, cid) or p[cc] not in ('Perfect', 'Strict'):
                            continue
                        try:
                            if float(p[cid]) >= 80: args.add(p[ci])
                        except ValueError: pass
        out[m.group(1)][int(m.group(2))] = (len(pl & args), len(pl))
    print('\n=== plasmid contigs carrying an ARG, as % of all plasmid contigs')
    print(f"{'matrix':22s}" + ''.join(f'{t:>9d}' for t in TIERS + [25]))
    bym = collections.defaultdict(list)
    for s, d in out.items():
        bym[plan.get(s, '?')].append(d)
    for mx in MATRICES:
        if mx not in bym: continue
        row = ''
        for t in TIERS + [25]:
            v = [100 * d[t][0] / d[t][1] for d in bym[mx] if t in d and d[t][1] >= 20]
            row += f'{st.median(v):8.2f}%' if v else '        -'
        print(f'{mx:22s}{row}')


profile_check(bracken, 'species abundances', 10)
profile_check(drugclass, 'ARG drug classes', 5)
plasmid_arg_rate()
