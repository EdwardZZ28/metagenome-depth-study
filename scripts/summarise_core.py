#!/usr/bin/env python3
"""Collect per-dataset metrics for the core set into one table."""
import csv, glob, os, re, sys

RDM = '/QRISdata/Q6636/depth_study_results'
RESULTS = RDM + '/stageA'
OUT = os.path.expanduser('~/depth_study/core_summary.tsv')
PLAN = os.path.expanduser('~/depth_study/sample_plan.tsv')

HIGH_RISK = re.compile(
    r'\b(NDM|KPC|OXA-48|OXA-23|IMP-|VIM-|MCR-|mcr-|VanA|vanA|VanB|vanB|CTX-M|TEM-|SHV-)', re.I)
ID_MIN = 80.0
COLS = ['sample','matrix','tier','seed','arg_subtypes','copies_16S','contigs','assembly_bp',
        'bins','plasmid','conjugative','virus','rgi_ps','rgi_id80','highrisk',
        'bracken_species','highrisk_genes']

matrix = {}
for r in csv.DictReader(open(PLAN), delimiter='\t'):
    matrix[r['sample']] = r['matrix']

def to_int(x):
    try: return int(x)
    except (TypeError, ValueError): return 0

def read_argoap(p):
    if not os.path.exists(p): return None, None
    L = [l for l in open(p).read().splitlines() if l.strip()]
    tot = 0.0
    for l in L[1:]:
        try: tot += float(l.split('\t')[-1])
        except ValueError: pass
    return len(L) - 1, round(tot, 4)

def read_genomad(d):
    s = glob.glob(f'{d}/genomad/*_summary')
    if not s: return 0, 0, 0
    b = os.path.basename(s[0]).replace('_summary', '')
    pl = cj = vi = 0
    fp = f'{s[0]}/{b}_plasmid_summary.tsv'
    if os.path.exists(fp):
        L = open(fp).read().splitlines(); pl = max(0, len(L) - 1)
        if L:
            h = L[0].split('\t')
            if 'conjugation_genes' in h:
                ci = h.index('conjugation_genes')
                for l in L[1:]:
                    f = l.split('\t')
                    if len(f) > ci and f[ci].strip() not in ('', 'NA'): cj += 1
    fv = f'{s[0]}/{b}_virus_summary.tsv'
    if os.path.exists(fv): vi = max(0, len(open(fv).read().splitlines()) - 1)
    return pl, cj, vi

def read_rgi(d):
    f = glob.glob(f'{d}/rgi/*.txt')
    if not f: return 0, 0, set()
    L = open(f[0]).read().splitlines()
    if not L: return 0, 0, set()
    h = L[0].split('\t')
    need = ('Cut_Off', 'Best_Identities', 'Best_Hit_ARO')
    if not all(c in h for c in need): return 0, 0, set()
    cc, cid, ca = (h.index(c) for c in need)
    n = n80 = 0; genes = set()
    for l in L[1:]:
        f2 = l.split('\t')
        if len(f2) <= max(cc, cid, ca) or f2[cc] not in ('Perfect', 'Strict'): continue
        n += 1
        try: idv = float(f2[cid])
        except ValueError: continue
        if idv >= ID_MIN:
            n80 += 1
            if HIGH_RISK.search(f2[ca]): genes.add(f2[ca])
    return n, n80, genes

def read_bracken(did):
    p = f'{RDM}/kraken/{did}.bracken_S.tsv'
    if not os.path.exists(p): return ''
    return max(0, len(open(p).read().splitlines()) - 1)

rows = []
for d in sorted(glob.glob(f'{RESULTS}/*')):
    did = os.path.basename(d)
    m = re.match(r'(.+)__t(\d+)__s(\d+)$', did)
    if not m: continue
    r = {'sample': m.group(1), 'tier': int(m.group(2)), 'seed': int(m.group(3))}
    r['matrix'] = matrix.get(r['sample'], 'unknown')
    st, cp = read_argoap(f'{d}/argoap/normalized_16S.subtype.txt')
    r['arg_subtypes'] = st if st is not None else ''
    r['copies_16S'] = cp if cp is not None else ''
    size = {}
    if os.path.exists(f'{d}/size.txt'):
        for t in open(f'{d}/size.txt').read().split():
            if '=' in t:
                k, v = t.split('=', 1); size[k] = v
    r['contigs'] = to_int(size.get('contigs'))
    r['assembly_bp'] = to_int(size.get('assembly_bp'))
    r['bins'] = to_int(size.get('bins'))
    r['plasmid'], r['conjugative'], r['virus'] = read_genomad(d)
    r['rgi_ps'], r['rgi_id80'], genes = read_rgi(d)
    r['highrisk'] = len(genes)
    r['highrisk_genes'] = ';'.join(sorted(genes))
    r['bracken_species'] = read_bracken(did)
    rows.append(r)

with open(OUT, 'w') as fh:
    fh.write('\t'.join(COLS) + '\n')
    for r in sorted(rows, key=lambda x: (x['matrix'], x['sample'], x['tier'], x['seed'])):
        fh.write('\t'.join(str(r.get(c, '')) for c in COLS) + '\n')
print(f'{len(rows)} datasets -> {OUT}')
