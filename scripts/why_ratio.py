#!/usr/bin/env python3
"""Why does the plasmid ARG-carriage ratio rise at low depth?

A. Is the same inflation seen over ALL contigs, not just plasmid contigs?
   If yes, coverage (not geNomad) drives it.
B. Does the inflation survive within a fixed contig-length band?
   If it vanishes, contig length drives it.
C. At 25 Gbp, do ARG-bearing plasmids have higher coverage than the rest?
   If yes, the coverage explanation is supported directly.
"""
import glob, os, re, csv, collections, statistics as st

RDM = '/QRISdata/Q6636/depth_study_results'
PLAN = os.path.expanduser('~/depth_study/sample_plan.tsv')
TIERS = [2, 5, 10, 15, 20, 25]
plan = {r['sample']: r['matrix'] for r in csv.DictReader(open(PLAN), delimiter='\t')}
WATER = ('influent', 'effluent', 'preserved_wastewater')

def arg_contigs(dd):
    out = set()
    rg = glob.glob(f'{dd}/rgi/*.txt')
    if not rg: return out
    L = open(rg[0]).read().splitlines()
    if not L: return out
    h = L[0].split('\t')
    if not all(c in h for c in ('Contig', 'Cut_Off', 'Best_Identities')): return out
    ci, cc, cid = (h.index(c) for c in ('Contig', 'Cut_Off', 'Best_Identities'))
    for l in L[1:]:
        p = l.split('\t')
        if len(p) <= max(ci, cc, cid) or p[cc] not in ('Perfect', 'Strict'): continue
        try:
            if float(p[cid]) >= 80: out.add(p[ci])
        except ValueError: pass
    return out

def plasmid_table(dd):
    """seq_name -> (length, coverage) from the geNomad plasmid summary."""
    g = glob.glob(f'{dd}/genomad/*_summary')
    if not g: return {}
    b = os.path.basename(g[0]).replace('_summary', '')
    fp = f'{g[0]}/{b}_plasmid_summary.tsv'
    if not os.path.exists(fp): return {}
    L = open(fp).read().splitlines()
    if not L: return {}
    h = L[0].split('\t')
    i_s = h.index('seq_name'); i_l = h.index('length')
    out = {}
    for l in L[1:]:
        p = l.split('\t')
        if len(p) <= max(i_s, i_l): continue
        cov = None
        m = re.search(r'multi=([0-9.]+)', p[i_s])     # megahit header carries coverage
        if m: cov = float(m.group(1))
        out[p[i_s]] = (int(p[i_l]), cov)
    return out

def all_contig_lengths(dd, did):
    """contig -> length from the .fai index written next to the assembly."""
    f = f'{dd}/{did}.contigs.fa.fai'
    out = {}
    if not os.path.exists(f): return out
    for l in open(f):
        p = l.split('\t')
        if len(p) >= 2:
            try: out[p[0]] = int(p[1])
            except ValueError: pass
    return out

A = collections.defaultdict(lambda: collections.defaultdict(list))   # all contigs
P = collections.defaultdict(lambda: collections.defaultdict(list))   # plasmid contigs
B = collections.defaultdict(lambda: collections.defaultdict(list))   # plasmids, 2-5 kb only
covs = []

for dd in sorted(glob.glob(f'{RDM}/stageA/*')):
    did = os.path.basename(dd)
    m = re.match(r'(.+)__t(\d+)__s1$', did)
    if not m: continue
    s, t = m.group(1), int(m.group(2))
    mx = plan.get(s, '?')
    args = arg_contigs(dd)
    pl = plasmid_table(dd)
    lens = all_contig_lengths(dd, did)
    if lens:
        A[mx][t].append(100 * len(args & set(lens)) / len(lens))
    if pl:
        P[mx][t].append(100 * len(args & set(pl)) / len(pl))
        band = {k for k, (L_, _) in pl.items() if 2000 <= L_ < 5000}
        if len(band) >= 20:
            B[mx][t].append(100 * len(args & band) / len(band))
    if t == 25 and pl:
        f = f'{dd}/depth.txt'
        if os.path.exists(f):
            L = open(f).read().splitlines()
            if L:
                h = L[0].split('\t')
                if 'contigName' in h and 'totalAvgDepth' in h:
                    ci, di = h.index('contigName'), h.index('totalAvgDepth')
                    for l in L[1:]:
                        p2 = l.split('\t')
                        if len(p2) <= max(ci, di) or p2[ci] not in pl: continue
                        try: covs.append((mx, p2[ci] in args, float(p2[di]), pl[p2[ci]][0]))
                        except ValueError: pass

def show(title, tab, fmt='8.3f'):
    print(f'\n=== {title}')
    print(f"{'matrix':22s}" + ''.join(f'{t:>9d}' for t in TIERS))
    for mx in ['influent', 'effluent', 'preserved_wastewater', 'river_water', 'soil']:
        if mx not in tab: continue
        row = ''.join(f'{st.median(tab[mx][t]):{fmt}}%' if tab[mx][t] else '        -'
                      for t in TIERS)
        print(f'{mx:22s}{row}')

show('A. ARG-bearing contigs as % of ALL contigs', A)
show('B. ARG-bearing plasmids as % of plasmid contigs (2-5 kb band only)', B)
show('   plasmid contigs, all lengths (for reference)', P)

print('\n=== C. mean coverage of plasmid contigs at 25 Gbp (jgi depth.txt)')
print(f"{'matrix':22s}{'n_ARG':>7s}{'cov_ARG':>10s}{'n_other':>9s}{'cov_other':>11s}{'ratio':>8s}")
bym = collections.defaultdict(lambda: ([], []))
for mx, has_arg, cov, L_ in covs:
    bym[mx][0 if has_arg else 1].append(cov)
for mx in ['influent', 'effluent', 'preserved_wastewater', 'river_water', 'soil']:
    if mx not in bym: continue
    a, b = bym[mx]
    if not a or not b: 
        print(f'{mx:22s}{len(a):7d}{"-":>10s}{len(b):9d}{"-":>11s}{"-":>8s}')
        continue
    ma, mb = st.median(a), st.median(b)
    print(f'{mx:22s}{len(a):7d}{ma:10.2f}{len(b):9d}{mb:11.2f}{ma/mb:8.2f}')
