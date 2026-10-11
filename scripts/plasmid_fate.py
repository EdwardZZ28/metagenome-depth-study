#!/usr/bin/env python3
"""Which plasmid contigs are lost at low depth, and do they differ from those retained?

For each sample, contigs called as plasmids at 25 Gbp are the reference set.
A reference plasmid counts as recovered at tier t if some plasmid contig assembled
at tier t aligns to it over >= MIN_COV of the reference length at >= MIN_ID identity.
"""
import glob, os, re, csv, collections, statistics as st

RDM = '/QRISdata/Q6636/depth_study_results/stageA'
PAF = '/scratch/user/uqnzhai/depth_study/work/plasmid_fate'
TIERS = [2, 5, 10, 15, 20]
MIN_COV = float(os.environ.get('MIN_COV', 0.80))   # fraction of reference covered
MIN_ID = 0.95       # matching bases / alignment block length

samples = [l.strip() for l in open(os.path.expanduser('~/depth_study/l1_samples.txt')) if l.strip()]

def ref_plasmids(sample):
    """name -> (length, has_conjugation_gene) from the 25 Gbp geNomad summary."""
    f = glob.glob(f'{RDM}/{sample}__t25__s1/genomad/*_summary/*_plasmid_summary.tsv')
    if not f:
        return {}
    lines = open(f[0]).read().splitlines()
    if not lines:
        return {}
    h = lines[0].split('\t')
    i_seq = h.index('seq_name'); i_len = h.index('length')
    i_cj = h.index('conjugation_genes') if 'conjugation_genes' in h else -1
    out = {}
    for l in lines[1:]:
        p = l.split('\t')
        if len(p) <= max(i_seq, i_len): continue
        cj = i_cj >= 0 and len(p) > i_cj and p[i_cj].strip() not in ('', 'NA')
        out[p[i_seq]] = (int(p[i_len]), cj)
    return out

def arg_contigs(sample, tier=25):
    """Contig names carrying at least one RGI hit at identity >= 80%."""
    f = glob.glob(f'{RDM}/{sample}__t{tier}__s1/rgi/*.txt')
    if not f:
        return set()
    lines = open(f[0]).read().splitlines()
    if not lines:
        return set()
    h = lines[0].split('\t')
    need = ('Contig', 'Cut_Off', 'Best_Identities')
    if not all(c in h for c in need):
        return set()
    i_c, i_cut, i_id = (h.index(c) for c in need)
    out = set()
    for l in lines[1:]:
        p = l.split('\t')
        if len(p) <= max(i_c, i_cut, i_id) or p[i_cut] not in ('Perfect', 'Strict'):
            continue
        try:
            if float(p[i_id]) < 80: continue
        except ValueError:
            continue
        out.add(p[i_c])   # RGI Contig column already holds the plain contig name
    return out

def recovered(sample, tier, refs):
    """Reference plasmid names recovered at this tier, by coverage of the reference."""
    p = f'{PAF}/{sample}__t{tier}.paf'
    if not os.path.exists(p):
        return None
    iv = collections.defaultdict(list)
    for line in open(p):
        f = line.split('\t')
        if len(f) < 12: continue
        tname, tlen, tstart, tend = f[5], int(f[6]), int(f[7]), int(f[8])
        nmatch, blocklen = int(f[9]), int(f[10])
        if blocklen and nmatch / blocklen >= MIN_ID:
            iv[tname].append((tstart, tend))
    out = set()
    for name, spans in iv.items():
        if name not in refs: continue
        spans.sort(); cov = 0; cur_s, cur_e = spans[0]
        for s, e in spans[1:]:
            if s > cur_e:
                cov += cur_e - cur_s; cur_s, cur_e = s, e
            else:
                cur_e = max(cur_e, e)
        cov += cur_e - cur_s
        if cov / refs[name][0] >= MIN_COV:
            out.add(name)
    return out

LENBINS = [(1000, 2000, '1-2 kb'), (2000, 5000, '2-5 kb'),
           (5000, 10000, '5-10 kb'), (10000, 10**9, '>10 kb')]

rows = []
for s in samples:
    refs = ref_plasmids(s)
    if not refs:
        print(f'{s}: no reference plasmids, skipped'); continue
    args = arg_contigs(s)
    for t in TIERS:
        rec = recovered(s, t, refs)
        if rec is None: continue
        for name, (length, cj) in refs.items():
            rows.append({'sample': s, 'tier': t, 'contig': name, 'length': length,
                         'conjugative': cj, 'has_arg': name in args,
                         'recovered': name in rec})

out = os.path.expanduser('~/depth_study/plasmid_fate.tsv')
with open(out, 'w') as fh:
    fh.write('sample\ttier\tcontig\tlength\tconjugative\thas_arg\trecovered\n')
    for r in rows:
        fh.write(f"{r['sample']}\t{r['tier']}\t{r['contig']}\t{r['length']}\t"
                 f"{int(r['conjugative'])}\t{int(r['has_arg'])}\t{int(r['recovered'])}\n")
print(f'{len(rows)} rows -> {out}')

def rate(sel, t):
    """Median across samples of the per-sample recovery rate."""
    per = collections.defaultdict(lambda: [0, 0])
    for r in rows:
        if r['tier'] != t or not sel(r): continue
        per[r['sample']][0] += r['recovered']; per[r['sample']][1] += 1
    v = [100 * a / b for a, b in per.values() if b >= 20]
    return st.median(v) if v else None

def show(title, groups):
    print(f'\n== {title} ==')
    print(f"{'group':18s}{'n_contigs':>10s}" + ''.join(f'{t:>8d}' for t in TIERS))
    for label, sel in groups:
        n = sum(1 for r in rows if r['tier'] == TIERS[0] and sel(r))
        cells = ''
        for t in TIERS:
            v = rate(sel, t)
            cells += f'{v:7.0f}%' if v is not None else '       -'
        print(f'{label:18s}{n:10d}{cells}')

show('recovery by contig length', [(lab, (lambda r, a=a, b=b: a <= r['length'] < b)) for a, b, lab in LENBINS])
show('recovery by ARG carriage', [('carries ARG', lambda r: r['has_arg']),
                                  ('no ARG', lambda r: not r['has_arg'])])
show('recovery by conjugation genes', [('conjugative', lambda r: r['conjugative']),
                                       ('non-conjugative', lambda r: not r['conjugative'])])
