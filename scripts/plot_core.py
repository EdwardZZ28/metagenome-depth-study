#!/usr/bin/env python3
"""Figures for the core set: recovery by matrix, and subsampling noise.

Usage:  python3 scripts/plot_core.py
Writes: figures/fig4_by_matrix.png, fig5_seed_noise.png
"""
import csv, collections, statistics as st, os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
rows = list(csv.DictReader(open(f'{HERE}/data/core_summary.tsv'), delimiter='\t'))
INT = ['arg_subtypes', 'contigs', 'assembly_bp', 'bins', 'plasmid', 'conjugative',
       'virus', 'rgi_ps', 'rgi_id80', 'highrisk', 'bracken_species']
for r in rows:
    r['tier'] = int(r['tier']); r['seed'] = int(r['seed'])
    for k in INT:
        r[k] = int(r[k]) if r[k] else 0
    r['copies_16S'] = float(r['copies_16S']) if r['copies_16S'] else 0.0

TIERS = [2, 5, 10, 15, 20]
MATRIX_ORDER = ['influent', 'effluent', 'preserved_wastewater', 'river_water', 'soil']
LABELS = {'influent': 'Wastewater influent', 'effluent': 'Wastewater effluent',
          'preserved_wastewater': 'Preserved wastewater', 'river_water': 'River water',
          'soil': 'Soil'}

by = collections.defaultdict(lambda: collections.defaultdict(dict))
for r in rows:
    if r['seed'] == 1:
        by[r['matrix']][r['sample']][r['tier']] = r

def rec(matrix, metric, tier):
    out = []
    for s, d in by[matrix].items():
        if tier in d and 25 in d and d[25][metric]:
            out.append(100 * d[tier][metric] / d[25][metric])
    return out

# ---- Figure 4: four metrics, five matrices ----
plt.rcParams.update({'font.size': 9, 'axes.linewidth': 1.2, 'font.family': 'DejaVu Sans'})
METRICS = [('copies_16S', 'ARG abundance (copies/16S)'),
           ('arg_subtypes', 'ARG subtype richness'),
           ('contigs', 'Contigs (>= 1 kb)'),
           ('plasmid', 'Plasmid contigs')]
colors = {'influent': '#1b4965', 'effluent': '#62b6cb', 'preserved_wastewater': '#5fa8d3',
          'river_water': '#f4a261', 'soil': '#b2182b'}
fig, axes = plt.subplots(2, 2, figsize=(9.2, 7.0))
for ax, (m, lab) in zip(axes.flat, METRICS):
    for mat in MATRIX_ORDER:
        if mat not in by: continue
        xs, ys = [], []
        for t in TIERS:
            v = rec(mat, m, t)
            if v:
                xs.append(t); ys.append(st.median(v))
        n = sum(1 for s, d in by[mat].items() if 25 in d)
        ax.plot(xs, ys, 'o-', color=colors[mat], lw=1.8, ms=5, mec='k', mew=0.6,
                label=f'{LABELS[mat]} (n={n})')
    ax.axhline(100, color='grey', ls=':', lw=1)
    ax.set_title(lab, fontsize=9.5)
    ax.set_xticks(TIERS); ax.set_ylim(0, 112)
    ax.set_xlabel('Depth (Gbp)'); ax.set_ylabel('Recovery vs 25 Gbp (%)')
    ax.grid(alpha=0.2)
axes.flat[0].legend(fontsize=7, loc='lower right', framealpha=0.95)
fig.suptitle('Depth requirement depends on the matrix, not just the metric', fontsize=11)
fig.tight_layout()
fig.savefig(f'{HERE}/figures/fig4_by_matrix.png', dpi=200)
plt.close(fig)

# ---- Figure 5: subsampling noise at 2 and 5 Gbp ----
g = collections.defaultdict(list)
for r in rows:
    if r['tier'] in (2, 5):
        g[(r['sample'], r['tier'], r['matrix'])].append(r)
fig, ax = plt.subplots(figsize=(6.8, 4.0))
mets = ['copies_16S', 'arg_subtypes', 'contigs', 'plasmid']
data = []
for m in mets:
    cv = []
    for k, v in g.items():
        if len(v) >= 3:
            x = [q[m] for q in v]
            if st.mean(x) > 0 and st.median(x) >= 5:
                cv.append(100 * st.stdev(x) / st.mean(x))
    data.append(cv)
bp = ax.boxplot(data, labels=['ARG\nabundance', 'ARG\nrichness', 'Contigs', 'Plasmid\ncontigs'],
                patch_artist=True, widths=0.55)
for p in bp['boxes']:
    p.set_facecolor('#8ecae6'); p.set_edgecolor('k')
for p in bp['medians']:
    p.set_color('k'); p.set_linewidth(1.6)
ax.set_ylabel('Coefficient of variation across 3 seeds (%)')
ax.set_title('Subsampling noise is small at 2 and 5 Gbp')
ax.grid(axis='y', alpha=0.25)
fig.tight_layout()
fig.savefig(f'{HERE}/figures/fig5_seed_noise.png', dpi=200)
plt.close(fig)
print('figures written')
