#!/usr/bin/env python3
"""Draw the three summary figures from data/L1_summary.tsv.

Usage:  python3 scripts/plot_L1.py
Writes: figures/fig1_recovery.png, fig2_complexity.png, fig3_highrisk.png
"""
import csv
import collections
import statistics as st
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
rows = list(csv.DictReader(open(f'{HERE}/data/L1_summary.tsv'), delimiter='\t'))
INT = ['arg_subtypes', 'contigs', 'assembly_bp', 'bins', 'plasmid',
       'conjugative', 'virus', 'rgi_ps', 'rgi_id80', 'highrisk']
for r in rows:
    r['tier'] = int(r['tier'])
    for k in INT:
        r[k] = int(r[k]) if r[k] else 0
    r['copies_16S'] = float(r['copies_16S']) if r['copies_16S'] else 0.0

bysamp = collections.defaultdict(dict)
for r in rows:
    bysamp[r['sample']][r['tier']] = r

TIERS = [2, 5, 10, 15, 20, 25]
METRICS = [('copies_16S', 'ARG abundance (copies/16S)'),
           ('arg_subtypes', 'ARG subtype richness'),
           ('rgi_id80', 'RGI ARG hits (identity >= 80%)'),
           ('contigs', 'Contigs (>= 1 kb)'),
           ('plasmid', 'Plasmid contigs'),
           ('conjugative', 'Conjugative plasmid contigs'),
           ('virus', 'Virus contigs'),
           ('bins', 'MetaBAT2 bins')]


def recovery(metric, tier):
    """Per-sample value at `tier` as a percentage of the same sample at 25 Gbp."""
    out = []
    for s, d in bysamp.items():
        if tier in d and 25 in d and d[25][metric]:
            out.append(100 * d[tier][metric] / d[25][metric])
    return out


plt.rcParams.update({'font.size': 9, 'axes.linewidth': 1.2, 'font.family': 'DejaVu Sans'})
fig, ax = plt.subplots(figsize=(7.2, 4.6))
colors = plt.cm.viridis([i / (len(METRICS) - 1) for i in range(len(METRICS))])
for (m, lab), c in zip(METRICS, colors):
    xs, ys = [], []
    for t in TIERS:
        v = recovery(m, t)
        if v:
            xs.append(t)
            ys.append(st.median(v))
    ax.plot(xs, ys, 'o-', color=c, label=lab, lw=1.8, ms=5, mec='k', mew=0.6)
ax.axhline(100, color='grey', ls=':', lw=1)
ax.set_xlabel('Sequencing depth (Gbp of clean paired reads)')
ax.set_ylabel('Recovery relative to 25 Gbp (%)')
ax.set_title('Depth dependence differs by metric (influent, 11 samples)')
ax.set_xticks(TIERS)
ax.set_ylim(0, 110)
ax.legend(fontsize=7.5, loc='lower right', framealpha=0.95)
ax.tick_params(direction='out', length=4)
fig.tight_layout()
fig.savefig(f'{HERE}/figures/fig1_recovery.png', dpi=200)
plt.close(fig)

fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.9))
xs = [d[25]['contigs'] / 1000 for s, d in bysamp.items() if 10 in d and 25 in d]
for ax, m, lab in [(axes[0], 'arg_subtypes', 'ARG subtype richness'),
                   (axes[1], 'contigs', 'Contigs')]:
    ys = [100 * d[10][m] / d[25][m] for s, d in bysamp.items()
          if 10 in d and 25 in d and d[25][m]]
    ax.scatter(xs, ys, s=46, c='#2c7fb8', edgecolors='k', linewidths=0.7, zorder=3)
    rho, p = spearmanr(xs, ys)
    ax.set_title(f'{lab}\nrho = {rho:+.2f}, p = {p:.3f}', fontsize=9)
    ax.set_xlabel('Community complexity\n(contigs at 25 Gbp, thousands)')
    ax.set_ylabel('Recovery at 10 Gbp (%)')
    ax.set_ylim(0, 100)
    ax.grid(alpha=0.25, zorder=0)
fig.suptitle('What 10 Gbp recovers depends on sample complexity', fontsize=10)
fig.tight_layout()
fig.savefig(f'{HERE}/figures/fig2_complexity.png', dpi=200)
plt.close(fig)

fig, ax = plt.subplots(figsize=(6.4, 3.9))
tot, ns, lab = [], [], []
for t in TIERS:
    v = [r for r in rows if r['tier'] == t and r['contigs'] > 0]
    if not v:
        continue
    tot.append(sum(len([g for g in r['highrisk_genes'].split(';') if g]) for r in v))
    ns.append(sum(1 for r in v if r['highrisk'] > 0))
    lab.append(f'{t}\n(n={len(v)})')
x = range(len(tot))
ax.bar([i - 0.19 for i in x], tot, 0.38, label='gene detections (total)',
       color='#b2182b', ec='k', lw=0.6)
ax.bar([i + 0.19 for i in x], ns, 0.38, label='samples with >= 1 detection',
       color='#fdae61', ec='k', lw=0.6)
ax.set_xticks(list(x))
ax.set_xticklabels(lab)
ax.set_xlabel('Sequencing depth (Gbp)')
ax.set_ylabel('Count')
ax.set_title('Clinically high-risk ARGs are missed at low depth')
ax.legend(fontsize=8)
ax.grid(axis='y', alpha=0.25)
fig.tight_layout()
fig.savefig(f'{HERE}/figures/fig3_highrisk.png', dpi=200)
plt.close(fig)
print('figures written to figures/')
