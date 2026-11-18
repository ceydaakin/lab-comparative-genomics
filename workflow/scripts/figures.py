"""Publication figures from the pipeline outputs (PNG 300 dpi + vector PDF). Genome order = core-genome tree.
Usage: figures.py --results results --focal ACC --out results/figures [names...]; a figure whose inputs are missing is skipped."""
import argparse, glob, os
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from Bio import Phylo

ap = argparse.ArgumentParser()
ap.add_argument("--results", default="results"); ap.add_argument("--focal", required=True); ap.add_argument("--out", required=True)
ap.add_argument("names", nargs="*")
A = ap.parse_args()
R = A.results; os.makedirs(A.out, exist_ok=True)
plt.rcParams.update({"font.size": 8, "axes.spines.top": False, "axes.spines.right": False,
                     "figure.facecolor": "#fcfcfb", "axes.facecolor": "#fcfcfb", "axes.edgecolor": "#52514e"})
CAT = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]  # validated categorical order
FOCAL = A.focal
G = [l.strip() for l in open(f"{R}/qc/genome_set.txt") if l.strip()]
meta = pd.read_csv(f"{R}/metadata.tsv", sep="\t").set_index("accession")
lab = lambda a: str(meta.strain.get(a, a)) + (" (MAG)" if "_bin" in str(meta.strain.get(a, "")) else "")
niche = meta.niche.reindex(G).fillna("Other")
order_n = [niche[FOCAL]] + [n for n in niche.value_counts().index if n != niche[FOCAL]]
order_n = [order_n[1]] + [order_n[0]] + order_n[2:] if len(order_n) > 1 else order_n   # most common niche blue, focal orange
NCOL = {n: CAT[i % len(CAT)] for i, n in enumerate(order_n)}
NMRK = {n: m for n, m in zip(order_n, ["o", "*", "s", "D", "^", "v", "P", "X"])}
tree = Phylo.read(f"{R}/phylogeny/core_fasttree.nwk", "newick"); tree.root_at_midpoint(); tree.ladderize()
ORDER = [t.name for t in tree.get_terminals() if t.name in G]
cg = pd.read_csv(f"{R}/ani/clonal_groups.tsv", sep="\t").set_index("accession").clonal_group
big = cg.value_counts()
CLONAL = set(cg[cg == big.index[0]].index) if len(big) and big.iloc[0] >= 3 else set()


def save(fig, name):
    fig.savefig(f"{A.out}/{name}.png", dpi=300, bbox_inches="tight"); fig.savefig(f"{A.out}/{name}.pdf", bbox_inches="tight"); plt.close(fig)


def bold(ax, order):
    for tl, a in zip(ax.get_yticklabels(), order):
        if a == FOCAL: tl.set_fontweight("bold")


def fig_busco():
    # bracket access: "duplicated" is also a DataFrame method name
    b = pd.read_csv(f"{R}/qc/busco_table.tsv", sep="\t").set_index("accession").sort_values("complete")
    fig, ax = plt.subplots(figsize=(7, max(3, 0.17 * len(b)))); y = np.arange(len(b))
    ax.barh(y, b["complete"] - b["duplicated"], color=CAT[0], height=0.7, label="complete, single copy")
    ax.barh(y, b["duplicated"], left=b["complete"] - b["duplicated"], color=CAT[2], height=0.7, label="complete, duplicated")
    ax.barh(y, b["fragmented"], left=b["complete"], color=CAT[3], height=0.7, label="fragmented")
    cut = float(open(f"{R}/qc/min_busco.txt").read())
    ax.axvline(cut, color="#e34948", lw=1, ls="--"); ax.text(cut + 0.1, len(b) - 0.5, f"{cut:g}% cut-off", color="#e34948", fontsize=7)
    ax.set_yticks(y); ax.set_yticklabels([lab(a) + ("" if inc else "  [excluded]") for a, inc in zip(b.index, b["included"])], fontsize=6.3); bold(ax, list(b.index))
    ax.set_xlim(max(0, b["complete"].min() - 3), 100.5); ax.set_xlabel("BUSCO completeness (%)"); ax.legend(frameon=False, loc="lower right", fontsize=7)
    ax.set_title("Genome completeness (BUSCO)", loc="left"); save(fig, "figS1_busco")


def fig_tree_ani():
    m = pd.read_csv(f"{R}/ani/ani_matrix.tsv", sep="\t", index_col=0).loc[ORDER, ORDER]
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(14, max(5, 0.24 * len(ORDER))), gridspec_kw={"width_ratios": [1, 1.25]})
    Phylo.draw(tree, axes=a1, do_show=False, label_func=lambda c: "", show_confidence=False)
    ymap = {t.name: i + 1 for i, t in enumerate(tree.get_terminals())}; depth = tree.depths(); xm = max(depth.values()) or 1
    for t in tree.get_terminals():
        n = niche.get(t.name, "Other")
        a1.scatter(depth[t], ymap[t.name], s=60 if t.name == FOCAL else 22, marker=NMRK.get(n, "o"), color=NCOL.get(n, "#b9b8b3"), zorder=3, edgecolor="#fcfcfb", lw=0.8)
        a1.text(depth[t] + xm * 0.03, ymap[t.name], lab(t.name), va="center", fontsize=7, fontweight="bold" if t.name == FOCAL else "normal")
    a1.set_xlim(0, xm * 1.6); a1.set_yticks([]); a1.set_ylabel(""); a1.spines["left"].set_visible(False); a1.set_xlabel("substitutions/site")
    a1.set_title("Roary core-genome alignment, FastTree 2 (GTR), midpoint-rooted", loc="left")
    for n in order_n: a1.scatter([], [], marker=NMRK[n], color=NCOL[n], label=n)
    a1.legend(frameon=False, loc="lower left", fontsize=7)
    lo = np.floor(m.where(m < 100).min().min() * 2) / 2
    im = a2.imshow(m.values, cmap="Blues", vmin=lo, vmax=100, aspect="auto")
    a2.set_xticks([]); a2.set_yticks(range(len(ORDER))); a2.set_yticklabels([lab(t) for t in ORDER], fontsize=6)
    a2.set_title("Pairwise ANI (skani), tree order", loc="left"); plt.colorbar(im, ax=a2, fraction=0.04, label="ANI %")
    save(fig, "fig1_tree_ani")


def fig_pan():
    h = pd.read_csv(f"{R}/pangenome/micropan/heaps.tsv", sep="\t").set_index("parameter").value
    rf = pd.read_csv(f"{R}/pangenome/micropan/rarefaction.tsv", sep="\t").set_index("Genome")
    pa = pd.read_csv(f"{R}/roary/gene_presence_absence.Rtab", sep="\t", index_col=0).values.T > 0
    rng = np.random.default_rng(1)
    core = np.mean([(pa[rng.permutation(len(pa))].cumprod(axis=0) > 0).sum(axis=1) for _ in range(1000)], axis=0)
    s = {l.split("\t")[0]: int(l.rstrip("\n").split("\t")[-1]) for l in open(f"{R}/roary/summary_statistics.txt")}
    x = rf.index.values[1:]; mean = rf.values[1:].mean(axis=1)
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(x, mean, color=CAT[0], lw=2); ax.plot(np.arange(1, len(core) + 1), core, color=CAT[1], lw=2)
    ax.text(x[-1], mean[-1], f"  pan-genome {s['Total genes']:,}", va="center"); ax.text(len(core), core[-1], f"  core {s['Core genes']:,}", va="center")
    ax.set_xlim(1, x[-1] * 1.3); ax.set_ylim(0); ax.grid(axis="y", color="#e5e4e0", lw=0.6)
    ax.set_xlabel("genomes (micropan rarefaction)"); ax.set_ylabel("Roary gene families")
    ax.set_title(f"Heaps' law alpha = {h['alpha']:.3f} ({'open' if h['alpha'] < 1 else 'closed'} pangenome)", loc="left"); save(fig, "fig2_pan_core")


def fig_pcoa():
    p = pd.read_csv(f"{R}/pangenome/micropan/pcoa_jaccard.tsv", sep="\t").set_index("genome")
    fig, ax = plt.subplots(figsize=(6.5, 5))
    for n in order_n:
        sub = p[[niche.get(a) == n for a in p.index]]
        ax.scatter(sub.PC1, sub.PC2, marker=NMRK[n], s=120 if FOCAL in sub.index else 40, color=NCOL[n], label=n, edgecolor="#fcfcfb", lw=1)
    for a in p.index:
        if a not in CLONAL and (niche.get(a) != order_n[0] or a == FOCAL):
            ax.annotate(lab(a), (p.PC1[a], p.PC2[a]), xytext=(5, 3), textcoords="offset points", fontsize=6.5)
    ax.axhline(0, color="#e5e4e0", lw=0.6); ax.axvline(0, color="#e5e4e0", lw=0.6)
    ax.set_xlabel(f"PCo1 ({p.var1.iloc[0]:.1f}%)"); ax.set_ylabel(f"PCo2 ({p.var2.iloc[0]:.1f}%)")
    ax.set_title("Gene content PCoA (prabclus Jaccard on Roary presence/absence)", loc="left"); ax.legend(frameon=False, fontsize=7); save(fig, "fig3_pcoa")


def fig_cog():
    t = pd.read_csv(f"{R}/cog/cog_by_set.tsv", sep="\t")
    w = t.pivot(index="cog_category", columns="set", values="percent_of_set").drop(index=["-"], errors="ignore")
    sets = [s for s in ["core", "soft-core", "shell", "cloud", "focal-specific"] if s in w.columns]
    w = w[sets]; w = w[w.max(axis=1) > 0.5].sort_values("core", ascending=False)
    fig, ax = plt.subplots(figsize=(10, 4.2)); x = np.arange(len(w)); bw = 0.8 / len(sets)
    for i, s in enumerate(sets):
        ax.bar(x + (i - (len(sets) - 1) / 2) * bw, w[s], width=bw * 0.9, color=CAT[i], label=s.replace("focal", lab(FOCAL)))
    ax.set_xticks(x); ax.set_xticklabels(w.index); ax.set_ylabel("% of genes in set"); ax.grid(axis="y", color="#e5e4e0", lw=0.6); ax.set_axisbelow(True)
    d = t.drop_duplicates("cog_category").set_index("cog_category").description
    ax.legend(frameon=False, ncol=len(sets), fontsize=7, loc="upper right"); ax.set_title("COG categories across the pangenome", loc="left")
    fig.text(0.01, -0.12, "  ".join(f"{k}: {d[k]}" for k in w.index), fontsize=5.6, wrap=True); save(fig, "fig4_cog")


def fig_crispr():
    cr = pd.read_csv(f"{R}/crispr/crispr_cas.tsv", sep="\t").set_index("accession").reindex(ORDER)
    types = [t for t in cr.cas_type.value_counts().index if t not in ("none", "fragment")]
    COL = {t: CAT[i % len(CAT)] for i, t in enumerate(types)}; COL.update({"fragment": "#e87ba4", "none": "#b9b8b3"})
    fig, ax = plt.subplots(figsize=(7, max(3, 0.2 * len(ORDER)))); y = np.arange(len(ORDER))[::-1]
    ax.barh(y, cr.spacers, color=[COL[c] for c in cr.cas_type], height=0.7)
    ax.scatter(np.full(len(y), -0.8), y, marker="s", s=18, color=[COL[c] for c in cr.cas_type], clip_on=False); ax.set_xlim(-1.5, None)
    ax.set_yticks(y); ax.set_yticklabels([lab(a) for a in ORDER], fontsize=6.5); bold(ax, ORDER)
    ax.legend(handles=[Patch(color=v, label=("Cas " + k) if k not in ("none", "fragment") else {"none": "no Cas genes", "fragment": "Cas fragment (<3 genes)"}[k]) for k, v in COL.items()],
              frameon=False, loc="lower right", fontsize=7)
    ax.set_xlabel("CRISPR spacers (MinCED); square = Cas type"); ax.grid(axis="x", color="#e5e4e0", lw=0.6); ax.set_axisbelow(True)
    ax.set_title("CRISPR-Cas systems per genome", loc="left"); save(fig, "fig5_crispr_cas")


def fig_cazy():
    fam = pd.read_csv(f"{R}/cazy/cazy_family_counts.tsv", sep="\t", index_col=0).reindex(columns=ORDER, fill_value=0)
    pres = fam > 0; var = fam[(pres.mean(axis=1) > 0) & (pres.mean(axis=1) < 1)]
    if var.empty: return
    fig, ax = plt.subplots(figsize=(max(6, len(var) * 0.28), max(3, 0.2 * len(ORDER))))
    im = ax.imshow(var.T.values, cmap="Blues", aspect="auto", vmin=0, vmax=max(2, var.values.max()))
    ax.set_xticks(range(len(var))); ax.set_xticklabels(var.index, rotation=90, fontsize=6.5)
    ax.set_yticks(range(len(ORDER))); ax.set_yticklabels([lab(a) for a in ORDER], fontsize=6.5); bold(ax, ORDER)
    ax.spines[:].set_visible(False); plt.colorbar(im, ax=ax, fraction=0.025, label="gene copies")
    ax.set_title(f"Variable CAZyme families (dbCAN; {len(var)} of {len(fam)} vary)", loc="left"); save(fig, "fig6_cazymes")


def fig_prophage():
    f = f"{R}/phastest/regions.tsv"
    if not os.path.exists(f): return
    r = pd.read_csv(f, sep="\t"); r = r[r.accession.isin(G)]
    C = {"intact": CAT[0], "questionable": CAT[3], "incomplete": "#c9c8c2"}
    cnt = r.pivot_table(index="accession", columns="completeness", values="region", aggfunc="count", fill_value=0).reindex(ORDER, fill_value=0)
    fig, ax = plt.subplots(figsize=(7, max(3, 0.2 * len(ORDER)))); y = np.arange(len(ORDER))[::-1]; left = np.zeros(len(ORDER))
    for yi, a in zip(y, ORDER):
        if a in CLONAL: ax.axhspan(yi - 0.5, yi + 0.5, color=CAT[0], alpha=0.06, lw=0, zorder=0)
    for c in C:
        v = cnt[c].values if c in cnt else np.zeros(len(ORDER))
        ax.barh(y, v, left=left, color=C[c], height=0.7, edgecolor="#fcfcfb", lw=1); left = left + v
    ax.set_yticks(y); ax.set_yticklabels([lab(a) for a in ORDER], fontsize=6.5); bold(ax, ORDER)
    ax.legend(handles=[Patch(color=C[c], label=c) for c in C], frameon=False, loc="lower right", fontsize=7)
    ax.set_xlabel("prophage regions (PHASTEST)"); ax.grid(axis="x", color="#e5e4e0", lw=0.6); ax.set_axisbelow(True)
    ax.set_title("Prophages per genome" + (" (shaded = largest clonal group)" if CLONAL else ""), loc="left"); save(fig, "fig8_prophages")


def fig_is():
    rows = []
    for a in G:
        f = f"{R}/isescan/{a}.tsv"
        if os.path.exists(f) and os.path.getsize(f):
            t = pd.read_csv(f, sep="\t")
            rows += [(a, fam, n) for fam, n in t.family.value_counts().items()]
    if not rows: return
    m = pd.DataFrame(rows, columns=["acc", "family", "n"]).pivot(index="acc", columns="family", values="n").reindex(ORDER).fillna(0)
    m = m[m.sum().sort_values(ascending=False).index]
    fig, ax = plt.subplots(figsize=(7, max(3, 0.2 * len(ORDER))))
    ax.imshow(m.values, cmap="Blues", aspect="auto")
    for i in range(m.shape[0]):
        for j in range(m.shape[1]):
            v = m.values[i, j]
            if v: ax.text(j, i, int(v), ha="center", va="center", fontsize=5.5, color="#0b0b0b" if v < m.values.max() * 0.6 else "#fcfcfb")
    ax.set_xticks(range(m.shape[1])); ax.set_xticklabels(m.columns, rotation=60, ha="right", fontsize=6.5)
    ax.set_yticks(range(len(ORDER))); ax.set_yticklabels([lab(a) + f"  ({int(m.loc[a].sum())})" for a in ORDER], fontsize=6.3); bold(ax, ORDER)
    ax.spines[:].set_visible(False); ax.set_title("Insertion-sequence families per genome (ISEScan)", loc="left"); save(fig, "fig9_IS_elements")



def fig_prokka():
    rows = []
    for a in ORDER:
        s = dict(l.strip().split(": ") for l in open(f"{R}/prokka/{a}/{a}.txt") if ": " in l)
        rows.append({"acc": a, "CDS": int(s.get("CDS", 0)), "tRNA": int(s.get("tRNA", 0)), "Mb": int(s["bases"]) / 1e6})
    t = pd.DataFrame(rows).set_index("acc")
    fig, axes = plt.subplots(1, 3, figsize=(11, max(3, 0.2 * len(ORDER))), sharey=True, gridspec_kw={"width_ratios": [2, 1, 1]})
    y = np.arange(len(ORDER))[::-1]
    for ax, col, c, xl in zip(axes, ["CDS", "tRNA", "Mb"], [CAT[0], CAT[2], CAT[3]], ["protein-coding genes", "tRNA genes", "genome size (Mb)"]):
        ax.barh(y, t[col], color=[CAT[1] if a == FOCAL else c for a in t.index], height=0.7)
        ax.set_xlabel(xl); ax.grid(axis="x", color="#e5e4e0", lw=0.6); ax.set_axisbelow(True)
    axes[0].set_yticks(y); axes[0].set_yticklabels([lab(a) for a in ORDER], fontsize=6.5); bold(axes[0], ORDER)
    axes[0].set_title(f"Prokka annotation per genome ({lab(FOCAL)} in orange)", loc="left"); save(fig, "fig11_prokka_annotation")


def fig_ani():
    m = pd.read_csv(f"{R}/ani/ani_matrix.tsv", sep="\t", index_col=0).loc[ORDER, ORDER]
    fig, ax = plt.subplots(figsize=(max(5, 0.24 * len(ORDER) + 3), max(4, 0.22 * len(ORDER) + 2)))
    lo = np.floor(m.where(m < 100).min().min() * 2) / 2
    im = ax.imshow(m.values, cmap="Blues", vmin=lo, vmax=100, aspect="auto")
    ax.set_xticks(range(len(ORDER))); ax.set_xticklabels([lab(a) for a in ORDER], rotation=90, fontsize=6)
    ax.set_yticks(range(len(ORDER))); ax.set_yticklabels([lab(a) for a in ORDER], fontsize=6); bold(ax, ORDER)
    idx = [i for i, a in enumerate(ORDER) if a in CLONAL]
    if idx: ax.add_patch(plt.Rectangle((min(idx) - 0.5, min(idx) - 0.5), len(idx), len(idx), fill=False, ec=CAT[1], lw=1.5))
    plt.colorbar(im, ax=ax, fraction=0.04, label="ANI %")
    ax.set_title(f"Pairwise ANI (skani); min {m.where(m < 100).min().min():.1f}%" + ("; box = largest clonal group" if idx else ""), loc="left")
    save(fig, "fig12_ani_heatmap")


def fig_bagel4():
    b = pd.read_csv(f"{R}/bagel4/bagel4_summary.tsv", sep="\t")
    cov = open(f"{R}/bagel4/coverage.txt").read().strip()
    cnt = b.groupby("accession").size().reindex(ORDER, fill_value=0)
    best = b.groupby("accession").best_hit.agg(lambda x: "; ".join(sorted(set(x)))).reindex(ORDER).fillna("")
    fig, ax = plt.subplots(figsize=(8, max(3, 0.2 * len(ORDER)))); y = np.arange(len(ORDER))[::-1]
    ax.barh(y, cnt.values, color=[CAT[1] if a == FOCAL else CAT[0] for a in ORDER], height=0.7)
    for yi, n, h in zip(y, cnt.values, best.values):
        if n: ax.text(n, yi, f"  {h}", va="center", fontsize=6.5)
    ax.set_yticks(y); ax.set_yticklabels([lab(a) for a in ORDER], fontsize=6.5); bold(ax, ORDER)
    ax.set_xlim(0, max(2, cnt.max() * 2.5)); ax.set_xlabel("BAGEL4 areas of interest (label = best core-peptide hit)")
    ax.set_title(f"Bacteriocins and RiPPs (BAGEL4); {cov}", loc="left", fontsize=8); save(fig, "fig7_bagel4")


def fig_rgi():
    d = pd.read_csv(f"{R}/rgi/rgi_all.tsv", sep="\t")
    cnt = d.groupby(["accession", "cut_off"]).size().unstack(fill_value=0).reindex(ORDER, fill_value=0)
    fig, ax = plt.subplots(figsize=(7, max(3, 0.2 * len(ORDER)))); y = np.arange(len(ORDER))[::-1]; left = np.zeros(len(ORDER))
    for c, col in [("Perfect", CAT[7]), ("Strict", CAT[3])]:
        v = cnt[c].values if c in cnt else np.zeros(len(ORDER))
        ax.barh(y, v, left=left, color=col, height=0.7, label=f"{c} hits"); left = left + v
    ax.set_yticks(y); ax.set_yticklabels([lab(a) for a in ORDER], fontsize=6.5); bold(ax, ORDER)
    ax.set_xlim(0, max(3, left.max() + 1)); ax.legend(frameon=False, loc="lower right", fontsize=7); ax.set_xlabel("CARD hits (RGI)")
    perfect = int((d.cut_off == "Perfect").sum()) if len(d) else 0
    ax.set_title(f"Antimicrobial resistance (RGI/CARD): {perfect} Perfect hits", loc="left"); save(fig, "fig14_card_rgi")


def fig_plsdb():
    p = pd.read_csv(f"{R}/plsdb/plsdb_hits.tsv", sep="\t")
    p = p[p.accession.isin(ORDER) & (p.call == "present")]
    fig, ax = plt.subplots(figsize=(max(6, 0.6 * max(1, p.plasmid.nunique()) + 4), max(3, 0.2 * len(ORDER))))
    if len(p):
        p = p.assign(label=p.plasmid + " (" + p.plasmid_host.fillna("?").str.replace("_", " ") + ", " + (p.plasmid_length / 1000).round(1).astype(str) + " kb)")
        m = p.pivot_table(index="accession", columns="label", values="hash_fraction", aggfunc="max").reindex(ORDER)
        im = ax.imshow(m.values, cmap="Blues", vmin=0.5, vmax=1, aspect="auto")
        ax.set_xticks(range(m.shape[1])); ax.set_xticklabels(m.columns, rotation=60, ha="right", fontsize=6.5)
        plt.colorbar(im, ax=ax, fraction=0.03, label="shared mash hashes")
    else:
        ax.text(0.5, 0.5, "no PLSDB plasmid detected", ha="center", transform=ax.transAxes); ax.set_xticks([])
    ax.set_yticks(range(len(ORDER))); ax.set_yticklabels([lab(a) for a in ORDER], fontsize=6.5); bold(ax, ORDER)
    ax.spines[:].set_visible(False); ax.set_title("PLSDB plasmids detected (mash screen)", loc="left"); save(fig, "fig15_plsdb")


ALL = ["busco", "tree_ani", "ani", "prokka", "pan", "pcoa", "cog", "crispr", "cazy", "prophage", "is", "bagel4", "rgi", "plsdb"]
for name in (A.names or ALL):
    globals()[f"fig_{name}"](); print("figure", name)
