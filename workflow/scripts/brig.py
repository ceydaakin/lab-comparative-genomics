"""BRIG-method circular comparison (Alikhan et al. 2011): BLASTn of every genome against the focal genome, one ring per genome
(identity >= 90% / 70-90% / < 70% in 1 kb bins), GC skew inside, focal prophages and CRISPR arrays outside.
Also writes focal-only regions (>= 5 kb at < 70% identity in >= 90% of genomes) as candidate genomic islands."""
import argparse, os, subprocess
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from Bio import Phylo, SeqIO

p = argparse.ArgumentParser()
p.add_argument("--focal", required=True); p.add_argument("--genome-dir", required=True); p.add_argument("--genomes", required=True)
p.add_argument("--tree", required=True); p.add_argument("--metadata", required=True); p.add_argument("--workdir", required=True)
p.add_argument("--prophages", default=None); p.add_argument("--crispr-gff", default=None); p.add_argument("--threads", type=int, default=4)
p.add_argument("--figure", required=True, help="output path without extension"); p.add_argument("--islands", required=True)
a = p.parse_args()
os.makedirs(a.workdir, exist_ok=True)
ref = f"{a.genome_dir}/{a.focal}.fna"
seq = "".join(str(r.seq).upper() for r in SeqIO.parse(ref, "fasta")); L = len(seq)
meta = pd.read_csv(a.metadata, sep="\t").set_index("accession")
G = {l.strip() for l in open(a.genomes) if l.strip()}
t = Phylo.read(a.tree, "newick"); t.root_at_midpoint(); t.ladderize()
order = [x.name for x in t.get_terminals() if x.name != a.focal and x.name in G]
subprocess.run(["makeblastdb", "-in", ref, "-dbtype", "nucl", "-out", f"{a.workdir}/ref"], check=True, capture_output=True)
# single-record coordinates: offset multi-contig focal assemblies so positions are global
offs, pos = {}, 0
for r in SeqIO.parse(ref, "fasta"): offs[r.id] = pos; pos += len(r.seq)
BIN = 1000; nb = L // BIN + 1; cov = {}
for g in order:
    o = f"{a.workdir}/{g}.tsv"
    subprocess.run(["blastn", "-query", f"{a.genome_dir}/{g}.fna", "-db", f"{a.workdir}/ref", "-evalue", "1e-10",
                    "-outfmt", "6 sseqid sstart send pident", "-num_threads", str(a.threads), "-out", o], check=True)
    best = np.zeros(nb)
    if os.path.getsize(o):
        b = pd.read_csv(o, sep="\t", header=None, names=["sid", "s", "e", "pid"])
        for sid, s, e, pid in zip(b.sid, b.s, b.e, b.pid):
            lo, hi = sorted((s + offs.get(sid, 0), e + offs.get(sid, 0)))
            sl = slice(lo // BIN, hi // BIN + 1); best[sl] = np.maximum(best[sl], pid)
    cov[g] = best

feats = []
if a.prophages and os.path.exists(a.prophages):
    r = pd.read_csv(a.prophages, sep="\t")
    for x in r[(r.accession == a.focal) & (r.completeness == "intact")].position:
        s, e = (int(v) for v in x.split("-")); feats.append(("intact prophage", s, e, "#6250d6"))
if a.crispr_gff and os.path.exists(a.crispr_gff):
    for line in open(a.crispr_gff):
        c = line.split("\t")
        if len(c) > 8 and c[2] == "repeat_region": feats.append(("CRISPR array", int(c[3]) + offs.get(c[0], 0), int(c[4]) + offs.get(c[0], 0), "#008300"))

win = 10000
skew = np.array([(seq[i:i + win].count("G") - seq[i:i + win].count("C")) / max(1, seq[i:i + win].count("G") + seq[i:i + win].count("C")) for i in range(0, L, win)])
fig = plt.figure(figsize=(9, 9)); fig.patch.set_facecolor("#fcfcfb")
ax = fig.add_subplot(111, projection="polar"); ax.set_facecolor("#fcfcfb"); ax.set_theta_zero_location("N"); ax.set_theta_direction(-1)
theta = np.linspace(0, 2 * np.pi, nb, endpoint=False); w = 2 * np.pi / nb; r0, dr = 3.0, 0.16
for i, g in enumerate(order):
    v = cov[g]
    ax.bar(theta, dr * 0.9, width=w, bottom=r0 + i * dr, color=np.where(v >= 90, "#2a78d6", np.where(v >= 70, "#9ec5f4", "#f0efec")), linewidth=0)
ts = np.linspace(0, 2 * np.pi, len(skew), endpoint=False)
ax.bar(ts, skew * 4, width=2 * np.pi / len(skew), bottom=2.2, color=np.where(skew >= 0, "#1baf7a", "#e87ba4"), linewidth=0)
rmax = r0 + len(order) * dr
for _, s, e, colr in feats:
    a0, a1 = 2 * np.pi * s / L, 2 * np.pi * e / L
    ax.bar((a0 + a1) / 2, 0.35, width=max(a1 - a0, 0.012), bottom=rmax + 0.05, color=colr, linewidth=0)
for mb in range(0, L // 1_000_000 + 1):
    ax.text(2 * np.pi * mb * 1e6 / L, rmax + 0.7, f"{mb} Mb", ha="center", va="center", fontsize=7, color="#52514e")
ax.set_ylim(0, rmax + 0.9); ax.axis("off")
ax.text(0, 0, f"{meta.strain.get(a.focal, a.focal)}\n{L / 1e6:.2f} Mb", ha="center", va="center", fontsize=11, fontweight="bold")
fig.legend(handles=[Patch(color="#2a78d6", label="BLASTn identity >= 90%"), Patch(color="#9ec5f4", label="70-90%"), Patch(color="#f0efec", label="< 70% / absent"),
                    Patch(color="#1baf7a", label="GC skew +"), Patch(color="#e87ba4", label="GC skew -"),
                    Patch(color="#6250d6", label="intact prophage"), Patch(color="#008300", label="CRISPR arrays")], loc="lower center", ncol=4, frameon=False, fontsize=7)
ax.set_title(f"BRIG-method comparison: {len(order)} genomes vs {meta.strain.get(a.focal, a.focal)} (core-tree order, inner to outer)", fontsize=8.5, pad=18)
plt.savefig(a.figure + ".png", dpi=300, bbox_inches="tight"); plt.savefig(a.figure + ".pdf", bbox_inches="tight")

absent = np.mean([cov[g] < 70 for g in order], axis=0) if order else np.zeros(nb)
isl, start = [], None
for i, f in enumerate(absent >= 0.9):
    if f and start is None: start = i
    if (not f or i == nb - 1) and start is not None:
        end = i if not f else i + 1
        if end - start >= 5: isl.append((start * BIN, end * BIN, end - start))
        start = None
pd.DataFrame(isl, columns=["start", "end", "length_kb"]).to_csv(a.islands, sep="\t", index=False)
print(f"{len(isl)} focal-only regions >= 5 kb")
