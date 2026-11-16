"""Pairwise ANI (skani triangle edge list) -> symmetric matrix + single-linkage clonal groups."""
import argparse
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform

p = argparse.ArgumentParser()
p.add_argument("--skani", required=True, help="skani triangle -E output")
p.add_argument("--genomes", required=True)
p.add_argument("--clonal", type=float, required=True)
p.add_argument("--matrix", required=True)
p.add_argument("--groups", required=True)
a = p.parse_args()

G = [l.strip() for l in open(a.genomes) if l.strip()]
e = pd.read_csv(a.skani, sep="\t")
acc = lambda path: path.split("/")[-1].removesuffix(".fna")
e["q"], e["r"] = e["Query_file"].map(acc), e["Ref_file"].map(acc)
m = pd.DataFrame(0.0, index=G, columns=G)
for q, r, v in zip(e.q, e.r, e.ANI):
    if q in m.index and r in m.columns:
        m.loc[q, r] = m.loc[r, q] = max(v, m.loc[q, r])
for g in G: m.loc[g, g] = 100.0
m.round(3).to_csv(a.matrix, sep="\t")
d = squareform((100 - m.values).clip(0), checks=False)
groups = fcluster(linkage(d, "single"), 100 - a.clonal, "distance")
pd.DataFrame({"accession": G, "clonal_group": groups}).to_csv(a.groups, sep="\t", index=False)
off = m.where(m < 100).stack()
print(f"min ANI {off.min():.2f}; {len(set(groups))} clonal groups at >= {a.clonal}%")
