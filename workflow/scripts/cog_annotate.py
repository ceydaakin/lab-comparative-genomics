"""COG categories for the Roary pangenome: rpsblast hits (vs NCBI CDD COG PSSMs) -> COG accession -> COG 2024 category.
Tallies categories per Roary compartment (core/soft-core/shell/cloud) and for focal-strain-specific genes."""
import argparse, collections
import pandas as pd
from Bio import SeqIO

p = argparse.ArgumentParser()
p.add_argument("--rpsblast", required=True); p.add_argument("--reference", required=True)
p.add_argument("--presence", required=True); p.add_argument("--cddid", required=True)
p.add_argument("--cog-def", required=True); p.add_argument("--cog-fun", required=True)
p.add_argument("--focal", required=True); p.add_argument("--genes-out", required=True); p.add_argument("--table-out", required=True)
a = p.parse_args()

h = pd.read_csv(a.rpsblast, sep="\t", header=None, names=["locus", "cdd", "evalue", "bits"])
h = h.sort_values("bits", ascending=False).drop_duplicates("locus")
# Roary reference headers are '>locus_tag group_name'; the presence/absence table is keyed by group_name
loc2gene = {r.id: r.description.split(" ", 1)[1] for r in SeqIO.parse(a.reference, "fasta")}
h["gene"] = h.locus.map(loc2gene)
h["pssm"] = h.cdd.str.split("|").str[-1].astype(int)
cdd = pd.read_csv(a.cddid, sep="\t", header=None, usecols=[0, 1], names=["pssm", "acc"])
cog = pd.read_csv(a.cog_def, sep="\t", header=None, usecols=[0, 1, 2], names=["acc", "cat", "cog_name"], encoding="latin-1")
h = h.merge(cdd, on="pssm", how="left").merge(cog, on="acc", how="left")
fun = {}
for line in open(a.cog_fun, encoding="latin-1"):
    x = line.rstrip("\n").split("\t")
    if len(x) == 4: fun[x[0]] = x[3]

pa = pd.read_csv(a.presence, low_memory=False)
genomes = pa.columns[pa.columns.get_loc("Avg group size nuc") + 1:]
frac = pa[genomes].notna().sum(axis=1) / len(genomes)
pa["category"] = pd.cut(frac, [0, .15, .95, .99, 1], labels=["cloud", "shell", "soft-core", "core"], include_lowest=True).astype(str)
focal = [c for c in genomes if c.startswith(a.focal)][0]
pa["focal_only"] = pa[focal].notna() & (pa[genomes].notna().sum(axis=1) == 1)
g = pa[["Gene", "Annotation", "category", "focal_only"]].merge(h[["gene", "acc", "cat", "cog_name"]], left_on="Gene", right_on="gene", how="left")
g.drop(columns="gene").to_csv(a.genes_out, sep="\t", index=False)

rows = []
for name, sub in [("core", g[g.category == "core"]), ("soft-core", g[g.category == "soft-core"]), ("shell", g[g.category == "shell"]),
                  ("cloud", g[g.category == "cloud"]), ("focal-specific", g[g.focal_only])]:
    c = collections.Counter(x for cats in sub.cat.fillna("-") for x in (cats if cats != "-" else ["-"]))
    for k in sorted(set(fun) | {"-"}):
        rows.append({"set": name, "cog_category": k, "description": fun.get(k, "no COG hit"), "genes": c.get(k, 0),
                     "percent_of_set": round(100 * c.get(k, 0) / len(sub), 2) if len(sub) else 0})
pd.DataFrame(rows).to_csv(a.table_out, sep="\t", index=False)
print("genes with a COG hit:", int(g.cat.notna().sum()), "of", len(g))
