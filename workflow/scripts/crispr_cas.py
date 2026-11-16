"""CRISPR arrays (MinCED GFFs) + Cas typing per genome.
Prokka's UniProt/HAMAP names miss type I-B Cas proteins (they come out 'hypothetical'), so Cas genes are found by an MMseqs2
hit (>= 50% id, >= 70% cov, done upstream) to Cas-named NCBI PGAP proteins; types by signature genes (Makarova et al. 2020)."""
import argparse, os
import pandas as pd

p = argparse.ArgumentParser()
p.add_argument("--minced", nargs="+", required=True, help="MinCED GFF per genome (<acc>.gff)")
p.add_argument("--hits", required=True, help="MMseqs2 m8: query=<acc>|<locus>, target=PGAP protein id")
p.add_argument("--cas-ref", required=True, help="FASTA of Cas-named PGAP proteins ('>id product')")
p.add_argument("--out", required=True); p.add_argument("--genes-out", required=True)
a = p.parse_args()

names = {l[1:].split()[0]: l.strip().split(" ", 1)[1] for l in open(a.cas_ref) if l.startswith(">")}
h = pd.read_csv(a.hits, sep="\t", header=None, names=["query", "target", "pid", "qcov"])
h["genome"] = h["query"].str.split("|").str[0]
h["product"] = h.target.map(names).str.replace("MAG: ", "", regex=False)
h.to_csv(a.genes_out, sep="\t", index=False)


def cas_type(products):
    s = " ".join(products).lower(); t = []
    if "cas3" in s or "cst" in s or "cas8b" in s or "i-b" in s: t.append("I-B")
    if "cas9" in s or "csn2" in s: t.append("II-A")
    if "csm" in s or "cas10" in s or "iii-a" in s: t.append("III-A")
    return " + ".join(t) or "fragment"


rows = []
for gff in a.minced:
    g = os.path.basename(gff).removesuffix(".gff")
    arrays = units = 0
    for line in open(gff):
        c = line.split("\t")
        if len(c) > 8 and c[2] == "repeat_region": arrays += 1
        if len(c) > 8 and c[2] == "repeat_unit": units += 1
    prods = list(h[h.genome == g]["product"])
    rows.append({"accession": g, "arrays": arrays, "spacers": max(units - arrays, 0), "cas_genes": len(prods),
                 "cas_type": "none" if not prods else (cas_type(prods) if len(prods) >= 3 else "fragment")})
pd.DataFrame(rows).to_csv(a.out, sep="\t", index=False)
