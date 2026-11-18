"""Merge per-genome RGI reports (<accession>.txt) into one table."""
import argparse, os
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--reports", nargs="+", required=True); ap.add_argument("--out", required=True)
a = ap.parse_args()
rows = []
for f in a.reports:
    t = pd.read_csv(f, sep="\t")
    acc = os.path.basename(f).removesuffix(".txt")
    rows += [{"accession": acc, "cut_off": r["Cut_Off"], "best_hit": r["Best_Hit_ARO"], "identity": r["Best_Identities"],
              "drug_class": r["Drug Class"], "mechanism": r["Resistance Mechanism"], "gene_family": r["AMR Gene Family"]}
             for _, r in t.iterrows()]
pd.DataFrame(rows, columns=["accession", "cut_off", "best_hit", "identity", "drug_class", "mechanism", "gene_family"]).to_csv(a.out, sep="\t", index=False)
print(f"{len(rows)} RGI hits from {len(a.reports)} genomes")
